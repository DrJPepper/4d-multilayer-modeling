import random
import pickle
from dataset import *

def hook_fn(module, inp, output):
    module.collapsed = output.var().item() < 1e-4

class _HeteroWrapper(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x_dict):
        out_dict = {}
        for node_type, x in x_dict.items():
            out_dict[node_type] = self.ff[node_type](x)
        return out_dict


class HeteroBatchNorm(_HeteroWrapper):
    def __init__(self, metadata, hidden_channels):
        super().__init__()
        # one BatchNorm per node type
        self.ff = nn.ModuleDict({
            node_type: pygnn.InstanceNorm(hidden_channels)
            for node_type in metadata[0]
        })


class HeteroActivation(_HeteroWrapper):
    def __init__(self, metadata, activation):
        super().__init__()
        self.ff = nn.ModuleDict({
            node_type: activation()
            for node_type in metadata[0]
        })


class HeteroDropout(_HeteroWrapper):
    def __init__(self, metadata, p):
        super().__init__()
        self.ff = nn.ModuleDict({
            node_type: nn.Dropout(p=p)
            for node_type in metadata[0]
        })


class HeteroPad(_HeteroWrapper):
    def __init__(self, metadata, hidden_channels):
        super().__init__()
        self.ff = ModuleDict({
            ntype: LazyLinear(hidden_channels)
            for ntype in metadata[0]
        })

class ScaledNN(torch.nn.Module):
    def __init__(self, gnn, edge_dim, index):
        super().__init__()
        i = index
        self.collapsed = False
        self.nn = torch.nn.Sequential(
            torch.nn.Linear(edge_dim, 64),
            torch.nn.ReLU(),
            torch.nn.LayerNorm(64),
            torch.nn.Linear(64, 128),
            torch.nn.ReLU(),
            torch.nn.Linear(128, (gnn.nnconv_hidden_channels if i==0 else
                                  gnn.layer_counts[i-1]) *
                                  gnn.layer_counts[i])
            )
        self.scale = torch.nn.Parameter(torch.tensor(0.1))  # learnable scalar

    def forward(self, edge_attr):
        return self.scale * self.nn(edge_attr)


class PerGraphNNConv(MessagePassing):
    def __init__(self, in_channels, out_channels, nn, aggr='mean'):
        super().__init__(aggr=aggr)
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.nn = nn
        
        # Mute the built-in MLP of NNConv by passing None
        self.conv = NNConv(in_channels, out_channels, nn=None, aggr=aggr)

    def forward(self, x, edge_index, edge_attr, batch):
        edge_batch = batch[edge_index[0]]

        weights = torch.empty(
            edge_attr.size(0),
            self.out_channels * self.in_channels,
            device=edge_attr.device,
            dtype=edge_attr.dtype
        )

        for g in edge_batch.unique():
            mask = (edge_batch == g)
            weights[mask] = self.nn(edge_attr[mask])

        weights = weights.view(-1, self.out_channels, self.in_channels)

        self._edge_weights = weights

        out = self.propagate(edge_index, x=x)
        return out

    def message(self, x_j, index, ptr, size_i):
        # x_j : [num_edges, in_channels]
        # weight: [num_edges, out_channels, in_channels]

        w = self._edge_weights
        out = torch.bmm(w, x_j.unsqueeze(-1)).squeeze(-1)
        return out

class HeteroGRUCell(nn.Module):
    def __init__(self, node_types, hidden_dim):
        super().__init__()
        
        # Create a separate GRUCell for each node type
        self.grus = nn.ModuleDict()
        for node_type in node_types:
            self.grus[node_type] = nn.GRUCell(
                input_size=hidden_dim, 
                hidden_size=hidden_dim
            )

    def forward(self, x_dict, h_dict):
        """
        x_dict: Dictionary mapping node_type -> tensor of spatial features (from GNN)
        h_dict: Dictionary mapping node_type -> tensor of hidden states (from t-1)
        """
        h_new_dict = {}
        
        for node_type, x in x_dict.items():
            # Get the previous hidden state for this specific node type
            h_prev = h_dict[node_type] if h_dict is not None else None
            
            h_new = self.grus[node_type](x, h_prev)
            
            h_new_dict[node_type] = h_new
            
        return h_new_dict

class MyNNConv(NNConv):
    def forward(self, x, edge_index, edge_attr, batch):
        out = super().forward(x, edge_index, edge_attr)
        return out

class HeteroGNN(nn.Module):

    def __init__(self, metadata, config, example_data):
        super().__init__()
        self.metadata = metadata
        self.dropout_int = config['dropout_int']
        self.dropout_ext = config['dropout_ext']
        self.nnconv_hidden_channels = 16
        self.hidden_channels = config['hidden_channels']
        self.out_channels = config['out_channels']
        self.n_layers = config['n_layers']
        self.do_batch_norm = config['batch_norm']
        self.conv = getattr(pygnn, config['conv'])
        self.convs = []
        self.h_dict_prevs = defaultdict(lambda: None)

        self.layer_counts = self.hidden_channels
        self.network = None

        edge_attr_dict = example_data.edge_attr_dict
        layers = []
        layers.append((HeteroPad(self.metadata, self.nnconv_hidden_channels),
                       'x_dict -> x_dict'))
        for i in range(self.n_layers):
            h = 1
            convs = {}
            ts = [('cell', 'neighbors', 'cell'), ('cell', 'contains', 'edge')]
            nns = {
                    ts[0]: ScaledNN(self, edge_attr_dict[ts[0]].size(-1), i),
                    ts[1]: ScaledNN(self, edge_attr_dict[ts[1]].size(-1), i)
                }
            for v in nns.values():
                v.register_forward_hook(hook_fn)
            nns[('edge', 'in', 'cell')] = nns[('cell', 'contains', 'edge')]
            for edge_type in self.metadata[1]:
                attr_nn = nns[edge_type]
                for m in attr_nn.modules():
                    if isinstance(m, torch.nn.Linear):
                        torch.nn.init.kaiming_uniform_(m.weight, nonlinearity='relu')
                        torch.nn.init.constant_(m.bias, 0.01)  # nonzero bias

                convs[edge_type] = MyNNConv(
                    in_channels=self.nnconv_hidden_channels if i==0 else self.layer_counts[i-1],
                    out_channels=self.layer_counts[i],
                    nn=attr_nn,
                    aggr='mean'
                )

            convs[('edge', 'in', 'cell')].nn.nn[0].register_forward_hook(hook_fn)
            conv_layer = HeteroConv(convs, aggr='mean')
            self.convs.append(convs)
            layers.append((conv_layer, 'x_dict, edge_index_dict, edge_attr_dict, batch -> x_dict'))
            if i != self.n_layers - 1:
                if self.do_batch_norm:
                    layers.append((HeteroBatchNorm(self.metadata, self.layer_counts[i]*h), 'x_dict -> x_dict'))
                layers.append((HeteroActivation(self.metadata, nn.PReLU), 'x_dict -> x_dict'))
                if i > 0:
                    layers.append((HeteroDropout(self.metadata, p=self.dropout_ext), 'x_dict -> x_dict'))
        self.network = pygnn.Sequential('x_dict, edge_index_dict, edge_attr_dict, batch', layers).to(DEVICE)

        hidden_dim = self.hidden_channels[-1]

        # HeteroGRU layer, not currently using
        #self.temporal_gru = HeteroGRUCell(self.metadata[0], hidden_dim).to(DEVICE)
        #self.het_relu = HeteroActivation(self.metadata, nn.PReLU)

        self.query_proj = nn.Linear(hidden_dim, hidden_dim)
        self.key_proj = nn.Linear(hidden_dim, hidden_dim)
        self.value_proj = nn.Linear(hidden_dim, hidden_dim)

        self.mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, self.out_channels),
        )


        self.att_pool_edge = AttentionalAggregation(
            nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, 1),
                nn.Sigmoid()
            )
        )

        self.att_pool_cell = AttentionalAggregation(
            nn.Sequential(
                nn.Linear(hidden_dim+1, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, 1),
                nn.Sigmoid()
            )
        )

        self.global_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim + 1, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, self.out_channels),
        )

    def reset_h_dicts(self):
        h_dict_prevs = defaultdict(lambda: None)

    def embed(self, data):
        edge_batch = None
        cell_batch = None
        try:
            edge_batch = data['edge'].batch
            cell_batch = data['cell'].batch
        except AttributeError:
            edge_batch = torch.tensor([0] * data.x_dict['edge'].shape[0]).to(DEVICE)
            cell_batch = torch.tensor([0] * data.x_dict['cell'].shape[0]).to(DEVICE)
        batch = {'edge': edge_batch, 'cell': cell_batch}
        x_dict = self.network(data.x_dict, data.edge_index_dict, data.edge_attr_dict, batch)
        return x_dict, edge_batch, cell_batch

    def attn(self, embedding, edge_batch, cell_batch):
        edge_x = embedding['edge']
        cell_x = embedding['cell']

        Q_flat = self.query_proj(edge_x)
        K_flat = self.key_proj(cell_x)
        V_flat = self.value_proj(cell_x)
        Q_dense, mask_edge = to_dense_batch(Q_flat, edge_batch) 
        K_dense, mask_cell = to_dense_batch(K_flat, cell_batch)
        V_dense, _         = to_dense_batch(V_flat, cell_batch)

        attn_weights = Q_dense @ K_dense.transpose(1, 2) / (self.hidden_channels[-1] ** 0.5)
        mask_combined = mask_cell.unsqueeze(1) # [B, 1, Max_Cells]
        attn_weights = attn_weights.masked_fill(~mask_combined, float("-1e20"))

        attn = F.softmax(attn_weights, dim=-1)
        context_dense = attn @ V_dense
        context_flat = context_dense[mask_edge] 

        edge_with_context = torch.cat([edge_x, context_flat], dim=-1)
        return edge_with_context

    def edge_mlp(self, context_embedding):
        out = self.mlp(context_embedding)
        return out

    def signal_mlp(self, signals, embedding, edge_batch, cell_batch):
        edge_x = embedding['edge']
        cell_x = embedding['cell']

        signals = signals[:, 0].view(-1, 1) / 10.0
        cell_x = torch.cat([cell_x, signals], dim=1)

        edge_pool = self.att_pool_edge(edge_x, index=edge_batch)
        cell_pool = self.att_pool_cell(cell_x, index=cell_batch)

        graph_embedding = torch.cat([edge_pool, cell_pool], dim=-1)

        signal_out = self.global_mlp(graph_embedding)
        return signal_out

    def forward(self, data, train=False):
        embedding, edge_batch, cell_batch = self.embed(data)
        context_embedding = self.attn(embedding, edge_batch, cell_batch)
        out = self.edge_mlp(context_embedding)
        signal_out = self.signal_mlp(data['cell'].x, embedding, edge_batch, cell_batch)

        if not train:
            out = out[data.train_mask]
            out = F.softmax(out, dim=1)
            out = out.argmax(dim=1)

            signal_out = F.softmax(signal_out, dim=1)
            signal_out = signal_out.argmax(dim=1)

        return out, signal_out

class ExplainerAdapter(torch.nn.Module):
    def __init__(self, original_model):
        super().__init__()
        self.model = original_model

    def forward(self, x_dict, edge_index_dict, edge_attr_dict=None, train_mask=None, *args, **kwargs):
        temp_data = HeteroData()

        for node_type, x in x_dict.items():
            temp_data[node_type].x = x
        for edge_type, edge_index in edge_index_dict.items():
            temp_data[edge_type].edge_index = edge_index

        if edge_attr_dict is not None:
            for edge_type, attr in edge_attr_dict.items():
                temp_data[edge_type].edge_attr = attr

        if train_mask is not None:
                temp_data.train_mask = train_mask

        return self.model(temp_data, train=True)[0]

class ExperimentGenerator(object):

    def __init__(self, dataset, config, rl=False):
        self.data_loader = DataLoader(dataset.data_list, batch_size=config["batch_size"], shuffle=True)
        self.model = HeteroGNN(dataset.data_list[0].metadata(), config, dataset[0]).to(DEVICE)

        self.optimizer = getattr(optim, config["optimizer"])(
            self.model.parameters(), lr=config["learning_rate"])
        self.scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(self.optimizer, T_max=50)
        weights = ft([min(config["max_weight_value"], i) for i in dataset.weights]).to(DEVICE)
        self.dataset = dataset
        config["class_weights"] = weights.tolist()
        self.criterion = nn.CrossEntropyLoss(weight=weights)

        self.model.train()
        self.data_len = len(self.data_loader)
        self.epoch = 0
        self.gamma = 1.0
        self.episodes = self.dataset.steps
        self.config = config
        self.overall_max = 0
        self.best_reward = 0
        self.best_model = None
        self.rl = rl
        self.edges = {}
        self.action_totals = {}
        self.cells = {}
        self.training_gui = True
        self.rest_length_step = self.dataset.cutoff / 1.05
        for i in self.dataset.cells.keys():
            self.cells[i] = 5
            #self.cells[i] = random.choice(range(0,11))
        for i in self.dataset.edges.keys():
            if self.rest_length_step is None:
                self.rest_length_step = self.dataset.edges[i][0]
            self.edges[i] = self.dataset.edges[i][0]

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def init_viz(self, correct=None):
        pld = {'action': 'init', 'scene': []}
        for i in self.dataset.cells.keys():
            c = self.dataset.cells[i]
            pld['scene'].append({'type': 'point', 'position': c.center,
                                 'opacity': 0.5, 'color':
                                 ratio_to_rgb(self.cells[i] / 5.0)})
        if correct:
            for i, c in correct.items():
                color = None
                if len(c) == 2:
                    if sum(c) == 0:
                        color = [1.0, 0.0, 0.0]
                    elif sum(c) == 1:
                        color = [1.0, 1.0, 0.0]
                    else:
                        color = [0.0, 1.0, 0.0]
                else:
                    if c[0]:
                        color = [0.0, 1.0, 0.0]
                    else:
                        color = [1.0, 0.0, 0.0]
                pld['scene'].append({'type': 'vector', 'color': color,
                                     'position':
                                     self.dataset.edge_positions[i]})
        else:
            for i in self.dataset.edge_positions.keys():
                pld['scene'].append({'type': 'vector', 'position':
                                     self.dataset.edge_positions[i]})
        try:
            if self.training_gui:
                requests.post("http://127.0.0.1:8000/update_scene", json=pld)
        except requests.exceptions.ConnectionError:
            if self.training_gui:
                print("Port 8000 closed, disabling calls to GUI server")
                self.training_gui = False

    def next_rl(self, print_status=False):
        
        episode_data = []
        episode_actions = []
        episode_rewards = []
        episode_targets = []

        # Turn off gradients completely for the acting phase
        self.model.eval()

        for i in self.cells.keys():
            self.cells[i] = 5
            #self.cells[i] = random.choice(range(0,11))
        for i in self.dataset.edges.keys():
            self.action_totals[i] = 0
        class_correct = 0
        total = 0
        t = 0
        stopped = False

        with torch.no_grad():
            zero_count = 0
            one_count = 0
            two_count = 0
            while not stopped:
                data = deepcopy(self.dataset.graph_map_e)
                for i in data.keys():
                    signals = []
                    rest_lengths = []
                    yl = []
                    for j in data[i]['cell'].indices:
                        signals.append(self.cells[j.item()] / 10.0)
                    data[i]['cell'].x[:,0] = ft(signals).to(DEVICE)

                    for j in data[i]['edge'].indices:
                        orig = self.edges[j.item()]
                        length = orig +\
                                self.action_totals[j.item()] *\
                                self.rest_length_step
                        length = min(max(length, orig * 0.2), orig * 2.0)
                        rest_lengths.append([length])
                        c = rest_lengths[-1][0]
                        n = self.dataset.targets[j.item()]
                        diff = c - n
                        if abs(diff) < self.dataset.cutoff:
                            yl.append(0)
                        elif diff < 0.0:
                            yl.append(1)
                        else:
                            yl.append(2)
                    data[i]['edge'].x = ((ft(rest_lengths) - self.dataset.data_mean) / self.dataset.data_std).to(DEVICE)
                    data[i]['edge'].y = torch.tensor(yl).type(torch.LongTensor).to(DEVICE)

                data_list = list(data.values())
                cell_inds = [i['cell'].indices[0].item() for i in data_list]
                batch_size = 16000
                data_loader = DataLoader(data_list, batch_size=batch_size, shuffle=False)

                w = 0
                for d in data_loader:
                    # Forward
                    out, signal_out = self.model(d, train=True)
                    
                    # Action selection
                    dist = torch.distributions.Categorical(logits=signal_out)
                    action = dist.sample()
                    
                    # Store raw data
                    episode_data.append(d)
                    episode_actions.append(action)
                    episode_targets.append(d['edge'].y)
                    
                    # Step environment & Get Reward
                    for i in range(len(cell_inds)):
                        s = action[i].item()

                        if s == 1:
                            self.cells[cell_inds[w*batch_size+i]] += 1.0
                        elif s == 2:
                            self.cells[cell_inds[w*batch_size+i]] -= 1.0
                        self.cells[cell_inds[w*batch_size+i]] =\
                            max(min(self.cells[cell_inds[w*batch_size+i]], 10.0), 0.0)
                    w += 1

                    with torch.no_grad():
                        classes = F.softmax(out[d.train_mask], dim=1).argmax(dim=1)
                    zero_count += (classes == 0).sum().item()
                    one_count += (classes == 1).sum().item()
                    two_count += (classes == 2).sum().item()
                    inds = d['edge'].indices[d.train_mask][:,0]
                    corrects = d['edge'].y[d.train_mask]
                    correct_list = (classes == corrects).type(torch.float32)
                    list_to_use = corrects
                    for i in range(len(classes)):
                        if list_to_use[i].item() == 1:
                            self.action_totals[inds[i].item()] += 1 
                        elif list_to_use[i].item() == 2:
                            self.action_totals[inds[i].item()] -= 1 
                    reward = correct_list.count_nonzero().item()
                    class_correct += reward
                    total += len(out[d.train_mask])
                    episode_rewards.append(reward)
                    correct_pairs = defaultdict(lambda: [])
                    for i, k in zip(inds.tolist(), correct_list.type(torch.bool).tolist()):
                        correct_pairs[i].append(k)

                stopped = t > self.episodes
                self.init_viz(correct_pairs)
                t += 1

        print(f'0: {zero_count} | 1: {one_count} | 2: {two_count} | Total: {zero_count+one_count+two_count}')
        self.model.reset_h_dicts()
        self.model.train()

        returns = []
        G = 0
        for r in reversed(episode_rewards):
            G = r + self.gamma * G
            returns.insert(0, G)
        returns = torch.tensor(returns)
        returns = (returns - returns.mean()) / (returns.std() + 1e-8)

        # Create a temporary loader to process the episode in chunks
        train_loader = zip(episode_data, episode_actions, episode_targets, returns)

        sl_total_loss = 0.0
        rl_total_loss = 0.0
        total_loss = 0.0

        for i, (d, action, target, Gt) in enumerate(train_loader):
            out, signal_out = self.model(d, train=True)

            sl_loss = self.criterion(out[d.train_mask], target[d.train_mask])

            dist = torch.distributions.Categorical(logits=signal_out)
            log_prob = dist.log_prob(action)
            entropy = dist.entropy().mean()

            pg_loss = (-log_prob * Gt - 0.05 * entropy).mean()
            total_loss = pg_loss + (sl_loss)

            (total_loss / len(episode_rewards)).backward()

            sl_total_loss += sl_loss.item()
            rl_total_loss += pg_loss.item()
            total_loss += total_loss.item()

        torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
        self.optimizer_rl.step()
        self.scheduler_rl.step()
        self.optimizer_rl.zero_grad()

        corr_percent = class_correct / float(total)
        sl_total_loss /= len(episode_rewards)
        rl_total_loss /= len(episode_rewards)
        total_loss /= len(episode_rewards)
        final_reward = sum(episode_rewards)
        print(f'Correct %: {corr_percent*100:.4f} | Reward: {final_reward} (Best: {self.best_reward}) | SL Loss (Unscaled): {sl_total_loss:.4f}')

        collapsed = True in [i[j].nn.collapsed for i in self.model.convs for j in i.keys()]
        if print_status:
            pass
        self.epoch += 1
        if not collapsed and final_reward > self.best_reward:
            self.best_model = deepcopy(self.model)
            self.best_reward = final_reward

        return sl_loss.item()

    def next(self, print_status=False):
        if self.epoch < self.config['epochs']:
            l_val = 0.0
            class_correct = 0
            total = 0
            if print_status:
                print(f"Epoch: {self.epoch}", end='')
            for data in self.data_loader:
                data_mod = data.to(DEVICE)
                if self.dataset.transform is not None:
                    data_mod = self.dataset.transform(data_mod)
                xs = data_mod['cell'].x[:,0].size()
                data_mod['cell'].x[:,0] = torch.round(torch.rand(xs))
                out = self.model.edge_mlp(
                        self.model.attn(self.model.embed(data_mod)[0],
                                                    data_mod['edge'].batch,
                                                    data_mod['cell'].batch))
                with torch.no_grad():
                    classes = F.softmax(out[data.train_mask], dim=1).argmax(dim=1)
                class_correct += (classes ==
                                  data['edge'].y[data.train_mask]).count_nonzero().item()
                total += len(out[data.train_mask])
                loss = self.criterion(out[data.train_mask], data['edge'].y[data.train_mask])
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                self.optimizer.step()
                self.optimizer.zero_grad()
                l_val += float(loss.item())
            self.model.reset_h_dicts()
            self.scheduler.step()
            corr_percent = class_correct / float(total)
            collapsed = True in [i[j].nn.collapsed for i in self.model.convs for j in i.keys()]
            if print_status:
                print(f" | Loss: {l_val} | Correct %: {corr_percent * 100:.2f}{' (Collapsed)' if collapsed else ''}, {f' (Max %: {self.overall_max})' if self.epoch > 50 else ''}")
            if not collapsed and corr_percent > self.overall_max:
                self.overall_max = corr_percent
                self.best_model = deepcopy(self.model)
            corr_percent = self.overall_max
            self.epoch += 1
            return l_val
        raise StopIteration()


def ratio_to_rgb(ratio: float) -> list:
    min_cap = 0.0
    max_cap = 2.0
    
    result = np.zeros(3)

    if np.isnan(ratio):
        print("WARNING: NaN in ratio_to_rgb", file=sys.stderr)
        return result

    if ratio > max_cap:
        result[:] = [255.0, 0.0, 0.0]
    elif ratio < min_cap:
        result[:] = [0.0, 0.0, 255.0]
    elif ratio > 1.0:
        other = (ratio - 1.0) * 255.0
        green = 255.0 - other
        result[:] = [other, green, 0.0]
    else:
        other = (1.0 - ratio) * 255.0
        green = 255.0 - other
        result[:] = [0.0, green, other]

    return list(result / 255.0)

'''
For curv

CONFIG = {
    # mlop labels
    "architecture": "GNN",
    "dataset": "curv_custom",

    # Model settings
    "dropout_int": 0.0,
    "dropout_ext": 0.0,
    "n_layers": 3,
    "hidden_channels": [94, 64, 22],
    "out_channels": 3,
    "batch_norm": True,

    # Training settings
    "learning_rate": 0.002,
    "epochs": 10000,
    "max_weight_value": 100,
    "batch_size": 8803,
    "optimizer": "RMSprop",
    "conv": "NNConv",
}
'''

CONFIG = {
    # Labels
    "architecture": "GNN",
    "dataset": "curv_custom",

    # Model settings
    "dropout_int": 0.0,
    "dropout_ext": 0.0,
    "n_layers": 2,
    "hidden_channels": [94, 64],
    "out_channels": 3,
    "batch_norm": False,

    # Training settings
    "learning_rate": 0.0004,
    "epochs": 10000,
    "max_weight_value": 3,
    "batch_size": 100,
    "optimizer": "RMSprop",
    "conv": "NNConv",
    "rl": False,
}

def train_model(rl=True):
    config = CONFIG

    ld = LocalDataset('./training_data.json', transform=None)
    ld.build_td()
    runner = ExperimentGenerator(ld, config, rl)
    try:
        for epoch in range(config["epochs"]):
            accuracy = runner.next(True)
    except KeyboardInterrupt:
        print()
    if rl:
        runner.model = runner.best_model
        runner.optimizer_rl = getattr(optim, config["optimizer"])([
                {'params': runner.model.network.parameters(), 'lr': config["learning_rate"]*1e-2},
                {'params': runner.model.query_proj.parameters(), 'lr': config["learning_rate"]},
                {'params': runner.model.key_proj.parameters(), 'lr': config["learning_rate"]},
                {'params': runner.model.value_proj.parameters(), 'lr': config["learning_rate"]},
                {'params': runner.model.mlp.parameters(), 'lr': config["learning_rate"]},
                {'params': runner.model.att_pool_edge.parameters(), 'lr': config["learning_rate"]},
                {'params': runner.model.att_pool_cell.parameters(), 'lr': config["learning_rate"]},
                {'params': runner.model.global_mlp.parameters(), 'lr': config["learning_rate"]}
            ])
        runner.scheduler_rl = torch.optim.lr_scheduler.CosineAnnealingLR(runner.optimizer_rl, T_max=50)
        with torch.no_grad():
            runner.model.global_mlp[0].weight[:, -1].fill_(0.5)
        # Freezing non signal part of model
        for name, param in runner.model.named_parameters():
            check = False
            for k in ["att_pool_", "global_mlp"]:
                check = (k in name) or check
        try:
            while True:
                accuracy = runner.next_rl(True)
        except KeyboardInterrupt:
            pass
    finish_rl(runner, ld)

def finish_rl(runner, ld):
    resp = input('\nSave current model? [y/N] ')
    if resp.lower() not in ['y', 'yes']:
        exit(0)
    torch.save(runner.best_model.state_dict(), 'graph_model_agency.pth')

def finish(runner, ld):
    if runner.best_model is None:
        runner.best_model = runner.model
    runner.best_model.eval()
    data_loader = DataLoader(ld.data_list, batch_size=8000, shuffle=False)
    with torch.no_grad():
        total = 0
        correct = 0
        for data in data_loader:
            out, _ = runner.best_model(data.to(DEVICE))
            target = data["edge"].y[data.train_mask].flatten()
            total += len(out)
            correct += (out == target).count_nonzero()
            print(f'\nCorrect running total: {correct} / {total} ({float(correct) / total})')
        total = 0
        correct = 0
        for i in range(len(list(ld.graph_map.values())[0])):
            for k, v in ld.graph_map.items():
                data = v[i].to(DEVICE)
                out, _ = runner.best_model(data)
                target = data["edge"].y[data.train_mask].flatten()
                total += len(out)
                correct += (out == target).count_nonzero()
        print(f'\nCorrect running total: {correct} / {total} ({float(correct) / total})')
        resp = input('\nSave current model? [y/N] ')
        if resp.lower() not in ['y', 'yes']:
            exit(0)
        torch.save(runner.best_model.state_dict(), 'best_graph_model.pth')
        torch.save(runner.model.state_dict(), 'last_graph_model.pth')

def run_model():
    config = CONFIG
    ld = LocalDataset('./training_data.json')
    dataset, y_list = ld.build_td()
    count = 10
    model = HeteroGNN(dataset[0].metadata(), config, dataset[0]).to(DEVICE)
    model.load_state_dict(torch.load('./graph_model.pth'))

    model.eval()
    correct = 0
    total = 0
    with torch.no_grad():
        for batch in ld.graph_map.values():
            batch = batch[0].to(DEVICE)
            out = model(batch, train=True)
            # assume one graph-level output per graph
            preds = out.argmax(dim=-1)
            correct += (preds == batch['edge'].y).sum().item()
            total += batch['edge'].y.size(0)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'run':
            run_model()
        else:
            train_model()
    else:
        train_model(True)
