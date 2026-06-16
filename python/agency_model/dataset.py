from util import *

def preprocess_edge_attr(edge_attr):
    return edge_attr

class FeatureNoise(T.BaseTransform):
    def __init__(self, sigma=0.01):
        self.sigma = sigma

    def forward(self, data):
        data['edge'].x = data['edge'].x + self.sigma * torch.randn_like(data['edge'].x)
        return data

class Cell():

    def __init__(self, neighbors, edge_map, edges, in_v_order, index, center):
        self.neighbors = neighbors
        self.real_neighbors = [i for i in neighbors if i != -1]
        self.edge_map = edge_map
        self.edges = edges
        self.index = index
        self.in_v_order = in_v_order
        self.on_edge = -1 in self.neighbors
        self.center = center

    def __repr__(self):
        return f'{"Edge " if self.on_edge else ""}Cell #{self.index}'


class LocalDataset(InMemoryDataset):

    def __init__(self, input_file_name, transform=None):
        super().__init__()
        self.td_json = json.load(open(input_file_name))
        self.cells = {}
        self.edges = {}
        self.edge_positions = {}
        self.targets = {}
        self.cell_map = {}
        self.graph_map = defaultdict(lambda: [])
        self.graph_map_e = defaultdict(lambda: [])
        self.edge_map = {}
        self.cell_indices = []
        self.edge_cell_count = 0
        self.process_order = []
        self.epochs = self.td_json['epochs']
        self.steps = len(self.epochs)
        self.transform = transform
        self.cutoff = self.td_json['cutoff']
        self.rest_length_step = self.td_json['rest_length_step']
        for c in self.td_json['cells']:
            self.cells[c['index']] = Cell(c['neighbors'], c['edge_map'],
                                          c['edges'], c['process_order'],
                                          c['index'], c['center'])
            self.edge_cell_count += self.cells[c['index']].on_edge
            self.cell_indices.append(c['index'])
        for e in self.td_json['edges']:
            self.edges[e['index']] = e['lengths']
            self.edge_positions[e['index']] = e['v1'] + e['v2']
            self.targets[e['index']] = e['target']
        self.cell_count = len(self.cells)
        self.weights = [1] * self.num_classes
        self.gen_processing_order()

    @property
    def num_classes(self):
        return 3

    def gen_processing_order(self):
        self.process_order = [set() for _ in range(2)]
        done = set()
        for c in self.cells.values():
            if c.on_edge:
                done.add(c)
                self.process_order[0].add(c)
                for n in c.real_neighbors:
                    nc = self.cells[n]
                    if not nc.on_edge:
                        self.process_order[1].add(nc)
                        done.add(nc)

        while len(done) != self.cell_count:
            self.process_order.append(set())
            for c in self.process_order[-2]:
                for n in c.real_neighbors:
                    nc = self.cells[n]
                    if nc not in done:
                        self.process_order[-1].add(nc)
                        done.add(nc)

    def build_io(self, cell, edge_lengths, e_indices, step):
        if type(cell) is int:
            cell = self.cells[cell]
        data = HeteroData()
        locations = []
        n_inds = []
        n_conns = []
        n_attrs = []

        locations = []
        n_inds = []
        n_conns = []
        n_attrs = []

        locations.append([-1.0] + o7(0))
        n_inds.append(cell.index)
        for i in range(6):
            if cell.neighbors[i] != -1:
                locations.append([-1.0] + o7(i+1))
                n_inds.append(cell.neighbors[i])

        cell_list = [cell] + [self.cells[i] for i in cell.neighbors if i != -1]
        
        for i in range(len(cell_list)):
            for j in range(len(cell_list)):
                c = cell_list[i]
                n = cell_list[j]
                inv_neighbors = {v: k+1 for k, v in zip(range(6), c.neighbors) if v != -1}
                if c != n and ((n.index in cell.neighbors and n.index in c.neighbors) or j == 0):
                    n_conns.append([i, j])
                    n_attrs.append(o7(inv_neighbors[n.index])[1:])

        data['cell'].x = ft(locations)

        edge_map = {}
        edges_to_add = list(set(cell.edges).union(*[self.cells[i].edges for i in cell.neighbors if i >= 0]))
        data.train_mask = torch.tensor([i in cell.edges for i in edges_to_add])
        apical = [0] * len(edges_to_add)
        for i in range(len(edges_to_add)):
            e = edges_to_add[i]
            edge_map[e] = i

        e_conns = []
        e_attrs = []
        cell_list = [cell] + [self.cells[i] for i in cell.neighbors if i != -1]
        for cn in range(len(cell_list)):
            c = cell_list[cn]
            for i in range(len(LAYERS)):
                for j in range(len(DIRS)):
                    l = LAYERS[i]
                    d = DIRS[j]
                    e_conns.append([cn, edge_map[c.edge_map[l][d]]])
                    e_attrs.append([0, 0] + o7(j+1)[1:])
                    e_attrs[-1][i] = 1
            i+=1

        data['cell', 'neighbors', 'cell'].edge_index = torch.tensor(n_conns).transpose(0, 1)
        data['cell', 'neighbors', 'cell'].edge_attr = preprocess_edge_attr(ft(n_attrs))
        data['cell', 'contains', 'edge'].edge_index = torch.tensor(e_conns).transpose(0, 1)
        data['cell', 'contains', 'edge'].edge_attr = preprocess_edge_attr(ft(e_attrs))
        data['edge', 'in', 'cell'].edge_index = data['cell', 'contains', 'edge'].edge_index.flip(0)
        data['edge', 'in', 'cell'].edge_attr = data['cell', 'contains', 'edge'].edge_attr

        data['edge'].indices = torch.tensor(e_indices).type(torch.LongTensor)
        data['cell'].indices = torch.tensor(n_inds).type(torch.LongTensor)
        if cell.index not in self.cell_map.keys():
            self.cell_map[cell.index] = deepcopy(data)

        data['edge'].x = ft(edge_lengths)
        return data

    def build_io_step(self, cell, step):
        e_list_x = []
        e_indices = []
        edges_to_add = list(set(cell.edges).union(*[self.cells[i].edges for i
                                                    in cell.neighbors
                                                    if i >= 0]))
        if cell.index not in self.edge_map.keys():
            self.edge_map[cell.index] = {}
            for i in range(len(edges_to_add)):
                self.edge_map[cell.index][edges_to_add[i]] = i

        for i in range(len(edges_to_add)):
            e = edges_to_add[i]
            e_indices.append([e])
            e_list_x.append([self.edges[e][step]])

        data = self.build_io(cell, e_list_x, e_indices, step)
        yl = []
        for e in edges_to_add:
            c = self.edges[e][step]
            n = self.targets[e]
            diff = c - n
            if abs(diff) < self.cutoff:
                yl.append(0)
            elif diff < 0.0:
                yl.append(1)
            else:
                yl.append(2)
        y = torch.tensor(yl).type(torch.LongTensor)
        data['edge'].y = y
        self.graph_map[cell.index].append(deepcopy(data))
        if not self.graph_map_e[cell.index]:
            data_e = deepcopy(data)
            del data_e['edge'].y
            del data_e['edge'].x
            data_e['cell'].x[:,0] = -1.0
            self.graph_map_e[cell.index] = data_e.to(DEVICE)
        return data, y

    def build_td(self):
        rows = self.cell_count * self.steps
        input_mat = np.zeros((rows, 60)).astype('float64')
        output_mat = np.zeros((rows, 12)).astype('float64')
        edge_map = np.zeros((rows, 12)).astype('int')
        counter = 0
        data_list = []
        y_list = []
        for step in range(self.steps):
            for s in self.process_order:
                for c in s:
                    data, y = self.build_io_step(c, step)
                    data_list.append(data.to(DEVICE))
                    y_list.append(y)
                    counter += 1

        out = torch.cat([d['edge'].x.flatten() for d in data_list])
        self.data_mean = out.mean().item()
        self.data_std = out.std().item()
        self.data_std = self.data_std if self.data_std else 1.0
        for i in range(len(data_list)):
            data_list[i]['edge'].x = (data_list[i]['edge'].x - self.data_mean) / self.data_std
        for i in self.graph_map.keys():
            for j in range(len(self.graph_map[i])):
                self.graph_map[i][j]['edge'].x = (self.graph_map[i][j]['edge'].x - self.data_mean) / (self.data_std)
        self.data, self.slices = self.collate(data_list)

        dl_test = DataLoader(data_list, batch_size=len(data_list))
        balance = [0] * self.num_classes
        for d in dl_test:
            dl = d['edge'].y[d.train_mask].tolist()
            for i in range(self.num_classes):
                balance[i] += dl.count(i)
        m = float(max(balance))
        for i in range(self.num_classes):
            if balance[i]:
                self.weights[i] = m / balance[i]
        print(self.weights)

        self.data_list = data_list
        return data_list, y_list

def main():
    from torch_geometric.explain import Explanation
    from torch_geometric.utils import to_networkx
    import networkx as nx
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from networkx.drawing.nx_agraph import graphviz_layout

    ld = LocalDataset('./training_data.json', transform=None)
    ld.build_td()

    data = ld.data_list[150]
    max_dim = max([data[nt].x.shape[1] for nt in data.node_types if hasattr(data[nt], 'x')])

    # Pad each node type's features to match max_dim
    for nt in data.node_types:
        if hasattr(data[nt], 'x'):
            feature = data[nt].x
            curr_dim = feature.shape[1]
            if curr_dim < max_dim:
                # Pad the last dimension: (pad_left, pad_right)
                pad_amount = max_dim - curr_dim
                data[nt].x = F.pad(feature, (0, pad_amount), "constant", 0)

    homo_data = data.to_homogeneous(node_attrs=['x'], edge_attrs=[])
    G = to_networkx(homo_data, node_attrs=['node_type'], edge_attrs=['edge_type'])

    node_labels = {}
    node_colors = []

    DIRS_SHORT = ['C', 'NE', 'SE', 'S', 'SW', 'NW', 'N']
    c_count = 0
    for node_id, attrs in G.nodes(data=True):
        if attrs['node_type'] == 0:
            ind = data['cell'].x[node_id].tolist()[1:].index(1.0)
            node_labels[node_id] = f'{DIRS_SHORT[ind]}'
            c_count += 1
            node_colors.append('red')
        else:
            nid = node_id-c_count
            x = data["edge"].x[nid,0].item() * ld.data_std + ld.data_mean
            node_labels[node_id] = f'{x:.2f}'
            if data.train_mask[nid]:
                y = ['N', 'E', 'S'][data["edge"].y[nid].item()]
                node_labels[node_id] = f'{node_labels[node_id]},{y}'
                node_colors.append('turquoise')
            else:
                node_colors.append('blue')

    edge_labels = {}
    count = 0

    for u, v, attrs in G.edges(data=True):
        if attrs['edge_type'] == 1:
            ee = data['cell', 'contains', 'edge']['edge_attr'][count].tolist()
            count += 1
            ab = ['A', 'B'][ee[:2].index(1.0)]
            dd = DIRS_SHORT[ee[2:].index(1.0)+1]
            edge_labels[(u, v)] = f'{ab},{dd}'

    highlight_color = 'gold'
    default_color = 'black'
    target_edge_type_index = 0

    edge_colors = []
    for u, v, attrs in G.edges(data=True):
        if attrs['edge_type'] == target_edge_type_index:
            edge_colors.append(highlight_color)
        else:
            edge_colors.append(default_color)

    plt.figure(figsize=(10, 6))
    pos = graphviz_layout(G, prog='neato')
    nx.draw_networkx_nodes(G, pos, node_size=2000, node_color=node_colors)
    nx.draw_networkx_edges(G, pos, width=1.5, edge_color=edge_colors)
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_color="white")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red', 
                                 bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    plt.axis('off')

    red_patch = mpatches.Patch(color='red', label='Cell')
    t_patch = mpatches.Patch(color='turquoise', label='Central Edge')
    blue_patch = mpatches.Patch(color='blue', label='Outer Edge')
    plt.legend(handles=[red_patch, t_patch, blue_patch])

    plt.show()

if __name__ == '__main__':
    main()
