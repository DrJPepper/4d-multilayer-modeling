import math

from model import *

import optuna
from optuna.trial import TrialState

def objective(dataset, trial):
    n_layers = trial.suggest_int("n_layers", 1, 5)
    config = {
        # Labels
        "architecture": "GNN",
        "dataset": "curv_custom",

        # Model settings
        "dropout_int": 0.0,#trial.suggest_float("dropout_int", 0.0, 0.5),
        "dropout_ext": 0.0,#trial.suggest_float("dropout_ext", 0.0, 0.5),
        "n_layers": n_layers,
        "hidden_channels": [trial.suggest_int(f"hidden_units_l{i}", 4, 128) for i in range(n_layers)],
        "out_channels": 3,
        "batch_norm": 0,#bool(trial.suggest_int("batch_norm", 0, 1)),

        # Training settings
        "learning_rate": trial.suggest_float("lr", 1e-5, 1e-1, log=True),
        "epochs": 300,
        "max_weight_value": trial.suggest_float("max_weight", 1, 5),
        "batch_size": trial.suggest_int("batch_size", 128, 15000),
        "optimizer": trial.suggest_categorical("optimizer", ["AdamW", "RMSprop", "SGD"]),
        #"conv": trial.suggest_categorical("conv", ["GCNConv", "GATConv", "GINConv, GatedGraphConv"]),
        "conv": trial.suggest_categorical("conv", ["GATConv"]),
    }

    try:
        runner = ExperimentGenerator(dataset, config)

        accuracy = math.inf
        with tqdm(range(config["epochs"]), leave=False, desc='Exp Progress', postfix='') as p:
            for epoch in p:
                accuracy = min(accuracy, runner.next(False))

                trial.report(accuracy, epoch)
                p.postfix = f'Accuracy: {accuracy}'

                # Handle pruning based on the intermediate value
                if trial.should_prune():
                    raise optuna.exceptions.TrialPruned()
    except torch.OutOfMemoryError:
        raise optuna.exceptions.TrialPruned()

    return accuracy

def main():
    ld = LocalDataset('./training_data.json')
    ld.build_td()

    study = optuna.create_study(direction="minimize")
    try:
        study.optimize(lambda x: objective(ld, x), n_trials=1000, timeout=None, n_jobs=1, show_progress_bar=True)
    except KeyboardInterrupt:
        pass

    pruned_trials = study.get_trials(deepcopy=False, states=[TrialState.PRUNED])
    complete_trials = study.get_trials(deepcopy=False, states=[TrialState.COMPLETE])

    print("Study statistics: ")
    print("  Number of finished trials: ", len(study.trials))
    print("  Number of pruned trials: ", len(pruned_trials))
    print("  Number of complete trials: ", len(complete_trials))

    print("Best trial:")
    trial = study.best_trial

    print("  Value: ", trial.value)

    print("  Params: ")
    for key, value in trial.params.items():
        print("    {}: {}".format(key, value))

if __name__ == '__main__':
    main()
