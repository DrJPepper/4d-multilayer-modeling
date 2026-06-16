# Tunes the spring constants of the physics simulation using Optuna

import math
import sys
import subprocess
import os
import shutil
import csv
from datetime import datetime

import optuna
from optuna.trial import TrialState

def logging_callback(study, trial):
    try:
        with open(logging_callback.csv_file, "a", newline="") as f:
            writer = csv.writer(f)
            p = trial.params
            writer.writerow([
                trial.number, 
                p["k_linear"],
                p["k_90"],
                p["k_120"],
                p["z_height_mult"],
                p["node_count"],
                trial.value
            ])
    except KeyError:
        pass

def objective(trial):
    fname_prefix = sys.argv[1]
    config = {
            "k_linear": trial.suggest_float("k_linear", 0.01, 25.0),
            "k_ang_90": trial.suggest_float("k_90", 0.01, 25.0),
            "k_ang_120": trial.suggest_float("k_120", 0.01, 25.0),
            "z_height_mult": trial.suggest_float("z_height_mult", 1.0, 4.0),
            "node_count": trial.suggest_int("node_count", 25, 85),
    }

    fname = f'TRIAL_{trial.number}_{fname_prefix}_run.toml'
    shutil.copyfile(f'{fname_prefix}.toml', fname)
    with open(fname, 'a') as f:
        f.writelines([f'{k} = {config[k]}\n' for k in config.keys() if k != "node_count"])
        f.write(f'\n[isolines]\ninitial_node_count = {config["node_count"]}\n')
        f.write(f'spline_bump_count = {math.ceil(config["node_count"]*.12)}\n')

    out_name = f'TRIAL_{trial.number}_out.txt'
    with open(out_name, 'w') as f: 
    	subprocess.run(['../build/agency_sim', '-lf', fname], stdout=f) 
    accuracy = 0.0
    try:
        with open(out_name) as f:
            lines = f.readlines()
            line = lines[-2]
            accuracy = float(line.split(': ')[1])
    except IndexError:
        raise optuna.exceptions.TrialPruned()

    trial.report(accuracy, 0)

    # Handle pruning based on the intermediate value.
    if trial.should_prune():
        raise optuna.exceptions.TrialPruned()

    return accuracy

def main():
    logging_callback.csv_file = f"{sys.argv[1]}_results_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S-%f')}.csv"
    with open(logging_callback.csv_file, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["trial_number", "k_linear", "k_ang_90", "k_ang_120", "z_height_mult", "node_count", "avg_chamfer_distance"])

    study = optuna.create_study(direction="minimize", pruner=optuna.pruners.NopPruner())
    study.optimize(objective, n_trials=1000, timeout=None, n_jobs=1, show_progress_bar=True, callbacks=[logging_callback])

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
