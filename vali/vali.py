#!/usr/bin/python3

import subprocess
from subprocess import DEVNULL
from tqdm import tqdm
import argparse
import os
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt
import pandas as pd

VALI_DIR='validation/'
TRAIN_DIR='crossvalidation/training/'
TEST_DIR='crossvalidation/test/'
MODEL_DIR="models/"
PREDICTIONS_DIR="predictions/"


def plot_summary_per_fold(path_to_summary: str, name_of_db: str):
    # Load the data
    data = pd.read_csv(path_to_summary, sep='\t', header=None)
    plot_title = path_to_summary.split(".")[0].split("/")[-1].replace("_SUMMARY_plot_ALL", "")
    plot_path = path_to_summary.split(".")[0] + ".png"

    # Calculate mean and standard deviation for each column
    means = data.mean()
    std_devs = data.std()
    labels = ['Q3', 'SOV', 'Q_H', 'Q_E', 'Q_C', 'SOV_H', 'SOV_E', 'SOV_C']

    # Scatter plot of means
    plt.scatter(range(len(means)), means, s=200, c='r', marker='_')

    # Error bars
    plt.errorbar(range(len(means)), means, yerr=std_devs, fmt='none', capsize=5, color='black')

    # X-axis ticks and labels
    plt.xticks(range(len(means)), labels, rotation=45)
    plt.ylim(0,100)
    plt.grid(axis='y',)

    plt.ylabel('Mean + Standard Deviation')
    plt.title("Crossvalidation " + plot_title + " on " + name_of_db)
    plt.tight_layout()
    plt.savefig(plot_path)
    # clear plot for nect plot
    plt.clf()


def create_folds(folds: int, file: str):
    # read file and create a list of tuples (pdb_id, sequence, secondary structure)
    seqs = []
    lines = []

    with open(file, 'r') as f:
        lines = f.readlines()

    for i in range(0, len(lines)):
        if lines[i].startswith('>'):
            pdb_id = lines[i].strip()
            seq = lines[i+1].strip()
            ss = lines[i+2].strip()
            seqs.append((pdb_id, seq, ss))

    kf = KFold(n_splits=folds, shuffle=True, random_state=42)

    for i, (train_index, test_index) in enumerate(kf.split(seqs)):
        train_lines = [seqs[idx] for idx in train_index]
        test_lines = [seqs[idx] for idx in test_index]

        with open(f"{TRAIN_DIR}/train_{i+1}.db", 'w') as f:
            for line in train_lines:
                f.write(line[0] + '\n' + line[1] + '\n' + line[2] + '\n')

        with open(f"{TEST_DIR}/test_{i+1}.db", 'w') as f:
            for line in test_lines:
                f.write(line[0] + '\n' + line[1] + '\n' + line[2] + '\n')

        with open(f"{TEST_DIR}/fasta_test_{i+1}.fasta", 'w') as f:
            for line in test_lines:
                f.write(line[0] + '\n' + line[1] + '\n')


def crossvalidate(folds: int, run_id: str, gor: list, name_of_db: str):
    with tqdm(gor) as t:
        for model in t:
            t.set_description(f'GOR {model}')
            paths_of_summary = []
            for fold in tqdm(range(1, folds + 1), desc='Folds'):
                command = ["java", "-jar", "JARS/trainPredict.jar",
                           "--model",  f"models/gor{model}_fold{fold}_{run_id}.mod",
                           "--seq", f"{TEST_DIR}/fasta_test_{fold}.fasta",
                           "--format", "txt",
                           "--db",  f"{TRAIN_DIR}/train_{fold}.db",
                           "--method", f"gor{model}",
                           "--modelT", f"models/gor{model}_fold{fold}_{run_id}.mod",
                           "--out", f"predictions/gor{model}_fold{fold}_{run_id}.prd"]

                subprocess.run(command, stdout=DEVNULL)

                command = ["java", "-jar", "JARS/evalGor.jar",
                           "-p", f"predictions/gor{model}_fold{fold}_{run_id}.prd",
                           "-r", f"{TEST_DIR}/test_{fold}.db",
                           "-s", f"validation/gor{model}_fold{fold}_{run_id}_SUMMARY.txt",
                           "-d", f"validation/gor{model}_fold{fold}_{run_id}_DETAILED.txt",
                           "-b"]

                subprocess.run(command, stdout=DEVNULL)

                paths_of_summary.append(f"validation/gor{model}_fold{fold}_{run_id}_SUMMARY_plot.scores")

            # concatenate all the summary files per model
            write_summary_to_file(paths_of_summary, name_of_db)



def crossvalidate_gor_v(folds: int, run_id: str, gor: list, name_of_db: str, ali: str):
    with tqdm(gor) as t:
        for model in t:
            t.set_description(f'GOR 5-{model}')
            paths_of_summary = []

            for fold in tqdm(range(1, folds + 1), desc='Folds'):
                command = ["java", "-jar", "JARS/trainPredict.jar",
                           "--model",  f"models/gor5_{model}_fold{fold}_{run_id}.mod",
                           "--maf", ali,
                           "--format", "txt",
                           "--db",  f"{TRAIN_DIR}/train_{fold}.db",
                           "--method", f"gor{model}",
                           "--modelT", f"models/gor5_{model}_fold{fold}_{run_id}.mod",
                           "--out", f"predictions/gor5_{model}_fold{fold}_{run_id}.prd"]

                subprocess.run(command, stdout=DEVNULL)

                command = ["java", "-jar", "JARS/evalGor.jar",
                           "-p", f"predictions/gor5_{model}_fold{fold}_{run_id}.prd",
                           "-r", "cb513.db",
                           "-s", f"validation/gor5_{model}_fold{fold}_{run_id}_SUMMARY.txt",
                           "-d", f"validation/gor5_{model}_fold{fold}_{run_id}_DETAILED.txt",
                           "-b"]

                subprocess.run(command, stdout=DEVNULL)

                paths_of_summary.append(f"validation/gor5_{model}_fold{fold}_{run_id}_SUMMARY_plot.scores")

            # concatenate all the summary files per model
            write_summary_to_file(paths_of_summary, name_of_db)


def write_summary_to_file(paths: list, name_of_db: str):

    summary_out_path = paths[0].split(".")[0].replace("fold1_", "") + "_ALL.scores"
    os.system(f"touch {summary_out_path}")
    with open(summary_out_path, 'w') as merged_file:
        # crate file (summary_out_path)
        # Iterate over each file path
        for file_path in paths:
            try:
                # Open each file and read its content
                with open(file_path, 'r') as file:
                    content = file.read()
                    # Write the content to the merged file
                    merged_file.write(content)
                    # Add a newline after each file's content
                    merged_file.write('\n')
            except FileNotFoundError:
                print(f"File {file_path} not found.")

    plot_summary_per_fold(summary_out_path, name_of_db)


def reset_workspace(clear: bool):
    # check if dirs exist
    if not os.path.exists(MODEL_DIR):
        os.mkdir(MODEL_DIR)
    if not os.path.exists(PREDICTIONS_DIR):
        os.mkdir(PREDICTIONS_DIR)
    if not os.path.exists(VALI_DIR):
        os.mkdir(VALI_DIR)
    if not os.path.exists(TRAIN_DIR):
        os.mkdir(TRAIN_DIR)
    if not os.path.exists(TEST_DIR):
        os.mkdir(TEST_DIR)

    if clear:
        # Clear content of cross-validation folder
        for filename in os.listdir(TRAIN_DIR):
            os.remove(os.path.join(TRAIN_DIR, filename))
        for filename in os.listdir(TEST_DIR):
            os.remove(os.path.join(TEST_DIR, filename))
        for filename in os.listdir(PREDICTIONS_DIR):
            os.remove(os.path.join(PREDICTIONS_DIR, filename))
        for filename in os.listdir(MODEL_DIR):
            os.remove(os.path.join(MODEL_DIR, filename))
        for filename in os.listdir(VALI_DIR):
            os.remove(os.path.join(VALI_DIR, filename))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Crossvalidat on dataset")
    parser.add_argument('--gor', type=int, nargs="+", required=False, help="Gor type to crossvalidate on")
    parser.add_argument('--id', type=str, required=True, help="Run id for the crossvalidation")
    parser.add_argument('--ali', type=str, required=False, help="Points to MultipleAlignments directory, required for gor 5")
    parser.add_argument('--folds', type=int, required=True, default=None, help="Number of folds for crossvalidation")
    parser.add_argument('--db', type=str, required=True, default=None, help="File to split into folds")
    parser.add_argument('--c', action='store_true', required=False, help="Clear previous workspace")
    args = parser.parse_args()
    run_id = args.id
    folds = args.folds
    gor = args.gor
    db = args.db
    ali = args.ali
    name_of_db = db.split(".")[-2].split("/")[-1]

    reset_workspace(args.c)  # reset workspace
    create_folds(folds, db)  # create folds
    crossvalidate(folds, run_id, gor, name_of_db)  # gor I III IV

    if args.ali:
        crossvalidate_gor_v(folds, run_id, gor, name_of_db, ali)  # for gor V
