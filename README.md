# GOR Crossvalidation
A script that automatically creates `n` folds of the provided `.db` file and
validated several types of **GOR**, by calculating `SOV` and `Q3` in each fold.
The scores get summarized (in `[model]_[run_id]_SUMMARY_plot_ALL.scores`) and plotted in 
(in `[gor_model]_[id]_SUMMARY_plot_ALL.png`).

# Usage
```py
python3 vali.py --id run_1 --gor 1 3 4 --folds 4 [--ali train/CB513MultipleAlignments/] --db cb513.db [--c]
```
- `--id`: name of validation run
- `--gor`: gor types
- `--folds`: number of folds
- `--ali`: gor V switch (points to a directory with secondary sequence alignments) 
           (currently only works with `cb513.db`)
- `--db`: path to db file
- `--c`: clear previous workspace 

# Example 
- `--id`: run_example
- `--gor`: 3
- `--folds`: 3
- `--db`: cb513.db

will result in the following file structure:
```
.
├── crossvalidation
│  ├── test
│  │  ├── fasta_test_1.txt      # fasta files for prediction (3 fold)
│  │  ├── fasta_test_2.txt
│  │  ├── fasta_test_3.txt
│  │  ├── test_1.txt            # validation files (3 fold)
│  │  ├── test_2.txt
│  │  └── test_3.txt
│  └── training
│     ├── train_1.txt           # training files (3 fold)
│     ├── train_2.txt
│     └── train_3.txt
├── models                      # trained models pre fold
│  ├── gor3_run_example_fold1.mod
│  ├── gor3_run_example_fold2.mod
│  └── gor3_run_example_fold3.mod
├── validation                  # results of crossvalidation
│  ├── gor3_fold1_run_example_SUMMARY.txt                   
│  ├── gor3_fold2_run_example_SUMMARY.txt                   
│  ├── gor3_fold3_run_example_SUMMARY.txt                   
│  ├── gor3_fold1_run_example_DETAILED.txt
│  ├── gor3_fold2_run_example_DETAILED.txt
│  ├── gor3_fold3_run_example_DETAILED.txt
│  ├── gor3_run_example_SUMMARY_plot_ALL.scores
│  └── gor3_run_example_SUMMARY_plot_ALL.png  # plot containing the mean and std dev over all folds
├── predictions
│  ├── gor3_fold1_run_example.prd             # predictions of each fold
│  ├── gor3_fold2_run_example.prd
│  └── gor3_fold3_run_example.prd
```
