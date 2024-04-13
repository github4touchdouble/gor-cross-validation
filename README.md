# Usage of vali.py
A script that automatically creates `n` folds of the provided `.db` file and
validated several types of **GOR**.

```py
vali.py --id run_1 --gor 1 3 4 --folds 4 [--ali train/CB513MultipleAlignments/] --db validation/cb513.db
```
- `--id`: name of validation run
- `--gor`: gor types
- `--folds`: number of folds
- `--ali`: gor V switch (points to a directory with seqondary sequence alignments)
- `--db`: path to db file


# Example 
 
- `--id`: run_example
- `--gor`: 3
- `--folds`: 3
- `--db`: validation/cb513.db

will result in the following file structure
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
├── JARS                        # relevant scripts
│  ├── evalGor.jar
│  ├── predict.jar
│  ├── train.jar
│  └── trainPredict.jar
├── models                      # trained models pre fold
│  ├── gor3_fold1.mod
│  ├── gor3_fold2.mod
│  └── gor3_fold3.mod
├── train                       # training files where we perform crossvalidation
│  ├── cb513.db
│  └── CB513MultipleAlignments  # directory for gor V
├── validation                  # results of crossvalidation
│  ├── cb513.db
│  ├── gor3_fold1_run_example_SUMMARY.txt                   
│  ├── gor3_fold2_run_example_SUMMARY.txt                   
│  ├── gor3_fold3_run_example_SUMMARY.txt                   
│  ├── gor3_fold1_run_example_DETAILED.txt
│  ├── gor3_fold2_run_example_DETAILED.txt
│  ├── gor3_fold3_run_example_DETAILED.txt
│  ├── gor3_run_example_SUMMARY_plot.scores           
│  └── gor3_run_example_SUMMARY_plot.png      # plot containing the mean and std dev over all folds
├── predictions
│  ├── gor3_fold1_run_example.prd             # predictions of each fold
│  ├── gor3_fold2_run_example.prd
│  └── gor3_fold3_run_example.prd
```
# Gor V
