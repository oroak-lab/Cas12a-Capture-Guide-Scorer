# Cas12a-Capture-Guide-Scorer
Python script to score guides for Cas12a Capture experiments.

This guide scoring script was tested with python v3.6.8. Install the following required libraries: numpy, scipy, pandas, sklearn.

# Example Install/Setup
```
python3 -m venv test
source test/bin/activate
pip install --upgrade pip
pip install --upgrade setuptools
pip install numpy
pip install scipy
pip install pandas
pip install sklearn
```

# Example Run
```
$ python --version
Python 3.6.8

$ python cas12a_cap_guide_scorer.py --seq TGTCTTTGGAGAGGCCCAACTGCAAGGTTGACCC --feature 100_selected_feats.txt --training JS_feature_table.csv
[2.79089647]
```
