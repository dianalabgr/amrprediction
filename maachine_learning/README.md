# Guidelines for Machine Learning

1. We must first apply k-mer detection to the data to convert the sequences derived from the annotation pipeline into numerical data. To do this, we run the following command:

```bash
cd Machine_Learning/build_dataset && python build_dataset.py && python kmer_build.py
```

2. For the k-mer dataset, we can run the XGBoost and Random Forest models by running the code in the results folder the bash scripts 

