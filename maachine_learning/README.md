# Guidelines for Machine Learning

### *k*-mer
1. We must first apply k-mer algoithm to the data to convert the sequences derived from the annotation pipeline into numerical data. To do this, we run the following command:

```bash
cd Machine_Learning/build_dataset && python build_dataset.py && python kmer_build.py
```

2. For the k-mer dataset, we can run the XGBoost and Random Forest models by running the code in the files `kmer_XGBoost.ipynb` and `kmer_random_forest.ipynb` respectively.

### One-hot Encoding
1. For the one-hot encoded dataset, we must first create the dataset by running the following command:

```bash
cd Machine_Learning && python one_hot_build.py
```

### The folder results_corrected contains some models for antibiotics szh as gentamycin, aztreonam, cefotexaxime with low accuarcy in the evaluation dataset, and i try to remove the species that they are not good with the specific antibiotics.
