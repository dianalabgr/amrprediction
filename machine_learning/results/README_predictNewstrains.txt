import xgboost as xgb
import numpy as np
import pickle
file_path = f"/mnt/raid1b/kdan_data/Paper/Machine_Learning/results/amikacin/models/amikacin_random_forest_kmer_3.pkl"
with open(file_path, 'rb') as f:
    model = pickle.load(f)
    
X_npz = np.load("kmer_5_flattened_GCF_000790735.1_AZPAE15041_genomic.fna_.npz", allow_pickle=True)

>>> print(X_npz.files)

X_dna = X_npz['all_kmer_counts_dna']
X_rrna = X_npz['all_kmer_counts_rrna']
X_protein = X_npz['all_kmer_counts_protein']

X_dna = X_dna.reshape(1, -1)
X_rrna = X_rrna.reshape(1, -1)
X_protein = X_protein.reshape(1, -1)

X = np.concatenate((X_dna, X_rrna, X_protein), axis=1)

# Get the number of features the model expects
#num_features = model.num_features()

#print(f"The model expects {num_features} features.")

#print(X.shape)

dtest = xgb.DMatrix(X)
y_pred = model.predict(dtest)
y_pred

X_npz = np.load("kmer_5_flattened_GCF_000794585.1_AZPAE14692_genomic.fna_.npz", allow_pickle=True)

kmer_5_flattened_GCF_001874995.1_ASM187499v1_genomic.fna_.npz
X_npz = np.load("kmer_5_flattened_GCF_001874995.1_ASM187499v1_genomic.fna_.npz", allow_pickle=True)
