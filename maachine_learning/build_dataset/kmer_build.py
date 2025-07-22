import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
import numpy as np
import pickle
from tqdm import tqdm




def build_kmer(file, kmer_size):
    X = pickle.load(open("/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/X_" + file +"_chars.pkl", 'rb'))
    headers=pickle.load(open("/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/headers_" + file +".pkl", 'rb'))
	
    print(len(X[0]))

    # find the characters in it and check if they are valid
    chars = set()
    for i in range(len(X)):
        for j in range(len(X[i])):
            for char in X[i][j]:
                chars.add(char)

    print("Characters found in the data:", chars)

    # switch data shape to match their respective part together
    dif_feat = []

    for i in range(len(X[0])):
        dif_feat.append([])
        for j in range(len(X)):
            dif_feat[i].append(X[j][i])
	
	#Find the number of assemblies
    num_assemblies=len(dif_feat[0])
    
    del X

    # create a vectorizer for the kmer
    vectorizer = CountVectorizer(analyzer='char', ngram_range=(kmer_size, kmer_size))

    print("Fitting Vectorizer")
    # find kmers for each feature
    kmer_prot = []
    k_mers=[]

    for i in tqdm(range(len(dif_feat))):
        kmer_prot.append(vectorizer.fit_transform(dif_feat[i]))
        k_mers.append(vectorizer.get_feature_names_out())
    
    # check the distribution buckets of the kmers
    kmer_sizes = [seq[i].shape[1] for seq in kmer_prot]


    print("Converting sparse matrix to list")
    # convert the sparse matrix (of vectorizer) to a list
    new_list = []
    new_list_range=[]
    for i in range(len(kmer_prot)):
        new_list.append(kmer_prot[i].toarray().astype(np.int16))
        new_list_range.append(len(kmer_prot[i].toarray()[0]))
    
    del kmer_prot

    # filter the list to remove the ones with shape 1 (sequence with only X)
    #filtered_list = [elem for elem in new_list if elem.shape[1] != 1]
    df = pd.DataFrame([new_list_range, k_mers], columns=headers)
	# Filter columns where the value is different than 1
    filtered_df = df.loc[:, df.iloc[0] != 1]
	# Get the columns that are retained
    filtered_columns = filtered_df.columns
    # Use these columns to get elements from another list
    filtered_list = [new_list[i] for i, col in enumerate(df.columns) if col in filtered_columns]

    
    print(len(filtered_list))
    print(len(headers))
	
    print("Bringing data to original shape")
    # bring data to original shape
    before_csr_list = [[[] for _ in range(len(filtered_list))] for _ in range(num_assemblies)]
    for i in range(len(filtered_list)):
        for j in range(num_assemblies):
            before_csr_list[j][i] = filtered_list[i][j]
    
    
	# Initialize lists to collect data for the new DataFrame
    indices = []
    kmers = []
    gene_names = []

	# Process the DataFrame to extract k-mers and their gene names
    index_counter = 0
    for gene in filtered_df.columns:
        kmer_list = filtered_df[gene][1]
        for kmer in kmer_list:
            kmers.append(kmer)
            indices.append(index_counter)
            gene_names.append(gene)
            index_counter += 1

	# Create the new DataFrame
    new_df = pd.DataFrame({
        'Index': indices,
        'K-mer': kmers,
        'Gene': gene_names
    })


    # concatenate all the arrays together to avoid padding
    result = []
    print("Concatenating data ")

    for entry in before_csr_list:
        concatenated_vector = np.concatenate(entry)
        result.append(concatenated_vector)
    
    del concatenated_vector
    del before_csr_list

    final_array = np.array(result)
    print("Saving data")
    np.save('/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/kmer_' + str(kmer_size) + '_flattened' + file + '_.npy', final_array)
    np.save('/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/kmer_' + str(kmer_size) + '_flattened' + file + '_stats_.npy', new_df)
    return 


build_kmer("dna", 4)
build_kmer("protein", 4)
build_kmer("rrna", 4)
build_kmer("dna", 5)
build_kmer("protein", 5)
build_kmer("rrna", 5)
