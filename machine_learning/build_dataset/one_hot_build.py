import numpy as np
import pickle
from scipy.sparse import csr_array
from tqdm import tqdm
from sklearn.preprocessing import OneHotEncoder

# X = pickle.load(open('../data/X_rrna_chars.pkl', 'rb'))

# chars = set()
# for i in range(len(X)):
#     for j in range(len(X[i])):
#         for char in X[i][j]:
#             chars.add(char)

# chars = set()
# for i in range(len(X)):
#     for j in range(len(X[i])):
#         for char in X[i][j]:
#             chars.add(char)

# dif_feat = []

# for i in range(len(X[0])):
#     dif_feat.append([])
#     for j in range(len(X)):
#         dif_feat[i].append(X[j][i])

# sample = dif_feat[1]

# flat_seqs = [item for sequence in sample for item in sequence]

# unique_items = list(set(flat_seqs))

# for i in range(len(sample)):
#     for elem in sample[i]:
#         if elem != 'X':
#             print(i)
#             break

# encoder = OneHotEncoder()
# encoder.fit(np.array(unique_items).reshape(-1, 1))

# final_sample = []
# for seq in sample:
#     seq = [char for char in seq if char in unique_items]
#     final_sample.append(encoder.transform(np.array(seq).reshape(-1, 1)).toarray())

# def encode_sublist(sublist):
#     # Step 1: Flatten the sequences within the sublist to create a single list of unique items
#     flattened_sequences = [item for sequence in sublist for item in sequence]
#     unique_items = list(set(flattened_sequences))
#     if len(unique_items) == 1:
#         return 
    
#     # Step 2: Fit a OneHotEncoder on the unique items
#     encoder = OneHotEncoder(handle_unknown='ignore')
#     encoder.fit(np.array(unique_items).reshape(-1, 1))
    
#     # Step 3: Transform each sequence within the sublist using the fitted encoder
#     def encode_sequence(sequence, encoder):
#         sequence = [char for char in sequence if char in unique_items]
#         encoded = encoder.transform(np.array(sequence).reshape(-1, 1))
#         return encoded.astype(np.int8)
    
#     encoded_sublist = [encode_sequence(sequence, encoder) for sequence in sublist]
#     return encoded_sublist



# encoded_data = []
# for i in tqdm(range(len(dif_feat))):
#     encoded_data.append(encode_sublist(dif_feat[i]))

# enc_sizes = [seq[0].shape[1] for seq in encoded_data if seq is not None]

# enc_dict = {}
# for value in enc_sizes:
#     if value in enc_dict:
#         enc_dict[value] += 1
#     else:
#         enc_dict[value] = 1

# new_list = []
# shape_flag = 0
# for i in range(len(encoded_data)):
#     temp_list = []
#     if(encoded_data[i] is None):
#         continue
#     for j in range(len(encoded_data[i])):
#         if(encoded_data[i][j] is None):
#             continue
#         temp_list.append(encoded_data[i][j].toarray())
#     new_list.append(temp_list)


# before_csr_list = [[[] for _ in range(len(new_list))] for _ in range(16176)]
# for i in range(len(new_list)):
#     for j in range(16176):
#         before_csr_list[j][i] = new_list[i][j]


# new_list = []
# for i in range(len(before_csr_list)):
#     temp_list = []
#     for j in range(len(before_csr_list[i])):
#         temp_list.append(csr_array(before_csr_list[i][j]))
#     new_list.append(temp_list)

# pickle.dump(before_csr_list, open('../data/one_hot_enc/enc_rrna.pkl', 'wb'))


def build_one_hot(X, file_out=None, compress=False):
    X = pickle.load(open(X, 'rb'))

    # find all the characters to be encoded
    chars = set()
    for i in range(len(X)):
        for j in range(len(X[i])):
            for char in X[i][j]:
                chars.add(char)
    print("Characters found in the data:", chars)

    #find all the unique items per number of sequence in the data 
    dif_feat = []

    for i in range(len(X[0])):
        dif_feat.append([])
        for j in range(len(X)):
            dif_feat[i].append(X[j][i])
        def encode_sublist(sublist):
            # Step 1: Flatten the sequences within the sublist to create a single list of unique items
            flattened_sequences = [item for sequence in sublist for item in sequence]
            unique_items = list(set(flattened_sequences))
            if len(unique_items) == 1:
                return 
            
            # Step 2: Fit a OneHotEncoder on the unique items
            encoder = OneHotEncoder(handle_unknown='ignore')
            encoder.fit(np.array(unique_items).reshape(-1, 1))
            
            # Step 3: Transform each sequence within the sublist using the fitted encoder
            def encode_sequence(sequence, encoder):
                sequence = [char for char in sequence if char in unique_items]
                encoded = encoder.transform(np.array(sequence).reshape(-1, 1))
                return encoded.astype(np.int8)
            
            encoded_sublist = [encode_sequence(sequence, encoder) for sequence in sublist]
            return encoded_sublist


    encoded_data = []
    for i in tqdm(range(len(dif_feat))):
        encoded_data.append(encode_sublist(dif_feat[i]))

    new_list = []
    for i in range(len(encoded_data)):
        temp_list = []
        if(encoded_data[i] is None):
            continue
        for j in range(len(encoded_data[i])):
            if(encoded_data[i][j] is None):
                continue
            temp_list.append(encoded_data[i][j].toarray())
        new_list.append(temp_list)

    before_csr_list = [[[] for _ in range(len(new_list))] for _ in range(len(X))]
    for i in range(len(new_list)):
        for j in range(len(X)):
            before_csr_list[j][i] = new_list[i][j]
        
    if compress:
        new_list = []
        for i in range(len(before_csr_list)):
            temp_list = []
            for j in range(len(before_csr_list[i])):
                temp_list.append(csr_array(before_csr_list[i][j]))
            new_list.append(temp_list)
        if file_out:
            pickle.dump(new_list, open(file_out, 'wb'))
        return new_list
    else:
        if file_out:
            pickle.dump(before_csr_list, open(file_out, 'wb'))
        return before_csr_list


build_one_hot('../data/X_rrna_chars.pkl', '../data/one_hot_enc/enc_rrna.pkl', compress=True)
build_one_hot('../data/X_dna_chars.pkl', '../data/one_hot_enc/enc_dna.pkl', compress=True)
build_one_hot('../data/X_protein_chars.pkl', '../data/one_hot_enc/enc_protein.pkl', compress=True)
