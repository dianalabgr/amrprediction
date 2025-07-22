# Description: XGBoost Forest Classifier for Antibiotic Resistance Prediction
# 
# To run the script, use the following command:
# python kmer_rf_xgboost.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/gentamycin/gentamicin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/assemblies.txt --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt

# second antibiogram path: /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/amoxicillin-clavulanicAcid

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
import numpy as np
import pickle
import xgboost as xgb
from sklearn.ensemble import RandomForestClassifier
from scipy.sparse import csr_matrix
from rich.progress import track
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import matplotlib as mpl
from enum import Enum
import shap
import argparse
from sklearn.model_selection import cross_val_score


print(os.getcwd())
# Save the model. Create the model folder if it does not exist
os.makedirs('./models', exist_ok=True)


# Save the results of the best model to a SVG file
def plot_metrics(y_test, y_pred, classes = ['S', 'R'], filename='metrics.svg', cmap='crest'):
    gs = plt.GridSpec(1, 2)

    cm = confusion_matrix(y_test, y_pred)
    # Change labels from 0, 1 to S and R
    cm_df = pd.DataFrame(cm, index=classes, columns=classes)

    fig = plt.figure(figsize=(12, 4))

    ax0 = fig.add_subplot(gs[0, 0])
    _ = ax0.annotate('A', xy=(-0.05, 1.02), xycoords="axes fraction", fontsize=16, fontweight='bold')
    _ = sns.heatmap(cm_df, annot=True, fmt='g', cmap=cmap, linewidths=0.5, annot_kws={"size": 14})  # Increase annotation size here
    _ = plt.title('Confusion Matrix', fontsize=14)  # Title font size
    _ = plt.xlabel('Predicted labels', fontsize=12)  # X-axis label font size
    _ = plt.ylabel('True labels', fontsize=12)  # Y-axis label font size
    _ = plt.xticks(fontsize=10)  # X-tick labels font size
    _ = plt.yticks(fontsize=10)  # Y-tick labels font size

    ax1 = fig.add_subplot(gs[0, 1])
    _ = ax1.annotate('B', xy=(-0.05, 1.02), xycoords="axes fraction", fontsize=16, fontweight='bold')

    clf_report = classification_report(y_test, y_pred, output_dict=True)
    keys_to_plot = [key for key in clf_report.keys() if key not in ('accuracy', 'macro avg', 'weighted avg')]
    df = pd.DataFrame(clf_report, columns=keys_to_plot).T

    rows, cols = df.shape
    mask = np.zeros(df.shape)
    mask[:,cols-1] = True

    ax1 = sns.heatmap(df, mask=mask, annot=True, cmap=cmap, fmt='.4g',
                      vmin=0.0, vmax=1.0, linewidths=2, linecolor='white',
                      annot_kws={"size": 12},  # Increase annotation size here
                      yticklabels=classes)

    # Then, let's add the support column by normalizing the colors in this column
    mask = np.zeros(df.shape)
    mask[:,:cols-1] = True    

    ax1 = sns.heatmap(df, mask=mask, annot=True, cmap=cmap, cbar=False,
                      linewidths=2, linecolor='white', fmt='.0f',
                      vmin=df['support'].min(), vmax=df['support'].sum(),         
                      norm=mpl.colors.Normalize(vmin=df['support'].min(),
                                                vmax=df['support'].sum()),
                      annot_kws={"size": 12},  # Increase annotation size here
                      yticklabels=classes)
            
    _ = plt.title("Classification Report", fontsize=14)  # Title font size
    _ = plt.xticks(rotation=0, fontsize=10)  # X-tick labels font size
    _ = plt.yticks(rotation=360, fontsize=10)  # Y-tick labels font size

    # Save as SVG
    plt.savefig(filename, dpi=300)



# Handle command line arguments
parser = argparse.ArgumentParser(description='XGBoost Classifier for Antibiotic Resistance Prediction')
parser.add_argument('--antibiogram', type=str, help='Path to the antibiogram file')
parser.add_argument('--cores', type=int, help='Number of cores to use for parallel processing')
parser.add_argument('--assemblies', type=str, help='Path to the assemblies txt file')
parser.add_argument('--output', type=str, help='Path to the output txt file')
parser.add_argument('--output_final', type=str, help='Path to the final output txt file')
args = parser.parse_args()

# Important Variables
njobs = args.cores
random_state = 42

class Algorithm(Enum):
    RANDOM_FOREST = 'random_forest'
    XGBOOST = 'xgboost'


# -------------------------
# Load the data
# -------------------------
antibiogram = pd.read_csv(args.antibiogram)
antibiogram['phenotype'].value_counts()

print(f"Preparing data for antibiotic: {antibiogram['antibiotic'].iloc[0]}")
print("-" * 40)

with open(args.output, 'w') as f:
    print(f"Preparing data for antibiotic: {antibiogram['antibiotic'].iloc[0]}", file=f)
    print("-" * 40, file=f)

# Remove antibiograms N phenotype
antibiogram = antibiogram[antibiogram['phenotype'] != 'N']
antibiogram = antibiogram[antibiogram['phenotype'] != 'I']
antibiogram = antibiogram[antibiogram['phenotype'] != 'NS']

antibiogram_unique = antibiogram.drop_duplicates(subset='id')

# Convert 'Element' to binary values
antibiotic_binary = antibiogram_unique['phenotype'].map({'S': 0, 'R': 1})
antibiotic_binary = antibiotic_binary.to_numpy().astype(np.int8)

filename = args.assemblies
# Open the file and read the lines into a list
files = os.listdir(filename)
assemblies = [os.path.basename(file).split('.fna')[0] for file in files]
# Find the indices of the elements in the sublist
indices = [assemblies.index(element) for element in antibiogram_unique['id']]

#New
# Resetting index for antibiogram_gentamycin_unique
antibiogram_unique.reset_index(drop=True, inplace=True)
antibiogram_unique['organism_name'] = antibiogram_unique['organism_name'].replace('E.coli and', 'Escherichia coli')
indices_antibiograms=antibiogram_unique.index.tolist()

#New
# Extract the genus (first word) from the organism_name column
antibiogram_unique['genus'] = antibiogram_unique['organism_name'].str.split().str[0]

# Find how many from each genus are there
genera_counts = antibiogram_unique['genus'].value_counts()

# Keep the indexes of each found genus
genera_indices = antibiogram_unique.groupby('genus').apply(lambda x: x.index.tolist())

# Report total number of rows
total_rows = antibiogram_unique.shape[0]


# Get the antibiotic from the first row
antibiotic_name = antibiogram_unique.loc[0, 'antibiotic']  # Assuming 'antibiotic' is the column name

# List of genera to include in the output
genera_list = ['Escherichia', 'Klebsiella', 'Acinetobacter', 'Staphylococcus', 'Pseudomonas', 'Enterobacter', 'Enterococcus']

# Example genera_counts data (replace with actual values from your dataset)
genera_counts = antibiogram_unique['genus'].value_counts()  # Replace df with your actual DataFrame

# Dictionary to store the results
summary_data = {
    'Total number': [total_rows],  # Total number of rows in your dataset
    'Antibiotic': [antibiotic_name]  # Add the antibiotic name as a row
}

# Populate the dictionary with genus counts, defaulting to 0 if the genus is not present
for genus in genera_list:
    summary_data[genus] = [genera_counts.get(genus, 0)]
    genus_res = antibiogram_unique.loc[antibiogram_unique["genus"]==genus,"phenotype"].value_counts()
    column_wanted=f"{genus}_NumberOfResistantStrains"
    summary_data[column_wanted] = [genus_res.get("R", 0)]
    
# Convert the dictionary to a DataFrame
summary_df = pd.DataFrame(summary_data)

summary_df.to_csv(f"./models/{antibiotic_name}_genera_summary.csv", index=False)



def load_data(kmer_size = 3):
    print(f"Loading data for kmer size {kmer_size}")

    # print to output file
    with open(args.output, 'a') as f:
        print(f"Loading data for kmer size {kmer_size}", file=f)
    
    # Make the final concatenated dataset 
    X_dna = np.load(f"/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/kmer_{kmer_size}_flatteneddna_.npy", allow_pickle=True)
    X_protein = np.load(f"/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/kmer_{kmer_size}_flattenedprotein_.npy", allow_pickle=True)
    X_rrna = np.load(f"/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/kmer_{kmer_size}_flattenedrrna_.npy", allow_pickle=True)

    X_rrna_subset = X_rrna[indices]
    X_dna_subset = X_dna[indices]
    X_protein_subset = X_protein[indices]
    
    del X_dna, X_protein, X_rrna

    X = np.concatenate((X_dna_subset, X_rrna_subset, X_protein_subset), axis=1)
    X = X.astype(np.int16)
    y = antibiotic_binary

     # Combine y and group labels for stratification
    stratification_labels = [f"{y_class}_{group}" for y_class, group in zip(y, antibiogram_unique['genus'])]
    
    # Perform the 80-20 split with both stratification and group constraints
    train_indices, test_indices = train_test_split(
        range(len(X)), 
        test_size=0.2, 
        random_state=42, 
        stratify=stratification_labels
    )
    
    # Use these indices to filter the original DataFrame
    #train_df = antibiogram_unique.iloc[train_indices]
    #test_df = antibiogram_unique.iloc[test_indices]
    
    # Optionally, you can retrieve specific columns if needed
    X_train, X_test = X[train_indices], X[test_indices]
    y_train, y_test = y[train_indices], y[test_indices]
    
    X_train = csr_matrix(X_train)
    X_test = csr_matrix(X_test)

    del X, y
    del X_dna_subset, X_protein_subset, X_rrna_subset

    return X_train, X_test, y_train, y_test, kmer_size, train_indices, test_indices


# -------------------------
# Hyperparameter Tuning
# -------------------------
param_grid = {
    'max_depth': [7, 9, 11, 13, 15, 17, 19, 21, 23],
    'num_boost_round': [300],
}

best_accuracy = 0
best_params = {}
best_kmer_size = 0

# Run XGBoost Classifier
for kmer in range(3, 6):
    # Load the data
    X_train, X_test, y_train, y_test, kmer_size, train_indices, test_indices = load_data(kmer)

    dtrain = xgb.DMatrix(X_train, label=y_train)
    dtest = xgb.DMatrix(X_test, label=y_test)

    del X_train
    del X_test

    for max_depth in track(param_grid['max_depth'], description='max_depth'):
        for num_boost_round in param_grid['num_boost_round']:
            params = {
                'eta': 0.1,
                'max_depth': max_depth,
                'objective': 'multi:softmax',
                'num_class': 2,
                'eval_metric': 'mlogloss'
            }
			# Assuming 'dtrain' is your full dataset in DMatrix format
            cv_results = xgb.cv(params=params,
                    dtrain=dtrain,
                    num_boost_round=num_boost_round,
                    nfold=5,  # 5-fold cross-validation
                    early_stopping_rounds=10,  # Early stopping based on CV results
                    verbose_eval=False)
            # Train the final model using all data with the optimal number of boosting rounds from CV
            best_num_boost_round = len(cv_results)
            model = xgb.train(params, dtrain, num_boost_round=best_num_boost_round, verbose_eval=False)
            y_pred = model.predict(dtest)
            
            # Calculate accuracy
            accuracy = accuracy_score(y_test, y_pred)

            # Print current parameters and accuracy
            print(f"max_depth: {max_depth}, actual_boost_rounds: {best_num_boost_round}, Accuracy: {accuracy}, Algorithm: xgboost")

            with open(args.output, 'a') as f:
                print(f"max_depth: {max_depth}, actual_boost_rounds: {best_num_boost_round}, Accuracy: {accuracy}, Algorithm: xgboost", file=f)
            
            # Check if the current model is the best so far
            if accuracy > best_accuracy:
                best_accuracy = accuracy
                best_params = {'algorithm': "xgboost", 'params': {'max_depth': max_depth, 'num_boost_round': best_num_boost_round}}
                best_kmer_size = kmer

    del dtrain, dtest
    # Print the best parameters and best accuracy
    print("-" * 40)
    print(f"Best parameters found for kmer {best_kmer_size}: {best_params}, Accuracy: {best_accuracy}\n")

    with open(args.output, 'a') as f:
        print("-" * 40, file=f)
        print(f"Best parameters found for kmer {best_kmer_size}: {best_params}, Accuracy: {best_accuracy}\n", file=f)


# Run Random Forest Classifier
param_grid = {
    'n_estimators': [100, 200, 400],
    'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
    'bootstrap': [True, False]
}


for kmer in range(3, 6):
    # Load the data for the current kmer size
    X_train, X_test, y_train, y_test, kmer_size, train_indices, test_indices = load_data(kmer)

    # Loop through the parameter grid
    for n_estimator in track(param_grid['n_estimators'], description=f"Progress for Kmer size {kmer}"):
        for max_depth in param_grid['max_depth']:
            for bootstrap in param_grid['bootstrap']:
                # Initialize the RandomForestClassifier with current parameters
                rf = RandomForestClassifier(n_estimators=n_estimator, max_depth=max_depth, bootstrap=bootstrap, n_jobs=njobs, verbose=0, random_state=random_state)
                
                # Perform cross-validation (5-fold) on the training data and compute average accuracy
                cv_scores = cross_val_score(rf, X_train, y_train, cv=5, scoring='accuracy')
                avg_cv_accuracy = np.mean(cv_scores)
                
                # Fit the model on the full training data (optional: only do this if you're fine-tuning on the full set after cross-validation)
                rf.fit(X_train, y_train)
                
                # Predict on the test data
                y_pred = rf.predict(X_test)
                
                # Calculate accuracy on the test data
                accuracy = accuracy_score(y_test, y_pred)
                
                # Print current parameters, CV accuracy, and test accuracy
                print(f"n_estimators: {n_estimator}, max_depth: {max_depth}, bootstrap: {bootstrap}, CV Accuracy: {avg_cv_accuracy}, Test Accuracy: {accuracy}, Algorithm: random_forest")
                
                with open(args.output, 'a') as f:
                    print(f"n_estimators: {n_estimator}, max_depth: {max_depth}, bootstrap: {bootstrap}, Accuracy: {accuracy}, Algorithm: random_forest", file=f)
                
                # Check if the current model is the best so far
                if accuracy > best_accuracy:
                    best_accuracy = accuracy
                    best_params = {'algorithm': "random_forest", 'params': {'n_estimators': n_estimator, 'max_depth': max_depth, 'bootstrap': bootstrap}}
                    best_kmer_size = kmer_size

    del X_train, X_test

    # Print the best parameters and best accuracy
    print(f"Best parameters found for kmer {best_kmer_size}: {best_params}, Accuracy: {best_accuracy}")

    with open(args.output, 'a') as f:
        print(f"Best parameters found for kmer {best_kmer_size}: {best_params}, Accuracy: {best_accuracy}", file=f)

print("-" * 40)
with open(args.output, 'a') as f:
    print("-" * 40, file=f)
    print(f"Best model : {best_params['algorithm']}\nBest parameters : {best_params}\nBest kmer: {best_kmer_size}", file=f)



# -------------------------
# Load the data with the best kmer size
# Fit the final model with the best parameters
# -------------------------
X_train, X_test, y_train, y_test, kmer_size, train_indices, test_indices = load_data(best_kmer_size)

if best_params['algorithm'] == "random_forest":
    # Fit the final model with the best parameters
    final_rf = RandomForestClassifier(**best_params['params'], n_jobs=njobs, verbose=0, random_state=random_state)
    final_rf.fit(X_train, y_train)
	
    # Predict on the test set with the final model
    y_pred = final_rf.predict(X_test)
    print("RANDOM FOREST")
    with open(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_RandomForest_kmer_{best_kmer_size}.pkl", 'wb') as f:
    	pickle.dump(final_rf, f)
    	
elif best_params['algorithm'] == "xgboost":
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dtest = xgb.DMatrix(X_test, label=y_test)

    del X_train
    # del X_test

    params = {
        'eta': 0.1,
        'max_depth': best_params['params']['max_depth'],
        'objective': 'multi:softmax',
        'num_class': 2,
        'eval_metric': 'mlogloss',
    }

    # fit XGBoost model with the best hyperparameters
    model = xgb.train(params, dtrain,
                        num_boost_round=best_params['params']['num_boost_round'],         # instead of best_params['params']['num_boost_round'],
                        verbose_eval=False)

    y_pred= model.predict(dtest)
    
    with open(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_xgboost_kmer_{best_kmer_size}.pkl", 'wb') as f:
    	pickle.dump(model, f)
    
# Print final performance metrics
print("\nFinal Model Performance:")
print("Algorithm: ", best_params['algorithm'])
print("Accuracy: ", accuracy_score(y_test, y_pred))
print("Precision: ", precision_score(y_test, y_pred))
print("Recall: ", recall_score(y_test, y_pred))
print("F1 Score: ", f1_score(y_test, y_pred))
print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
print("Classification Report:", classification_report(y_test, y_pred))

with open(args.output, 'a') as f:
    print("\nFinal Model Performance:", file=f)
    print("Algorithm: ", best_params['algorithm'], file=f)
    print("Accuracy: ", accuracy_score(y_test, y_pred), file=f)
    print("Precision: ", precision_score(y_test, y_pred), file=f)
    print("Recall: ", recall_score(y_test, y_pred), file=f)
    print("F1 Score: ", f1_score(y_test, y_pred), file=f)
    print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred), file=f)
    print("Classification Report:", classification_report(y_test, y_pred), file=f)




with open(args.output_final, 'a') as f:
    print("Best Algorithm: ", best_params['algorithm'], file=f)
    print("Best kmer size: ", best_kmer_size, file=f)
    print("Accuracy: ", accuracy_score(y_test, y_pred), file=f)
    print("Precision: ", precision_score(y_test, y_pred), file=f)
    print("Recall: ", recall_score(y_test, y_pred), file=f)
    print("F1 Score: ", f1_score(y_test, y_pred), file=f)
    print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred), file=f)
    print("Classification Report:", classification_report(y_test, y_pred), file=f)




plot_metrics(y_test, y_pred, filename=f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}.svg")


#NEW
antibiograms_test = antibiogram_unique.loc[test_indices]
antibiograms_test.reset_index(drop=True, inplace=True)

# Get the list of columns with at least one non-zero value
columns_to_keep = summary_df.columns[(summary_df != 0).any(axis=0)].tolist()

# Define columns to exclude
exclude_columns = ['Total number', 'Antibiotic']

# Get the columns with at least one non-zero value, excluding specified columns
columns_to_keep = summary_df.loc[:, (summary_df != 0).any(axis=0)].columns.difference(exclude_columns).tolist()

metrics_dict = {}
# Loop through each genus in the list
for genus in columns_to_keep:
    # Find the indices where the 'organism_name' column contains the genus
    indices_genus = antibiograms_test[antibiograms_test['genus'].str.contains(genus, case=False)].index.tolist()
    
    # If no indices found, skip this genus
    if not indices_genus:
        continue
    

    y_pred_genus = y_pred[indices_genus]
    y_test_genus = y_test[indices_genus]
    print(y_test_genus.shape)
    print(y_pred_genus.shape)
        # Extract the corresponding values from y_pred and y_test
    if len(y_pred_genus) == 0 :
        print(genus)
        continue 
    # Confusion matrix and classification report
    # Check if there is only one unique label in both y_test_genus and y_pred_genus
    if len(set(y_test_genus)) == 1 and len(set(y_pred_genus)) == 1:
        print("Skipping confusion matrix and classification report as only one label is present.")
    else:
    # Calculate confusion matrix and classification report if there are multiple labels
        conf_matrix = confusion_matrix(y_test_genus, y_pred_genus)
        print(conf_matrix)

        class_report = classification_report(y_test_genus, y_pred_genus, zero_division=1, output_dict=True)
        print(class_report)
    
    
        # Calculate the metrics for this genus
        accuracy = accuracy_score(y_test_genus, y_pred_genus)
        number = len(y_test_genus)
        # Store the metrics in the dictionary
        # Store the metrics in the dictionary, including confusion matrix and classification report
        metrics_dict[genus] = {
            'Number_Assemblies': number,
            'accuracy': accuracy,
            'confusion_matrix': conf_matrix.tolist(),  # Convert to list for easier display/storage
            'classification_report': class_report
        }



        plot_metrics(y_test_genus,y_pred_genus,filename=f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}_{genus}.svg")
        with open(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}_{genus}_output.txt", 'a') as f:
        	print("Best Algorithm: ", best_params['algorithm'], file=f)
        	print("Best kmer size: ", best_kmer_size, file=f)
        	print("Accuracy: ", accuracy_score(y_test_genus,y_pred_genus), file=f)
        	print("Precision: ", precision_score(y_test_genus,y_pred_genus), file=f)
        	print("Recall: ", recall_score(y_test_genus,y_pred_genus), file=f)
        	print("F1 Score: ", f1_score(y_test_genus,y_pred_genus), file=f)
        	print("Confusion Matrix:\n", confusion_matrix(y_test_genus,y_pred_genus), file=f)
        	print("Classification Report:", classification_report(y_test_genus,y_pred_genus), file=f)



    #New
# Function to flatten the classification report for better insertion into the DataFrame
def flatten_classification_report(report):
    flat_report = {}
    for key, value in report.items():
        if isinstance(value, dict):
            for metric, score in value.items():
                flat_report[f'{key}_{metric}'] = score
        else:
            flat_report[key] = value
    return flat_report

# Prepare data for the DataFrame
flattened_data = []

# Loop through the metrics_dict to prepare the data for DataFrame creation
for genus, metrics in metrics_dict.items():
    # Flatten the classification report so each metric becomes a column
    flattened_report = flatten_classification_report(metrics['classification_report'])
    
    # Prepare a dictionary for each genus containing all metrics
    genus_data = {
        'Number_Assemblies': metrics['Number_Assemblies'],  # Include the number of assemblies
        'accuracy': metrics['accuracy'],                   # Include the accuracy
        'confusion_matrix': metrics['confusion_matrix']     # Include the confusion matrix (optional, can remove if not needed)
    }
    
    # Merge the flattened report with the general metrics
    genus_data.update(flattened_report)
    
    # Append the genus-specific data to the flattened data list
    flattened_data.append(genus_data)
    
print("Now i will save the genera metrics")
print(metrics_dict.keys())
# Convert the list of dictionaries into a DataFrame, with genera as the index
metrics_df = pd.DataFrame(flattened_data, index=metrics_dict.keys())
print(metrics_df)

# Save
metrics_df.to_csv(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}_genus_metrics.csv")


# -------------------------
# Feature Importance function
# -------------------------
def feature_importance(model, X_test, best_params, best_kmer_size):
    # Get feature names
    features_dna=np.load(f'/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/kmer_{best_kmer_size}_flatteneddna_stats_.npy', allow_pickle=True)
    features_protein=np.load(f'/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/kmer_{best_kmer_size}_flattenedprotein_stats_.npy', allow_pickle=True)
    features_rrna=np.load(f'/mnt/raid1b/kdan_data/Paper/Machine_Learning/build_dataset/data/kmer_{best_kmer_size}_flattenedrrna_stats_.npy', allow_pickle=True)

    features_dna=pd.DataFrame(features_dna)
    features_dna['type']='promoter'
    features_rrna=pd.DataFrame(features_rrna)
    features_rrna['type']='rrna'

    features_rrna.iloc[:,2] = features_rrna.iloc[:,2].apply(lambda x: x[0])

    features_protein=pd.DataFrame(features_protein)
    features_protein['type']='protein'

    features = np.concatenate((features_dna, features_rrna, features_protein), axis=0)

    # Convert to DataFrame
    features_df = pd.DataFrame(features, columns=['Index', 'Sequence', 'Gene', 'Type'])
    features_df = features_df[['Sequence', 'Gene', 'Type']]

    # split the Gene column and keep the rest of the columns
    features_df['Gene'] = features_df['Gene'].str.split('_NODE', expand=True)[0].replace('>', '', regex=True)

    # concatenate Gene, Sequence and Type columns and make a list. If Gene is nan, replace it with 'nan'
    features_df['Gene_Sequence_Type'] = features_df['Gene'].fillna('nan') + '_' + features_df['Sequence'] + '_' + features_df['Type']

    print(f"NAs: {features_df['Gene_Sequence_Type'].isnull().sum()}")
    features_df['Index'] = range(len(features_df))


    features_list = list(features_df['Gene_Sequence_Type'])


    if best_params['algorithm'] == "random_forest":
        model=final_rf
        importances = final_rf.feature_importances_
        feature_importances = pd.DataFrame(importances, columns=['Score'])
        feature_importances = feature_importances.sort_values(by='Score', ascending=False)
        features_importance = pd.concat([feature_importances, features_df], axis=1)
        
        features_importance = features_importance[['Sequence', 'Gene', 'Type', 'Gene_Sequence_Type', 'Score']]

        # remove features with score less than 0.001 %
        features_importance[['Score_percentage']] = features_importance[['Score']] / features_importance[['Score']].sum() * 100
        #features_importance = features_importance[features_importance['Score_percentage'] > 0.001]
    
        features_importance['Index'] = features_importance.index
        data = features_importance[['Score_percentage']]
        data.index = features_importance.index
        
    elif best_params['algorithm'] == "xgboost":
        features_importance = model.get_score(importance_type='weight')
    
        keys = list(features_importance.keys())
        values = list(features_importance.values())

        data = pd.DataFrame(data=values, index=keys, columns=["Score"]).sort_values(by = "Score", ascending=True)

        data['Score_percentage'] = (data['Score'] / data['Score'].sum()) * 100

        # get 'data' row names, remove the 'f' prefix and convert to integer
        data.index = data.index.str[1:].astype(int)

#Check for dupliactes rows with same index 
        duplicates = features_df['Index'].duplicated()

        duplicate_rows = features_df[duplicates]
        duplicate_rows
        
        
        features_importance = features_df[features_df['Index'].isin(data.index)]
        features_importance['Score_percentage'] = features_importance.index.map(lambda idx: data['Score_percentage'].get(idx, 0))
        
    
    # SHAP (SHapley Additive exPlanations) values analysis
    explainer = shap.TreeExplainer(model, feature_names=features_list)

    # convert X_test to a pandas DataFrame
    X_test_df = pd.DataFrame(X_test.toarray())
    #X_test_df = pd.DataFrame(X_test)
    shap_values = explainer.shap_values(X_test_df)

		# Assuming shap_values_rf_class_1 is already computed and features_df is your feature list
    shap_values_rf_class_1 = shap_values[:, :, 1]  # Extract the second value (index 1) for each feature

# Convert the SHAP values for class 1 to a DataFrame
    shap_values_rf_df = pd.DataFrame(shap_values_rf_class_1, columns=features_list)

# Calculate the mean absolute SHAP value for each feature (for importance)
    mean_abs_shap_values_rf = shap_values_rf_df.abs().mean()

# Calculate the mean SHAP value for each feature (for direction)
    mean_shap_direction_rf = shap_values_rf_df.mean()

# Determine if the majority of SHAP values for each feature are positive or negative
# If mean direction is positive, the feature contributes positively on average
# If mean direction is negative, the feature contributes negatively on average
    shap_direction_label = mean_shap_direction_rf.apply(lambda x: 'Positive' if x > 0 else 'Negative')

# Create a DataFrame to store feature names, their importance, and direction
    features_importance_shap = pd.DataFrame({
    	'Feature': features_list,
    	'Mean_Abs_SHAP_Value': mean_abs_shap_values_rf,
    	'Mean_SHAP_Direction': mean_shap_direction_rf,  # Mean SHAP direction (numeric)
    	'Direction_Label': shap_direction_label          # Label: Positive or Negative
    })

# Sort the DataFrame by SHAP value importance (descending order)
    features_importance_shap = 	features_importance_shap.sort_values(by='Mean_Abs_SHAP_Value', ascending=False)

# Filter the DataFrame to keep only features with SHAP importance greater than or equal to 0.00001
    features_importance_shap_filtered = features_importance_shap[features_importance_shap['Mean_Abs_SHAP_Value'] >= 0.00001]

# Reset the index to make the index a column
    features_importance_shap_filtered = features_importance_shap_filtered.reset_index(drop=True)

	# Display the filtered DataFrame with the SHAP value direction
    #print(features_importance_shap_filtered)



# Step 2: Merge the two DataFrames based on the index (from feature_importance_filtered) and 'gene_sequence_type' (from features_df_filtered)
# Assuming that the reset index from feature_importance_filtered corresponds to 'gene_sequence_type' in features_df_filtered
    merged_df_rf = pd.merge(features_importance_shap_filtered, 
                     features_importance, 
                     left_on='Feature',  # 'index' from feature_importance_filtered (after resetting)
                     right_on='Gene_Sequence_Type',
                    how='left')
# 'gene_sequence_type' from features_df_filtered
	# Step 2: Find features in features_importance_shap_filtered that do not match with features_importance
    missing_features = merged_df_rf[merged_df_rf['Index'].isna()]['Feature']

# Step 3: For missing features, complete the columns from the dataframe features_df
# Filter the features_df to find the corresponding missing features
    additional_rows = features_df[features_df['Gene_Sequence_Type'].isin(missing_features)].copy()

# Add a new column 'Score_percentage' with value 0 for these additional rows
    additional_rows['Score_percentage'] = 0

# Step 3: Use loc to replace rows in merged_df_rf where 'Index' is NaN
# First, match the index of the rows to replace in merged_df_rf
    na_indices = merged_df_rf[merged_df_rf['Index'].isna()].index

# Then, replace these rows with the corresponding rows from additional_rows
    merged_df_rf.loc[na_indices, additional_rows.columns] = additional_rows.values

# Display the merged DataFrame
    #print(merged_df_rf)


    merged_df_rf.to_csv(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}_feature_importance_filtered.csv")
    merged_df_rf_names = merged_df_rf['Gene_Sequence_Type'].tolist()

# Open a file in write mode ('w'), or create it if it doesn't exist
    with open(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_feature_names_filtered_{best_params['algorithm']}_kmer_{best_kmer_size}.txt", 'w') as f:
# Iterate through the list and write each item to the file
    	for item in merged_df_rf_names:
    		f.write(f"{item}\n")  # Add a newline character after each item

#Me vasi to filtro poyu evala parapano (Thelei alagi auto)
    shap_values_filtered = shap_values[:, merged_df_rf['Index'].tolist(), :]
    np.save(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}_shap_values.npy", shap_values_filtered)

# 3. Keep only the X_test values with the highest importance
    X_test_df_filtered = X_test_df.iloc[:, merged_df_rf['Index'].tolist()]
    X_test_df_filtered.to_csv(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}_X_test_filtered.csv")

# Get the list of columns with at least one non-zero value
    columns_to_keep = summary_df.columns[(summary_df != 0).any(axis=0)].tolist()

# Define columns to exclude
    exclude_columns = ['Total number', 'Antibiotic']

# Get the columns with at least one non-zero value, excluding specified columns
    columns_to_keep = summary_df.loc[:, (summary_df != 0).any(axis=0)].columns.difference(exclude_columns).tolist()


#Create for each genus the shap values 

    for genus in columns_to_keep:
    # Find the indices where the 'genus' column contains the genus (case-insensitive)
        indices_genus = antibiograms_test[antibiograms_test['genus'].str.contains(genus, case=False)].index.tolist()
    
    # If no indices found, skip this genus
        if not indices_genus:
            continue
    
    # Extract the corresponding values from y_pred, y_test, and X_test for the genus
        y_pred_genus = y_pred[indices_genus]
        y_test_genus = y_test[indices_genus]
        X_test_genus_df = X_test_df.iloc[indices_genus, :]
    
        # Compute the SHAP values for this genus subset
        shap_values_genus = explainer.shap_values(X_test_genus_df)


		# Assuming shap_values_rf_class_1 is already computed and feautures_df_list is your feature list
        shap_values_rf_class_1 = shap_values_genus[:, :, 1]  # Extract the second value (index 1) for each feature

	# Convert the SHAP values for class 1 to a DataFrame
        shap_values_rf_df = pd.DataFrame(shap_values_rf_class_1, columns=features_list)

	# Calculate the mean absolute SHAP value for each feature (for importance)
        mean_abs_shap_values_rf = shap_values_rf_df.abs().mean()

	# Calculate the mean SHAP value for each feature (for direction)
        mean_shap_direction_rf = shap_values_rf_df.mean()

	# Determine if the majority of SHAP values for each feature are positive or negative
	# If mean direction is positive, the feature contributes positively on average
	# If mean direction is negative, the feature contributes negatively on average
        shap_direction_label = mean_shap_direction_rf.apply(lambda x: 'Positive' if x > 0 else 'Negative')

	# Create a DataFrame to store feature names, their importance, and direction
        features_importance_genus_shap = pd.DataFrame({
        	'Feature': features_list,
        	'Mean_Abs_SHAP_Value': mean_abs_shap_values_rf,
        	'Mean_SHAP_Direction': mean_shap_direction_rf,  # Mean SHAP direction (numeric)
        	'Direction_Label': shap_direction_label          # Label: Positive or Negative
        })

	# Sort the DataFrame by SHAP value importance (descending order)
        features_importance_genus_shap = 	features_importance_genus_shap.sort_values(by='Mean_Abs_SHAP_Value', ascending=False)

	# Filter the DataFrame to keep only features with SHAP importance greater than or equal to 0.00001
        features_importance_genus_shap_filtered = features_importance_genus_shap[features_importance_genus_shap['Mean_Abs_SHAP_Value'] >= 0.00001]

	# Reset the index to make the index a column
        features_importance_genus_shap_filtered = features_importance_genus_shap_filtered.reset_index(drop=True)

		# Display the filtered DataFrame with the SHAP value direction
		#print(features_importance_genus_shap_filtered)

	
	# Step 2: Merge the two DataFrames based on the index (from feature_importance_filtered) and 'gene_sequence_type' (from features_df_filtered)
	# Assuming that the reset index from feature_importance_filtered corresponds to 'gene_sequence_type' in features_df_filtered
        merged_df_rf = pd.merge(features_importance_genus_shap_filtered, 
                         features_importance, 
                         left_on='Feature',  # 'index' from feature_importance_filtered (after resetting)
                         right_on='Gene_Sequence_Type',
                        how='left')
	# 'gene_sequence_type' from features_df_filtered

	# Display the merged DataFrame
		#print(merged_df_rf)
# 'gene_sequence_type' from features_df_filtered
	# Step 2: Find features in features_importance_shap_filtered that do not match with features_importance
        missing_features = merged_df_rf[merged_df_rf['Index'].isna()]['Feature']

# Step 3: For missing features, complete the columns from the dataframe features_df
# Filter the features_df to find the corresponding missing features
        additional_rows =     features_df[features_df['Gene_Sequence_Type'].isin(missing_features)].copy()

# Add a new column 'Score_percentage' with value 0 for these additional rows
        additional_rows['Score_percentage'] = 0

# Step 3: Use loc to replace rows in merged_df_rf where 'Index' is NaN
# First, match the index of the rows to replace in merged_df_rf
        na_indices = merged_df_rf[merged_df_rf['Index'].isna()].index

# Then, replace these rows with the corresponding rows from additional_rows
        merged_df_rf.loc[na_indices, additional_rows.columns] = additional_rows.values

	
        merged_df_rf.to_csv(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}_feature_importance_filtered_{genus}.csv")
        merged_df_rf_names = merged_df_rf['Gene_Sequence_Type'].tolist()

	# Open a file in write mode ('w'), or create it if it doesn't exist
        with open(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_feature_names_filtered_{best_params['algorithm']}_kmer_{best_kmer_size}_{genus}.txt", 'w') as f:
    # Iterate through the list and write each item to the file
        	for item in merged_df_rf_names:
        		f.write(f"{item}\n")  # Add a newline character after each item

	#Me vasi to filtro poyu evala parapano (Thelei alagi auto)
        shap_values_genus_filtered = shap_values_genus[:, merged_df_rf['Index'].tolist(), :]
        np.save(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}_shap_values_{genus}.npy", shap_values_genus_filtered)
	
	# 3. Keep only the X_test values with the highest importance
        X_test_genus_df_filtered = X_test_genus_df.iloc[:, merged_df_rf['Index'].tolist()]
        X_test_genus_df_filtered.to_csv(f"./models/{antibiogram['antibiotic'].iloc[0].split(' ')[0]}_{best_params['algorithm']}_kmer_{best_kmer_size}_X_test_filtered_{genus}.csv")

    print("Feature importance analysis completed")
    with open(args.output, 'a') as f:
        print("Feature importance analysis completed", file=f)

	

if best_params['algorithm'] == "random_forest":
    feature_importance(final_rf, X_test, best_params, best_kmer_size)
elif best_params['algorithm'] == "xgboost":
    feature_importance(model, X_test, best_params, best_kmer_size)
