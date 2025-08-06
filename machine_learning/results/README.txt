Run from conda python3.12.5
cd folder_antibiotic
conda activate py312
mkdir models
python ../kmer_rf_xgboost_final3.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/amoxicillin-clavulanicAcid/amoxicillin_clavulanic_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
python ../kmer_rf_xgboost_final3.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/gentamicin/gentamicin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
(START RUNNING 11/11/2024 12:28pm)

#python kmer_rf_xgboost_final.py --antibiogram ../amoxicillin_clavulanic_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/assemblies.txt --output ./kmer_rf_xgboost_output.txt

#python kmer_rf_xgboost_final3.py --antibiogram ../amoxicillin_clavulanic_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt

#python ../kmer_rf_xgboost_final3.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/amoxicillin-clavulanicAcid/amoxicillin_clavulanic_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
