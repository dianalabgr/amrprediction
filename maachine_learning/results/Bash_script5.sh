mkdir amikacin 
cd amikacin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/amikacin/amikacin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir amoxicillin-clavulanicAcid 
cd amoxicillin-clavulanicAcid
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/amoxicillin-clavulanicAcid/amoxicillin_clavulanic_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir ampicillin 
cd ampicillin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/Ampicillin/ampicillin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir ampicillin-sulbactam 
cd ampicillin-sulbactam
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/ampicillin-sulbactam/ampicillin_sulbactam_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir azithromycin 
cd azithromycin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/azithromycin/azithromycin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir aztreonam 
cd aztreonam
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/aztreonam/aztreonam_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir cefazolin 
cd cefazolin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/cefazolin/cefazolin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir cefepime 
cd cefepime
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/cefepime/cefepime_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..



