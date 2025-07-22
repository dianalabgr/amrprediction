mkdir piperacillin_tazobactam 
cd piperacillin_tazobactam
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/piperacillin_tazobactam/piperacillin_tazobactam_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir sulfisoxazole 
cd sulfisoxazole
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/sulfisoxazole/sulfisoxazole_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir tetracycline 
cd tetracycline
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/tetracycline/tetracycline_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir tigecycline
cd tigecycline
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/tigecycline/tigecycline_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir tobramycin 
cd tobramycin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/tobramycin/tobramycin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir trimethoprim 
cd trimethoprim
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/trimethoprim/trimethoprim_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir trimethoprim_sulfamethoxazole 
cd trimethoprim_sulfamethoxazole
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/trimethoprim_sulfamethoxazole/trimethoprim_sulfamethoxazole_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir vancomycin
cd vancomycin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/vancomycin/vancomycin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

