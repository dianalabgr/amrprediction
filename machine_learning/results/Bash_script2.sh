mkdir chloramphenicol 
cd chloramphenicol
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/chloramphenicol/chloramphenicol_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir ciprofloxacin 
cd ciprofloxacin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/ciprofloxacin/ciprofloxacin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir clindamycin 
cd clindamycin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/clindamycin/clindamycin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir colistin 
cd colistin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/colistin/colistin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir doripenem 
cd doripenem
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/doripenem/doripenem_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir ertapenem 
cd ertapenem
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/ertapenem/ertapenem_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir erythromycin 
cd erythromycin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/erythromycin/erythromycin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir gentamicin 
cd gentamicin
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/gentamicin/gentamicin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir imipenem 
cd imipenem
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/imipenem/imipenem_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..
