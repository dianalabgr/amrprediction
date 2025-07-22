#mkdir cefotaxime 
#cd cefotaxime 
#python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/cefotaxime/cefotaxime_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt


#cd .. 
mkdir cefotetan 
cd cefotetan 
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/cefotetan/cefotetan_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

mkdir cefoxitin 
cd cefoxitin 
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/cefoxitin/cefoxitin_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd .. 

mkdir ceftazidime 
cd ceftazidime 
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/ceftazidime/ceftazidime_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd .. 

mkdir cefpodoxime
cd cefpodoxime
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/cefpodoxime/cefpodoxime_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd .. 

mkdir ceftriaxone 
cd ceftriaxone 
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/ceftriaxone/ceftriaxone_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd .. 

mkdir cefuroxime 
cd cefuroxime
python ../kmer_rf_xgboost_final3_correct.py --antibiogram /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/cefuroxime/cefuroxime_QCed_metadata.csv --cores 30 --assemblies /mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/filtered_assemblies/ --output ./kmer_rf_xgboost_output.txt --output_final ./models/kmer_rf_xgboost_output_final.txt
cd ..

