This step uses pre-trained machine learning models to predict antibiotic resistance profiles for newly annotated assemblies.
📥 Step 1: Download Models

Download models_and_shap_values.zip from Zenodo and extract it inside the predict_newAssemblies/ folder:

predict_newAssemblies/
├── models_and_shap_values/
│   ├── trained_models/
│   ├── shap_values/
│   └── ...

🚀 Step 2: Launch Jupyter Notebook

Start Jupyter Notebook:

jupyter notebook

Then open the notebook located in predict_newAssemblies/ and run all cells to:

    Load the trained models

    Process the k-mer encoded assemblies

    Predict the antibiotic resistance profile

    Output the in silico antibiogram for each sample
