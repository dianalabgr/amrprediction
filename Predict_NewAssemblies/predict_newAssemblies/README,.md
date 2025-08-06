To generate in silico antibiograms for the new assemblies using pre-trained models:
1. Download Models

Download models_and_shap_values.zip from Zenodo and unzip it into the predict_newAssemblies/ directory.

predict_newAssemblies/
├── models_and_shap_values/
│   ├── trained_models/
│   ├── shap_values/
│   └── ...

2. Launch Jupyter Notebook

From your terminal, start Jupyter:

jupyter notebook

Open the notebook located in the predict_newAssemblies/ folder and run all cells to:

    Load the trained machine learning models

    Generate predictions for the assemblies

    Output the predicted antibiograms
