spg_fresh_blob_202104
==============================

Lagrangian particle tracking of Subpolar Gyre freshening using Viking20x and Ocean Parcels

Code for the particle tracking and analysis for the Fox et al. 2022 paper "Exceptional freshening and cooling in the eastern subpolar North Atlantic caused by reduced Labrador Sea surface heat loss".

Particle tracking uses OceanParcels. All analysis done in jupyter notebooks.

### Code to run the particle tracking

    notebooks/037_afox_RunParcels*
    The large trajectory output files are stored at GEOMAR: https://hdl.handle.net/20.500.12085/830c72af-b5ca-44ac-8357-3173392f402b

### Extract pathways, times and distances

    notebooks/exploratory/115_afox*.ipynb
    
Extracts the pathways, transit times and distances from the full track dataset. Stores in a reduced form with start and end points

### Labrador Sea water mass analysis

    notebooks/exploratory/2*.ipynb
    
Calculates the water mass budget in density space from the VIKING20X data: inflow, outflow, volume change, transformations due to surface fluxes of heat and freshwater.

### Sums and means of track metrics

    notebooks/exploratory/3*.ipynb
    
Calculates sums and means from the outputs of notebooks/exploratory/115_afox*.ipynb by source and pathway.

### Plots for the paper

    notebooks/exploratory/4*.ipynb




Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
