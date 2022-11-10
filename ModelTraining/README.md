# Model Training

- The model was trained on the mega dataset created from the 4 datasets.
- But first, I removed the outliers using IQR method based on the molecular weight column. 
- This removed 494 compounds
- I tried diffrent models, but the model that showed the best performance was a stack model that conists of:
        1. DecisionTreeRegressor(max_depth=9)
        2. HistGradientBoostingRegressor(max_depth=21)
        3. RandomForestRegressor(n_estimators=29)

The model evalauation is:
        1. R^2 equals 0.775817
        2. MAE equals 0.731621
        3. RMSD equals 1.094537