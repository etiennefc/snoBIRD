#!/usr/bin/python3
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV, RandomizedSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

""" Tune the hyperparameters of each models (Logistic regression, Support
    vector classifier, Random Forest and Gradient boosting classifier)
    before even training them, using GridSearchCV with stratified k-fold on the
    cross-validation (X_cv) set."""

X_cv = pd.read_csv(snakemake.input.X_cv, sep='\t', index_col='gene_id')
y_cv = pd.read_csv(snakemake.input.y_cv, sep='\t').drop(columns=['gene_id'])
rs = snakemake.params.random_state
output = snakemake.output.best_hyperparameters


# Instantiate the model defined by the 'models' wildcard
if snakemake.wildcards.simple_models == "logreg":
    model = LogisticRegression(random_state=rs, max_iter=500)
elif snakemake.wildcards.simple_models == "svc":
    model = svm.SVC(random_state=rs)
elif snakemake.wildcards.simple_models == "rf":
    model = RandomForestClassifier(random_state=rs)
elif snakemake.wildcards.simple_models == "knn":
    model = KNeighborsClassifier()
else:
    model = GradientBoostingClassifier(random_state=rs)

# Configure the cross-validation strategy (StratifiedKFold where k=3)
cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=rs)

# Define search space (i.e. in which range of params the GridSearch happens) per model
space = snakemake.params.hyperparameter_space

# Execute the gridsearch per model on the CV set
search = GridSearchCV(estimator=model, param_grid=space,
                        cv=cv, scoring="f1_macro")
search.fit(X_cv, y_cv.values.ravel())
print(snakemake.wildcards.simple_models)
print(search.best_score_)
print(search.best_params_)

# Return the best hyperparameters found by GridSearchCV, and the accuracy of each model
# fitted on the CV set with these hyperparameters into a dataframe
params_dict = search.best_params_
params_dict['accuracy_cv'] = search.best_score_
params_df = pd.DataFrame(params_dict, index=[0])
params_df.to_csv(output, sep='\t', index=False)