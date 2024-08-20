#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)  # ignore all future warnings
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
import pickle

""" Train (fit) each model on the training set using the best
    hyperparameters found by hyperparameter_tuning_cv_simple_models. Pickle
    these fitted models (into .sav files) so that they can be reused after without the
    need to retrain them all over again."""

# Get best hyperparameters per model
hyperparams_df = pd.read_csv(snakemake.input.best_hyperparameters, sep='\t')
hyperparams_df = hyperparams_df.drop('accuracy_cv', axis=1)
rs = snakemake.params.random_state

def df_to_params(df):
    """ Convert a one-line dataframe into a dict of params and their value. The
        column name corresponds to the key and the value corresponds to
        the value of that param (ex: 'max_depth': 2 where max_depth was the column
        name and 2 was the value in the df)."""
    cols = list(df.columns)
    params = {}
    for col in cols:
        value = df.loc[0, col]
        params[col] = value
    return params

hyperparams = df_to_params(hyperparams_df)


# Get training set
X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id')
X_train = X_train[['box_score_norm', 'structure_mfe_norm', 'terminal_stem_mfe_norm', 'length_norm']]
y_train = pd.read_csv(snakemake.input.y_train, sep='\t').drop(columns=['gene_id'])

# Instantiate the model defined by the 'simple_models' wildcard using the best hyperparameters
# specific to each model (logreg, svc, rf, gbm, knn)
if snakemake.wildcards.simple_models == "logreg":
    model = LogisticRegression(C=hyperparams['C'], solver=hyperparams['solver'],
                                random_state=rs, max_iter=500)
elif snakemake.wildcards.simple_models == "svc":
    model = svm.SVC(C=hyperparams['C'], degree=hyperparams['degree'],
                    gamma=hyperparams['gamma'], kernel=hyperparams['kernel'],
                    random_state=rs)
elif snakemake.wildcards.simple_models == "rf":
    model = RandomForestClassifier(max_depth=hyperparams['max_depth'],
                min_samples_leaf=hyperparams['min_samples_leaf'],
                min_samples_split=hyperparams['min_samples_split'],
                n_estimators=hyperparams['n_estimators'], random_state=rs)
elif snakemake.wildcards.simple_models == "knn":
    model = KNeighborsClassifier(n_neighbors=hyperparams['n_neighbors'],
                weights=hyperparams['weights'],
                leaf_size=hyperparams['leaf_size'], p=hyperparams['p'])
else:
    model = GradientBoostingClassifier(loss=hyperparams['loss'],
                max_depth=hyperparams['max_depth'],
                min_samples_leaf=hyperparams['min_samples_leaf'],
                min_samples_split=hyperparams['min_samples_split'],
                n_estimators=hyperparams['n_estimators'], random_state=rs)

# Train model and save training accuracy to df
model.fit(X_train, y_train.values.ravel())
print(snakemake.wildcards.simple_models)
print(model.score(X_train, y_train.values.ravel()))

acc = {}
acc[snakemake.wildcards.simple_models+'_training_accuracy'] = model.score(X_train, y_train.values.ravel())
acc_df = pd.DataFrame(acc, index=[0])
acc_df.to_csv(snakemake.output.training_accuracy, sep='\t', index=False)

# Pickle the model as a .sav file ('wb' for write in binary)
pickle.dump(model, open(snakemake.output.pickled_trained_model, 'wb'))