#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)  # ignore all future warnings
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, accuracy_score, precision_score, recall_score
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
import pickle

""" Train (fit) each model on the training set using the best
    hyperparameters found by hypertuning_sno_pseudo_simple_models. Pickle
    these fitted models (into .sav files) so that they can be reused after without the
    need to retrain them all over again."""

fold_num = snakemake.wildcards.fold_num

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


# Get training set and filter to keep only expressed C/D and C/D pseudogenes
X_train = pd.read_csv(snakemake.input.X_train_scaled, sep='\t').drop(columns=['gene_id']).reset_index(drop=True)  # already filtered during training
y_train = pd.read_csv(snakemake.input.y_train, sep='\t').drop(columns=['gene_id'])
y_train = y_train[y_train['target'] != 0].reset_index(drop=True)


# Train over given fold (fold_num) in stratified 10-fold CV
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=rs)
fold_dict = {str(fold_index+1): [train_index, test_index]
            for fold_index, (train_index, test_index) in
            enumerate(skf.split(X_train, y_train))}

# Get train and test (CV) sets for given fold
print(X_train)
print(y_train)
train_index = fold_dict[fold_num][0]
test_index = fold_dict[fold_num][1]
X_fold_train = X_train.loc[train_index]
y_fold_train = y_train.loc[train_index]
X_fold_test = X_train.loc[test_index]
y_fold_test = y_train.loc[test_index]


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
model.fit(X_fold_train, y_fold_train.values.ravel())

# Test model on test set (heldout tenth of all training example)
y_pred = model.predict(X_fold_test)
print(y_pred)

acc = {}
acc['model'] = snakemake.wildcards.simple_models
acc['training_accuracy'] = accuracy_score(y_fold_test.target, y_pred)
acc['precision'] = precision_score(y_fold_test.target, y_pred, average='macro')
acc['recall'] = recall_score(y_fold_test.target, y_pred, average='macro')
acc['training_f1_score'] = f1_score(y_fold_test.target, y_pred, average='macro')
acc_df = pd.DataFrame(acc, index=[0])
print(acc_df)
acc_df.to_csv(snakemake.output.train_metrics, sep='\t', index=False)

# Pickle the model as a .sav file ('wb' for write in binary)
pickle.dump(model, open(snakemake.output.model, 'wb'))
