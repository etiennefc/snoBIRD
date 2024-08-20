#!/usr/bin/python3
import pandas as pd

# snoRNA pseudogenes are not considered here 
# for the sake of the confusion matrix

snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t')
snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t')
infernal_rfam = pd.read_csv(snakemake.input.infernal_rfam, sep='\t')

output_snoreport = snakemake.output.matrix_snoreport
output_snoscan = snakemake.output.matrix_snoscan
output_infernal_rfam = snakemake.output.matrix_infernal_rfam

def confusion(df, y_true, y_pred, true_label, false_label, output):
    # Return the confusion matrix of a binary 
    # classification results in df in columns 
    # y_true and y_pred
    tp = len(df[(df[y_true] == true_label) & (df[y_pred] == true_label)])
    tn = len(df[(df[y_true] == false_label) & (df[y_pred] == false_label)])
    fp = len(df[(df[y_true] == false_label) & (df[y_pred] == true_label)])
    fn = len(df[(df[y_true] == true_label) & (df[y_pred] == false_label)])
    df = pd.DataFrame([[tp, fp], [fn, tn]], columns=['positives', 'negatives'], index=['positives', 'negatives'])
    df.columns = pd.MultiIndex.from_product([["Actual values"], df.columns])
    df.index = pd.MultiIndex.from_product([["Predicted values"], df.index])
    df.to_csv(output, sep='\t')

confusion(snoreport, 'target', 'snoreport2_prediction', 'expressed_CD_snoRNA', 'other', output_snoreport)
confusion(snoscan, 'target', 'snoscan_prediction', 'expressed_CD_snoRNA', 'other', output_snoscan)
confusion(infernal_rfam, 'target', 'infernal_rfam_prediction', 'expressed_CD_snoRNA', 'other', output_infernal_rfam)
