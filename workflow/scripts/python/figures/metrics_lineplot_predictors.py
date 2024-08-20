#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, recall_score, matthews_corrcoef

#snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t')
#snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t')
#infernal_rfam = pd.read_csv(snakemake.input.infernal_rfam, sep='\t')
snoreport = pd.read_csv('~/Desktop/Etienne/cd_predictor/workflow/results/predictions/snoreport2/fixed_length_194nt/test_predictions.tsv', sep='\t')
snoscan = pd.read_csv('~/Desktop/Etienne/cd_predictor/workflow/results/predictions/snoscan/fixed_length_194nt/test_predictions.tsv', sep='\t')
infernal_rfam = pd.read_csv('~/Desktop/Etienne/cd_predictor/workflow/results/predictions/infernal_rfam/fixed_length_194nt/test_predictions.tsv', sep='\t')
snoBIRD = pd.read_csv('results/predictions/snoBIRD/transformer/194/3e-5_3e-6_32_4_data_aug_1_ratio/transformer_2_classes_LR_schedule_test_predictions_194nt_fold_8.tsv', sep='\t')
snoBIRD['snoBIRD_prediction'] = snoBIRD['y_pred'].replace({0: 'other', 1: 'expressed_CD_snoRNA'})

# with snoRNA_pseudogene included
#metrics_transformer = [['accuracy', 0.9661, 'snoBIRD'], 
#                ['precision', 0.9728, 'snoBIRD'], 
#                ['recall', 0.8165, 'snoBIRD'],
#                ['f1_score', 0.8878, 'snoBIRD']] 
# without snoRNA_pseudogene_included
#metrics_transformer = [['accuracy', 0.9769, 'snoBIRD'], 
#                ['precision', 0.9697, 'snoBIRD'], 
#                ['recall', 0.8653, 'snoBIRD'],
#                ['f1_score', 0.9145, 'snoBIRD']] 
#transformer = pd.DataFrame(metrics_transformer, columns=['score_name', 'score_value', 'predictor'])


output_all = snakemake.output.dotplot_all
output_cd = snakemake.output.dotplot_expressed_cd
output_pseudo = snakemake.output.dotplot_pseudo
color_dict = snakemake.params.predictors_colors

# Merge dfs
df = snoreport.merge(snoscan[['gene_id', 'snoscan_prediction']], 
                        how='left', on='gene_id')
df = df.merge(infernal_rfam[['gene_id', 'infernal_rfam_prediction']], 
                        how='left', on='gene_id')
df = df.merge(snoBIRD[['gene_id', 'snoBIRD_prediction']], 
                        how='left', on='gene_id')

# Predict on 1 window per sno (no data aug, which is less representative for snoBIRD 
# because it uses consecutive windows to predict a C/D)
#df = df[~df['gene_id'].str.contains(r'_[AB]\d$')]

# Compute the different metrics for each existing CD predictor
titles = ['Test metrics on all C/D', 'Test metrics on expressed C/D only', 'Test metrics on C/D pseudogenes only']
for i, subset in enumerate([output_all, output_cd, output_pseudo]):
    dfs = []
    tempo_df = df.copy()
    if i == 0:  # consider all sno as "expressed"
        tempo_df['gene_biotype'] = tempo_df['gene_biotype'].replace('snoRNA_pseudogene', 'expressed_CD_snoRNA')
    elif i == 1:
        tempo_df = tempo_df[tempo_df['gene_biotype'] != 'snoRNA_pseudogene']
    elif i == 2:
        tempo_df = tempo_df[tempo_df['gene_biotype'] != 'expressed_CD_snoRNA']
        # Convert pseudo to "expressed for the sake of metrics computation"
        tempo_df['gene_biotype'] = tempo_df['gene_biotype'].replace('snoRNA_pseudogene', 'expressed_CD_snoRNA')
    tempo_df['target'] = tempo_df['target'].replace('snoRNA_pseudogene', 'expressed_CD_snoRNA')
    for model in color_dict.keys():
        tempo_df[f'{model}_prediction'] = tempo_df[f'{model}_prediction'].replace('snoRNA_pseudogene', 'expressed_CD_snoRNA')
        accuracy = accuracy_score(tempo_df['target'], tempo_df[f'{model}_prediction'])
        precision, recall, f1_score, _ = precision_recall_fscore_support(
                                        tempo_df['target'], tempo_df[f'{model}_prediction'], 
                                        pos_label='expressed_CD_snoRNA', average='binary')
        #mcc = matthews_corrcoef(tempo_df['target'], tempo_df[f'{model}_prediction'])
        temp_df = pd.DataFrame([['accuracy', accuracy, model], 
                                ['precision', precision, model], 
                                ['recall', recall, model], 
                                ['f1_score', f1_score, model]],
                                columns=['score_name', 'score_value', 'predictor'])
        if model == 'infernal_rfam':
            print(tempo_df[(tempo_df['target'] != 'expressed_CD_snoRNA') & (tempo_df[f'{model}_prediction'] == 'expressed_CD_snoRNA')][['target', f'{model}_prediction', 'gene_id']])
        print(model, titles[i])
        print('accuracy: ', accuracy)
        print('precision: ', precision)
        print('recall: ', recall)
        print('f1_score: ', f1_score)
        dfs.append(temp_df)

    final_df = pd.concat(dfs).reset_index(drop=True)
    #print(final_df)

    # Create graph
    ft.lineplot(final_df, 'score_name', 'score_value', 'predictor', 
                'Metrics', 'Metrics value', titles[i], color_dict, subset, markersize=30, linewidth=16)

