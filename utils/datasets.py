import pandas as pd
import numpy as np

def load_all_data():
    fhr = pd.read_csv('./annotations/fhr-annotations-raw.2Mar25.tsv', sep='\t', index_col=0)
    fhr = fhr[fhr['risk'] != -1]
    all_label = fhr['risk'].map({0:'SR', 1:'GHR', 2:'FHR'})

    all_cn_arm = pd.read_csv('./matrices/broad_cn_matrix_integer.tsv', sep='\t', index_col=0)
    all_cn_seg = pd.read_csv('./matrices/segment_cn_matrix_uncorrelated.tsv', sep='\t', index_col=0)
    all_rna = pd.read_csv('./matrices/gene_exp_logcpm_proteincoding.tsv', sep='\t', index_col=0)

    annot = pd.read_csv('./annotations/ensgid-autosomal-proteincoding.txt', sep='\t', index_col=0)
    gene_name_rna = annot.loc[all_rna.columns, 'Gene name']
    duplicated = (gene_name_rna.duplicated(keep='first')) & (~gene_name_rna.isna())
    missing = gene_name_rna.isna()
    gene_name_id_rna = np.where(missing | duplicated, all_rna.columns, gene_name_rna)
    all_rna.columns = gene_name_id_rna

    all_mut = pd.read_csv('./matrices/gene_mut_matrix.tsv', sep='\t', index_col=0)
    all_ig = pd.read_csv('./matrices/canonical_ig_translocations.tsv', sep='\t', index_col=0)

    return {
        'label': all_label,
        'rna': all_rna,
        'cn_arm': all_cn_arm,
        'cn_seg': all_cn_seg,
        'ig': all_ig,
        'mut': all_mut
    }

def load_dev_data(all_data):
    
    # subset entire dataset to the model development dataset
    # i.e. exclude test set samples
    
    all_label = all_data['label']
    all_rna = all_data['rna']
    all_cn_arm = all_data['cn_arm']
    all_cn_seg = all_data['cn_seg']
    all_mut = all_data['mut']
    all_ig = all_data['ig']

    dev_id = pd.read_csv('./data/splits/dev.txt', sep=',')['PUBLIC_ID']

    dev_label = all_label.loc[dev_id]
    dev_rna = all_rna.loc[dev_id]
    dev_cn_arm = all_cn_arm.loc[dev_id]
    dev_cn_seg = all_cn_seg.loc[dev_id]
    dev_mut = all_mut.loc[dev_id, (all_mut.loc[dev_id].sum(axis=0) > 0)]
    dev_ig = all_ig.loc[dev_id]
    return {
        'label': dev_label,
        'rna': dev_rna,
        'cn_arm': dev_cn_arm,
        'cn_seg': dev_cn_seg,
        'mut': dev_mut,
        'ig': dev_ig
    }

def load_split_data(dev_data, shuffle=10, fold=5):
    # subset model development dataset into 80-20 train and valid subsets
    # shuffle is 1-10 indicating the number of repeats of the 5-fold cross-validation
    # fold is 1-5 as per the 5-fold cross-validation
    # 415/414
    train_id = pd.read_csv(f'./data/splits/{shuffle}/train_{fold}.txt')['PUBLIC_ID']
    # 103/104
    valid_id = pd.read_csv(f'./data/splits/{shuffle}/valid_{fold}.txt')['PUBLIC_ID']
    
    return {
        'train_id': train_id,
        'valid_id': valid_id,
        'train_label': dev_data['label'].loc[train_id],
        'valid_label': dev_data['label'].loc[valid_id],
        'train_rna': dev_data['rna'].loc[train_id],
        'valid_rna': dev_data['rna'].loc[valid_id],
        'train_cn_arm': dev_data['cn_arm'].loc[train_id],
        'valid_cn_arm': dev_data['cn_arm'].loc[valid_id],
        'train_cn_seg': dev_data['cn_seg'].loc[train_id],
        'valid_cn_seg': dev_data['cn_seg'].loc[valid_id],
        'train_ig': dev_data['ig'].loc[train_id],
        'valid_ig': dev_data['ig'].loc[valid_id],
        'train_mut': dev_data['mut'].loc[train_id],
        'valid_mut': dev_data['mut'].loc[valid_id]        
    }

# usage example
# all_data = load_all_data()
# dev_data = load_dev_data(all_data)
# split_data = load_split_data(dev_data, shuffle=10, fold=1)