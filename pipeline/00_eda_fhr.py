import pandas as pd
from scipy.stats import ttest_ind

file_trtresp = '/home/jiageng/Documents/fhr/data/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_STAND_ALONE_TRTRESP.tsv'
df_trtresp = pd.read_csv(file_trtresp, sep='\t', encoding='cp1252')

# Display the first few rows of the DataFrame
for c,v in zip(df_trtresp.columns, df_trtresp.iloc[13, :]):
    print(c, v)
    
    
file_clin = '/home/jiageng/Documents/fhr/data/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_PER_PATIENT.tsv'
df_clin = pd.read_csv(file_clin, sep='\t')

print(df_clin['fhr'].fillna(0).value_counts()) # 359
fhr_patients = df_clin.loc[df_clin['fhr']==1,'PUBLIC_ID']
print(df_clin.fillna(0).groupby(['fhr','D_PT_pdflag'])[['fhr','D_PT_pdflag']].value_counts())

print(df_clin.fillna(0).groupby(['fhr'])['D_PT_ttfpdw'].max())

file_surv = '/home/jiageng/Documents/fhr/data/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_STAND_ALONE_SURVIVAL.tsv'

df_surv = pd.read_csv(file_surv, sep='\t')


# refractory to induction therapy
names_nonresp = df_surv.query('pdflag1 == 1')['PUBLIC_ID'].unique()
names_nonresp.__len__() # 420
# other ways to define refractory to induction therapy
df_trtresp.query('therftrt == 1.0 and line == 1 and fresp == "Partial Response"') # 595
df_surv['censtf1'].value_counts() # 632
df_surv['censfrsp'].value_counts() # 112

# relapse within 12 months.
names_relapse = df_surv.query('ttpfs < 365 and censpfs==1')['PUBLIC_ID'].unique()
names_relapse.__len__() # 189
(df_surv['ttpfs'] < 365).value_counts() # 189
((df_surv['pfscdy'] < 365) & (df_surv['censpfs'] == 1)).value_counts() # 189 
(df_surv['ttfpd'] < 365).value_counts() # 110

df_surv['PUBLIC_ID'].isin(fhr_patients).value_counts() # 359
df_surv_fhr = df_surv.loc[df_surv['PUBLIC_ID'].isin(fhr_patients),:]
df_surv_notfhr = df_surv.loc[~df_surv['PUBLIC_ID'].isin(fhr_patients),:]

pvalues = pd.Series(index=df_surv.columns, dtype=float)

for column in df_surv.columns:
    if df_surv[column].dtype in ['int64', 'float64']:
        t_stat, p_val = ttest_ind(df_surv_fhr[column].dropna(), df_surv_notfhr[column].dropna(), equal_var=False)
        pvalues[column] = p_val

pvalues[~pvalues.isna()].sort_values().head(10)

(df_surv_fhr['pdflag']==0).sum() # 0
(df_surv_fhr['censtf1']==0).sum() # 0
(df_surv_fhr['censrdur']==0).sum() # 0
(df_surv_fhr['censpfs']==0).sum() # 0
(df_surv_fhr['censt2line']==0).sum() # 50
(df_surv_fhr['pdflag1']==0).sum() # 74


fhr_patients_reconstr = \
df_surv.loc[
    ((df_surv['pdflag']==1) &\
    (df_surv['censtf1']==1) &\
    (df_surv['censpfs']==1) &\
    (df_surv['censrdur']==1) &\
    (df_surv['censt2line']==1))
,'PUBLIC_ID']

overlap_patients = set(fhr_patients).intersection(set(fhr_patients_reconstr))
print(f"Number of overlapping patients: {len(overlap_patients)}/{len(fhr_patients_reconstr)} out of {len(fhr_patients)}")


# df_surv_fhr.loc[df_surv_fhr['ttpfs']==1086,:]
# for key, value in df_surv_fhr.loc[df_surv_fhr['ttpfs'] == 1086, :].items():
#     print(key, value.values)
    
    
pvalues_patient = pd.Series(index=df_clin.columns, dtype=float)
df_clin_fhr = df_clin.loc[df_clin['fhr']==1,:]
df_clin_notfhr = df_clin.loc[df_clin['fhr']!=1,:]
for column in df_clin.columns:
    if df_clin[column].dtype in ['int64', 'float64']:
        t_stat, p_val = ttest_ind(df_clin_fhr[column].dropna(), df_clin_notfhr[column].dropna(), equal_var=False)
        pvalues_patient[column] = p_val

pvalues_patient[~pvalues_patient.isna()].sort_values().head(10)


file_clin = '/home/jiageng/Documents/fhr/data/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_PER_PATIENT.tsv'
df_clin = pd.read_csv(file_clin, sep='\t').assign(fhr=lambda x: x['fhr'].fillna(0))
df_clin.groupby(['fhr','D_PT_pdflag'])[['D_PT_pddy','D_PT_ttfpdw']].agg(['min','max','count'])

# refractory to induction therapy
df_trtresp.query('therftrt == 1.0 and line == 1 and fresp == "Partial Response"')['PUBLIC_ID']

# time to progression < 12 months
df_surv.query('pfscdy < 365 and censpfs==1')['PUBLIC_ID']

# translocation-cyclin D classifications
file_tc = '/home/jiageng/Documents/fhr/data/SeqFISH Files_MMRF_CoMMpass_IA16a_RNAseq_Canonical_Ig_Translocations.txt'
df_tc = pd.read_csv(file_tc, sep='\t')\
    .assign(PUBLIC_ID=lambda x: x['Specimen_ID'].str.extract(r'(MMRF_\d{4})'))\
    .filter(regex='PUBLIC_ID|_Call')

# t(4;14) FGFR3 overexpression or t(14;16) c-MAF overexpression
names_tc = df_tc.query('RNASeq_FGFR3_Call == 1 or RNASeq_MAF_Call == 1')['PUBLIC_ID'].unique() 
names_tc.__len__() # 116

# bi-allelic TP53 deletion
df_fish = '/home/jiageng/Documents/fhr/data/SeqFISH Files_MMRF_CoMMpass_IA22_exome_gatk_cna_seqFISH.tsv'
df_fish = pd.read_csv(df_fish, sep='\t').assign(PUBLIC_ID=lambda x: x['SAMPLE'].str.extract(r'(MMRF_\d{4})'))
df_fish.groupby(['SeqExome_Cp_17p13_50percent'])['SeqExome_Cp_17p13'].agg(['min','max'])

import matplotlib.pyplot as plt

# Plot histogram of SeqExome_Cp_17p13
feature = 'SeqExome_Cp_1q21' # SeqExome_Cp_17p13
plt.hist(df_fish[feature].dropna(), bins=30, edgecolor='black')
plt.xlabel(feature)
plt.ylabel('Frequency')
plt.title(f'Histogram of {feature}')
plt.savefig(f'histogram_{feature}.png')
plt.close()


# TP53 mutations
file_mut = '/home/jiageng/Documents/fhr/data/IGV Downloads_MMRF_CoMMpass_IA22_exome_vcfmerger2_IGV_All_Canonical_NS_Variants.mut'
df_mut = pd.read_csv(file_mut, sep='\t')
df_mut = df_mut[df_mut['sample'].str.endswith('BM_CD138pos')].assign(PUBLIC_ID=lambda x: x['sample'].str.extract(r'(MMRF_\d{4})'))

names_mut_tp53 = df_mut.query('chr=="chr17" and GENE=="TP53"')['PUBLIC_ID'].unique()

# TP53 deletion. Either is the same.
names_del_tp53 = df_fish.loc[-0.2 > df_fish['SeqExome_Cp_17p13'],'PUBLIC_ID'].unique() # 101
names_del_tp53 = df_fish.loc[df_fish['SeqExome_Cp_17p13_20percent']==1,'PUBLIC_ID'].unique() # 103

# Deep TP53 deletion
# There are no individuals with -2 > log2FC 

# TP53 LOH
file_loh = '/home/jiageng/Documents/fhr/data/Loss of Heterozygosity Files_MMRF_CoMMpass_IA22_exome_gatk_baf.seg'
df_loh = pd.read_csv(file_loh, sep='\t').assign(PUBLIC_ID=lambda x: x['SAMPLE'].str.extract(r'(MMRF_\d{4})'))
TP53_start=7661779
TP53_end=7687564
LEN=(TP53_end-TP53_start+1)/2
# select where (Start, End) overlaps with (TP53_start, TP53_end) by at least 5000 positions
df_cn_tp53 = df_loh.query("Chromosome=='chr17' and ((Start < @TP53_start and End + @LEN > @TP53_start) or (Start < @TP53_end - @LEN and End >= @TP53_end) or (Start > @TP53_start and End < @TP53_end and (End - Start) >= @LEN))") # 1196

names_loh_tp53 = df_cn_tp53.loc[df_cn_tp53['Segment_Mean'] < 0.25, 'PUBLIC_ID'].unique() # 104

# Biallelic TP53 events
# not sure whether to include deletion + LOH as a bi-allelic event
biallelic_tp53_type1 = set(names_mut_tp53).intersection(set(names_del_tp53))
biallelic_tp53_type2 = set(names_mut_tp53).intersection(set(names_loh_tp53))
# biallelic_tp53_type3 = set(names_del_tp53).intersection(set(names_loh_tp53))
biallelic_tp53 = biallelic_tp53_type1.union(biallelic_tp53_type2) #.union(biallelic_tp53_type3)

print(f"Number of patients with known bi-allelic loss of TP53: {len(biallelic_tp53)}")
df_cn_tp53['PUBLIC_ID'].unique().__len__() # 974
df_fish['PUBLIC_ID'].unique().__len__() # 892
df_mut['PUBLIC_ID'].unique().__len__() # 974
print(len(biallelic_tp53)/974*100)

# 1q gain and ISS3
names_iss3 = df_clin.loc[df_clin['D_PT_iss']==3,'PUBLIC_ID'].unique() # 311
names_gain1q = df_fish.loc[df_fish['SeqExome_Cp_1q21'] > 0.2,'PUBLIC_ID'].unique() # 385
names_gain1q_iss3 = set(names_iss3).intersection(set(names_gain1q)) # 110


# FHR calling pipeline
len(names_relapse) # 189
len(names_nonresp) # 420

fhr_ghr = set(names_relapse).union(set(names_nonresp))
len(fhr_ghr) # 523
ghr = set(names_tc).union(set(names_gain1q_iss3)).union(set(biallelic_tp53)).intersection(fhr_ghr)
len(ghr) # 125
fhr = fhr_ghr.difference(ghr)
len(fhr) # 398

