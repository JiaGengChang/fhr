data = read.delim('./data/SeqFISH Files_MMRF_CoMMpass_IA16a_LongInsert_Canonical_Ig_Translocations.txt',sep='\t')

call.cols = c('SeqWGS_WHSC1_CALL','SeqWGS_CCND3_CALL','SeqWGS_MYC_CALL','SeqWGS_MAFA_CALL','SeqWGS_CCND1_CALL','SeqWGS_CCND2_CALL','SeqWGS_MAF_CALL','SeqWGS_MAFB_CALL')
sample.id = data$Study_Visit_iD
public.id = sub('^(MMRF_[0-9]+)_[1-9]_.*','\\1',sample.id)

res = apply(data[!duplicated(public.id),call.cols],2,as.integer)

new.cols = c('WHSC1','CCND3','MYC','MAFA','CCND1','CCND2','MAF','MAFB')
colnames(res) = new.cols
rownames(res) = public.id[!duplicated(public.id)]

write.table(res,'./matrices/canonical_ig_translocations.tsv',sep='\t',quote=F)
