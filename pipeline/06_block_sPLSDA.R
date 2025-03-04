library(mixOmics)
source('./utils/datasets.R')
source('./utils/mixo_utils.R')

load_common_data()

X = list(rna=common.data.rna, 
         seg=common.data.cn.seg, 
         arm=common.data.cn.arm)
Y = common.data.label$risk

design = matrix(1, ncol=length(X), nrow=length(X), dimnames=list(names(X),names(X)))
diag(design) = 0

list.keepX = list(rna=c(772),
                  seg=c(20,272),
                  arm=c(20,42))

BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores()-1)

varnames = c(names(list.keepX),'Y')
design = matrix(c(0,1,1,1,1,1,
                  1,0,0,0,0,1,
                  1,0,0,0,0,1,
                  1,0,0,0,0,1,
                  1,0,0,0,0,1,
                  1,1,1,1,1,0),nrow = 6, ncol=6, dimnames = list(varnames,varnames))

design='null'
design='full' 

tune.block.splsda.results = tune.block.splsda(X, Y,
                                              ncomp=c(3),
                                              test.keepX=list.keepX,
                                              design = design,
                                              validation="Mfold",
                                              folds=5,
                                              nrepeat=10,
                                              BPPARAM=BPPARAM)

plot(tune.block.splsda.results, col = color.jet(3))

ncomp = tune.block.splsda.results$choice.ncomp$ncomp
keepX = lapply(tune.block.splsda.results$choice.keepX,function(...)...[1:ncomp])

experiment = cv(ncomp, keepX)
smry = summary(experiment)
print(smry)

# RNA, CNSEG, CNARM, IG (no MUT), ncomp=2, keepX maximum
#     precision.mu precision.sd recall.mu  recall.sd     f1.mu      f1.sd
# SR     0.7405089   0.04793903 0.5012346 0.06499304 0.5955008 0.05652765
# GHR    0.6684943   0.07220454 0.6720513 0.06012949 0.6675045 0.05274122
# FHR    0.1487167   0.04204881 0.4402198 0.15914775 0.2205017 0.06653067

# RNA, CNSEG, CNARM, ncomp=3, Block PLSDA
#     precision.mu precision.sd recall.mu  recall.sd     f1.mu      f1.sd
# SR     0.7774167   0.05509319 0.3765432 0.07090956 0.5030079 0.06936986
# GHR    0.5664309   0.05168621 0.7684487 0.05774509 0.6503679 0.04276787
# FHR    0.1673751   0.05038363 0.5006593 0.13403064 0.2487552 0.06657288

