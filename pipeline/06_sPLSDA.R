library(mixOmics)
source('./utils/datasets.R')
source('./utils/mixo_utils.R')

load_common_data()

X = common.data.rna
Y = common.data.label$risk
tune.splsda.results = tune.splsda(X, Y, ncomp=3,
                                  test.keepX = c(1:3),
                                  validation="Mfold",
                                  folds=5,
                                  nrepeat=10,
                                  cpus=3)

plot(tune.splsda.results, col = color.jet(10))

ncomp=tune.splsda.results$choice.ncomp$ncomp
keepX=tune.splsda.results$choice.keepX
# RNA 8
# comp1  comp2  comp3  comp4  comp5  comp6  comp7  comp8  comp9 comp10 
# 772    772     40     20    772      8      5    600     70      2 
# MUT 2
# comp1  comp2  comp3  comp4  comp5  comp6  comp7  comp8  comp9 comp10 
# 32      3      5      1      1      1     32      2      3      1 
# CN ARM 3
# comp1  comp2  comp3  comp4  comp5  comp6  comp7  comp8  comp9 comp10 
# 20     42     42     42     42      1      3     20      5      3 
# CN SEG 3
# comp1  comp2  comp3  comp4  comp5  comp6  comp7  comp8  comp9 comp10 
# 20    272    272      3    272     50     50     20      2      2 
# IG 2
# comp1 comp2 comp3 
# 3     3     1 


experiment = cv(ncomp=2, keepX=c(3,3))
smry = summary(experiment)
print(smry)
# RNA 772
#     precision.mu precision.sd recall.mu  recall.sd     f1.mu      f1.sd
# SR     0.7806906   0.03356064 0.7491358 0.06121937 0.7628835 0.03528237
# GHR    0.8027771   0.06140978 0.6681410 0.08045259 0.7249106 0.04766826
# FHR    0.1765375   0.05458421 0.3002198 0.10728442 0.2196039 0.06577481

# mut 32
#     precision.mu precision.sd recall.mu  recall.sd     f1.mu      f1.sd
# SR     0.6836178   0.01426483 0.9637037 0.02032585 0.7997323 0.01355117
# GHR    0.8690924   0.08105110 0.2639615 0.05290953 0.4021459 0.06375447
# FHR    0.8118889   0.16187598 0.4409890 0.11045387 0.5620743 0.11017349

# Arm CN 42
#     precision.mu precision.sd recall.mu  recall.sd     f1.mu      f1.sd
# SR     0.7798818   0.03424128 0.5822222 0.08997683 0.6623719 0.06192417
# GHR    0.5527620   0.04737535 0.7442308 0.06577808 0.6328044 0.04496498
# FHR    0.1305104   0.08019212 0.1912088 0.12587781 0.1490906 0.08945834

# Seg CN 272
#     precision.mu precision.sd recall.mu  recall.sd     f1.mu      f1.sd
# SR     0.8015815   0.04596309 0.5600000 0.07469131 0.6559126 0.05415355
# GHR    0.5822356   0.05549459 0.7228846 0.06652352 0.6433229 0.05052709
# FHR    0.1795043   0.05276687 0.3732967 0.12827489 0.2394880 0.06951835

# IG 3
#     precision.mu precision.sd recall.mu  recall.sd     f1.mu      f1.sd
# SR     0.7408388   0.02623810 0.7654321 0.03831669 0.7522020 0.02233259
# GHR    1.0000000   0.00000000 0.5785769 0.07086316 0.7305360 0.05691723
# FHR    0.1763335   0.05794164 0.3634066 0.12833540 0.2364213 0.07786549
