# baseline survival curves by risk group
library(survival)
library(survminer)
library(RColorBrewer)

fhr = read.delim('./annotations/fhr-annotations.tsv',row.names=1)
fhr['risk'] = factor(fhr$risk,levels=c(0,1,2),labels=c("SR","GHR","FHR"))
clin = read.delim('./data/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_STAND_ALONE_SURVIVAL.tsv',row.names=1)
clin2 = read.delim('./data/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_PER_PATIENT.tsv',row.names=1)

data = cbind(clin, risk=fhr[rownames(clin),'risk'])
data = data[!is.na(data['risk']), ] # 884

fit = survfit(Surv(oscdy,censos) ~ risk, data)

# Plot survival curve
line.colors = rev(brewer.pal(3,'Set1'))
png(filename = "./assets/os.png")
ggsurvplot(fit, data, palette=line.colors, risk.table=TRUE, pval=T, pval.method=T, test.for.trend=T, conf.int=T,risk.table.y.text=F)
dev.off()

ggsurvplot(fit, data, palette=line.colors, pval=T, pval.method=T, test.for.trend=T)

# Stratify by transplant or not
transplant = clin2[rownames(data),"bmtx_type"]
transplant = replace(transplant, transplant=="", "None")
transplant = factor(transplant)
dim(data[transplant=="None",]) # 402
dim(data[transplant=="Stem cell, Autologous",]) # 480
dim(data[transplant=="Stem cell, Allogenic",]) # 2

fit.no.transplant = survfit(Surv(oscdy,censos) ~ risk, data[transplant=="None",])
fit.transplant = survfit(Surv(oscdy,censos) ~ risk, data[transplant=="Stem cell, Autologous",])

ggsurvplot(fit.no.transplant,data[transplant=="None",],pval=T,pval.method=T,test.for.trend=T) + ggtitle("No transplant")
ggsurvplot(fit.transplant,data[transplant=="Stem cell, Autologous",],pval=T,pval.method=T,test.for.trend=T) + ggtitle("Transplant")

# Multi-state model
# 0 = alive, 1 = dead from cancer, 2 = dead from other causes
event = clin2[rownames(data),'D_PT_dthreas']
event = replace(event, is.na(event), 0)
event = factor(event,0:2,c('censored','cancer','non-cancer'))

fit2 = survfit(Surv(oscdy,event) ~ 1, data)
plot(fit2,col=line.colors[2:3])
legend("bottomright",fit2$states[2:3],lty=1,lwd=2,col=line.colors[2:3],bty='n')

age = factor(clin2[rownames(data),'D_PT_age']<63, levels=c(T,F), labels=c('Age<63','Age>=63'))
gender = factor(clin2[rownames(data),'D_PT_gender'],levels=1:2,labels=c('Male','Female'))
ISS = factor(clin2[rownames(data),'D_PT_iss'],levels=1:3, labels=c('I','II','III'))

# baseline demographic characteristics
# and their significance
chisq.test(table(age, data$risk))
chisq.test(table(gender, data$risk))
chisq.test(table(ISS, data$risk) )
chisq.test(table(transplant, data$risk))
