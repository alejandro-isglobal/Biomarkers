#####################
#retrieve study data#
#####################

library(GEOquery)

gsm <- getGEO("GSE117468",destdir ="./data")

#obtain phenotype data
phenobb <- pData(phenoData(gsm[[1]]))

#obtain clinical data  
id <- rownames(phenobb)
age <- as.numeric(phenobb$"age:ch1")
bmi <- as.numeric(phenobb$"bmi:ch1")
race <- factor(phenobb$"race:ch1")
effbase <- as.numeric(phenobb$"pasi:ch1")
tissue <- phenobb$"tissue:ch1"
t <- factor(factor(phenobb$"treatment:ch1", 
            labels = c("brodalumab","brodalumab", "placebo", NA)), levels=c("placebo", "brodalumab"))
visit <- phenobb$"visit:ch1"
patient <- phenobb$"patientid:ch1"
  
#select subjects with base line, nonlesional tissue and european ancestry data
selBLN <- visit=="BL" & tissue=="non-lesional skin" #& race=="WHITE"
age <- age[selBLN]
bmi <- bmi[selBLN]
t <- t[selBLN]
id <- id[selBLN]
effbase <- effbase[selBLN]
  
pheno <- data.frame(age, bmi, patient=patient[selBLN], t)
rownames(pheno) <- id
  
#obtain PASI at week 12
effend <- eff[visit=="W12" & tissue=="non-lesional skin"] 
names(effend) <- patient[visit=="W12" & tissue=="non-lesional skin"]
effend <- effend[as.character(pheno$patient)]
  
#add effects
pheno <- cbind(pheno, 
               eff = as.factor(effend < effbase), #repose in PASI
               effdif = (effbase-effend)/effbase, #level of repose in PASI
               effbase = effbase, # PASI at baseline
               effend = effend) # PASI at week 12

#store clinical data, store in phenodat
phenodat <- pheno[complete.cases(pheno),]

#obtain annotation data, store in genesIDs
genesIDs <- fData(gsm[[1]])

#obtain transcriptomic data, store in expr
expr <- exprs(gsm[[1]])
expr <- expr[,rownames(phenodat)]
  
genesid <- sapply(strsplit(genesIDs$"Gene Symbol", "/"), function(x) x[1])
names(genesid) <- rownames(genesIDs)
genesentrez <- genesIDs$ENTREZ_GENE_ID
names(genesentrez) <- rownames(genesIDs) 

save(phenodat, expr, genesid, file="GSE117468.Rdata")

#################################################
#Perform transcriptome-wide interaction analysis#
#################################################

library(sva)
library(limma)
#load("GSE117468.Rdata")

##intearaction between treatment and improvement in PASI: t*eff

#compute SVAs
mod0 <- model.matrix( ~  t + eff  + age + bmi, data = phenodat)
mod <- model.matrix( ~ t:eff + t + eff  + age + bmi, data = phenodat)
ns <- num.sv(expr, mod, method="be")
ss <- sva(expr, mod, mod0, n.sv=ns)$sv
modss <- cbind(mod, ss)
  
#estimate associations
fit <- lmFit(expr, modss)
fit <- eBayes(fit)
  
#volcano plot
pdf("volcano.pdf")
 volcanoplot(fit, highlight=11, coef="tbrodalumab:effTRUE", 
            names=genesid[rownames(fit$coefficients)], cex=0.1)
dev.off()

tt <- topTable(fit, number=Inf, coef="tbrodalumab:effTRUE")

#Select significant associations
trascriptname <- rownames(tt)
sigGenespso <- trascriptname[tt$adj.P.Val<0.05]

tt <- data.frame(Gene= genesid[sigGenespso], tt[sigGenespso,])

tt[,c(2:4,7)] <- format(tt[,c(2:4,7)], digits=3) 
tt[,5:6] <- format(tt[,5:6], digits=3, scientific=TRUE) 

write.table(tt, file="tableS1.txt", quote=FALSE, sep="\t")

toptrans <- rownames(tt)[c(1,2,4)]

################
#plot top genes#
################

library(vioplot)


pdf("genes.pdf")
par(mfrow=c(1,3))

top <- rownames(tt)[1]
tr <- log(expr[top,])
res <- summary(lm(tr~modss[,-c(1,2,3,6)]))$residuals
et <- modss[,6]
fc <- factor(modss[,6], labels=c("\n Placebo \n OR \n No response", "\n Brodalumab  \n AND \n Response"))
vioplot(res~fc, xlab="", main="NR4A2", ylab="log-fold change residuals", cex.names=0.5)

top <- rownames(tt)[2]
tr <- log(expr[top,])
res <- summary(lm(tr~modss[,-c(1,2,3,6)]))$residuals
et <- modss[,6]
fc <- factor(modss[,6], labels=c("\n Placebo \n OR \n No response", "\n Brodalumab  \n AND \n Response"))
vioplot(res~fc, xlab="", main="IGH", ylab="log-fold change residuals", cex.names=0.5)


top <- rownames(tt)[8]
tr <- log(expr[top,])
res <- summary(lm(tr~modss[,-c(1,2,3,6)]))$residuals
et <- modss[,6]
fc <- factor(modss[,6], labels=c("\n Placebo \n OR \n No response", "\n Brodalumab  \n AND \n Response"))
vioplot(res~fc, xlab="", main="EGR4 ", ylab="log-fold change residuals", cex.names=0.5)


dev.off()

############
#enrichment#
############
library(clusterProfiler)

mappedgenesIds <- genesIDs[rownames(tt), "ENTREZ_GENE_ID"]
mappedgenesIds <- unique(unlist(strsplit(mappedgenesIds, " /// ")))

#run enrichment in GO
GO <- enrichGO(gene = mappedgenesIds, 'org.Hs.eg.db', ont="MF", pvalueCutoff=0.05, pAdjustMethod="BH")


GO <- data.frame(ID=GO$ID, Description=GO$Description, 
                 Padj=format(GO$p.adjust, digits=3, sientific=TRUE), GeneRatio=GO$GeneRatio)

write.table(GO, file="enriched.txt",
            quote=FALSE, sep="\t", row=FALSE)

GO

########################
#Random causal modeling#
########################
library(teff)

#Prepare data, features: trascription data, teff: treatment, effect and covariates
teffdata <- modss[,-c(1,6)]
colnames(teffdata)[1:2] <- c("t", "eff") 
colnames(teffdata)[5:ncol(teffdata)] <- paste0("cov",5:ncol(teffdata))

psoriasis <- list(features=t(expr), teffdata=teffdata)

#save(psoriasis, sigGenespso, file="psoriasis.RData")

#grow forest in 80% of individuals
#estimate individual treatment effect in the 20% left

#load("psoriasis.RData")

pso <- teff::profile(psoriasis, featuresinf=sigGenespso, seed=1234, dup=TRUE)

# compare estimates of treatment effects on PASI improvement to
# observed levels of PASI improvement
d <- phenodat[pso$subsids,"effdif"]
p <- pso$predictions
t <- pso$treatment+1

summary(lm(log(p/(1-p))~t))


library(drc)

dd <- d[t==2]
dp <- p[t==2]
metB<-drm(dd*100~dp, fct=LL.3())

dd <- d[t==1]
dp <- p[t==1]
metP<-drm(dd*100~dp, fct=LL.3())



pdf("pred.pdf")

plot(metB, log = "", pch=16, col="red", ylim=c(-100,100), xlim=c(0.2,0.45), ylab="% of PASI improvement at week 12",
     xlab="Probability of response to brodalumab at baseline")

plot(metP, log = "", pch=16, col="black", ylim=c(-100,100), xlim=c(0.2,0.45), add=TRUE)

points(p[t==1], d[t==1]*100, pch=16)
legend("bottomleft", legend=c("Brodalumab", "Placebo"), 
       pch=16, col=c("red","black"), bty="n")

dev.off()

noEffect(metB)
noEffect(metP)



####
d <- phenodat$effdif
p <- phenodat$effbase 
t <- as.numeric(phenodat$t)

summary(lm(p~t))


library(drc)

dd <- d[t==2]
dp <- p[t==2]
metB<-drm(dd*100~dp, fct=LL.3())

dd <- d[t==1]
dp <- p[t==1]
metP<-drm(dd*100~dp, fct=LL.3())



pdf("predBase.pdf")

plot(metB, log = "", pch=16, col="red", ylim=c(-100,100), xlim=c(0,50), ylab="% of PASI improvement at week 12",
     xlab="PASI at baseline")

plot(metP, log = "", pch=16, col="black", ylim=c(-100,100), xlim=c(0,50), add=TRUE)

points(p[t==1], d[t==1]*100, pch=16)
legend("bottomleft", legend=c("Brodalumab", "Placebo"), 
       pch=16, col=c("red","black"), bty="n")

dev.off()

noEffect(metB)
noEffect(metP)



plot(phenodat$effdif, phenodat$effbase, col= 1*(phenodat$t=="brodalumab")+1, pch=16)
     
cor.test(phenodat$effdif[phenodat$t=="brodalumab"], phenodat$effbase[phenodat$t=="brodalumab"])

cor.test(phenodat$effdif[phenodat$t=="placebo"], phenodat$effbase[phenodat$t=="placebo"])

, phenodat$effbase