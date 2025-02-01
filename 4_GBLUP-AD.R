# Estimate genetic parameters using a Aditive Dominant  GBLUP model -------------------------
 
rm(list=ls())
gc()

# Library -----------------------------------------------------------------

packages <- c("tidyverse", "asreml", "lme4", "car", "snpReady", "AGHmatrix",
              "ASRgenomics", "redxl", "data.table", "gridExtra", "PupillometryR",
              "readxl")

# Use pacman to load and install packages if needed
pacman::p_load(char = packages)

# Phenotypic data ---------------------------------------------------------

pheno <- read_excel("data_raw/Dados_phen_3anos.xlsx")
head(pheno); dim(pheno)
str(pheno)

pheno[1:19] = lapply(pheno[1:19], as.factor) #transform in factors
str(pheno)

names(pheno)
pheno = pheno %>% select(IND, FAM, Mae, Pai, bloco, dap, h_total, secao_ind, vol_ind, IMA)
names(pheno) = c("IND", "FAM", "MOTHER", "FATHER", "REP", "DAP", "HT", "SEC_IND", "VOL", "IMA")

# Relationship Matrix -----------------------------------------------------------------
G.inv = read_rds("data/matrices/G.inv.rds")
D.inv = read_rds("data/matrices/D.inv.rds")

# Modelo Aditivo ----------------------------------------------------------
pheno$GENA = as.factor(pheno$IND)

asreml.options(ai.sing=TRUE, threads=-1)
ma <- asreml(
  fixed = DAP ~ REP,
  random = ~ vm(GENA, G.inv), 
  na.action = na.method(x = "exclude", y = "exclude"), 
  workspace = "4gb",
  data = pheno,
  ai.sing=TRUE
)
summary(ma)$varcomp

asreml.options(ai.sing=TRUE, threads=-1)
ma2 <- asreml(
  fixed = DAP ~ REP,
  random = ~ vm(GENA, G.inv) + FAM, 
  na.action = na.method(x = "exclude", y = "exclude"), 
  workspace = "4gb",
  data = pheno,
  ai.sing=TRUE
)

summary(ma2)$varcomp

# Modelo Aditivo Dominante ----------------------------------------------------------
library(asreml)

pheno$GEND = as.factor(pheno$IND)

asreml.options(ai.sing=TRUE, threads=-1)
md <- asreml(
  fixed = DAP ~ REP,
  random = ~ vm(GENA, G.inv) +     # aditivo
    vm(GEND, D.inv),               # dominante
  na.action = na.method(x = "exclude", y = "exclude"), 
  workspace = "4gb",
  data = pheno,
  ai.sing=TRUE
)

summary(md)$varcomp


asreml.options(ai.sing=TRUE, threads=-1)
md2 <- asreml(
  fixed = DAP ~ REP,
  random = ~ vm(GENA, G.inv) +     # aditivo
             vm(GEND, D.inv) +     # dominante
             FAM,                  # family
  na.action = na.method(x = "exclude", y = "exclude"), 
  workspace = "4gb",
  data = pheno,
  ai.sing=TRUE
)

summary(md)$varcomp
summary(md2)$varcomp
# AIC ---------------------------------------------------------------------
library(asremlPlus)
mods.asr <- list(ma, ma2, md, md2)
aic = infoCriteria(mods.asr, IClikelihood = "REML"); aic

# fixedDF varDF NBound      AIC      BIC    loglik
# 1       0     2      0 2790.049 2799.182 -1393.024
# 2       0     3      0 2653.663 2667.363 -1323.831
# 3       0     3      0 2775.948 2789.648 -1384.974
# 4       0     4      0 2626.443 2644.710 -1309.222

# Results md -----------------------------------------------------------------

# BLUES e BLUP ------------------------------------------------------------

#BLUEs
BLUE = data.frame(summary(md, coef=TRUE)$coef.fixed); dim(BLUE)

write.table(BLUE, "results/BLUES-AD-GBLUP.txt")

# # Efeito de familia
# BLUP = data.frame(summary(md, coef=TRUE)$coef.random); dim(BLUP)
# pattern = "vm\\(GEND, D\\.inv\\)_|vm\\(GENA, G\\.inv\\)_"
# indices_fam <- grep(pattern, rownames(BLUP))
# length(indices_fam)
# fam_eff <- BLUP[-indices_fam, ]
#write.csv(fam_eff, "results/fam_eff_md.csv")

# Efeito aditivo
BLUP = data.frame(summary(md, coef=TRUE)$coef.random); dim(BLUP)
pattern <- "vm\\(GENA"
indices_a <- grep(pattern, rownames(BLUP))
length(indices_a)
a <- BLUP[indices_a, ] ; dim(a)
rownames(a) = gsub("^vm\\(GENA, G\\.inv\\)_", "", rownames(a))

# Efeito dominancia
BLUP = data.frame(summary(md, coef=TRUE)$coef.random); dim(BLUP)
pattern <- "vm\\(GEND"
indices_d <- grep(pattern, rownames(BLUP))
length(indices_d)
d <- BLUP[indices_d, ] ; dim(d)
rownames(d) = gsub("^vm\\(GEND, D\\.inv\\)_", "", rownames(d))

# eBLUP ------------------------------------------------------------

eBLUP = a %>%
  rownames_to_column("rowname") %>%
  inner_join(d %>% rownames_to_column("rowname"), by = "rowname") %>%
  column_to_rownames("rowname") ; dim(eBLUP)

eBLUP = eBLUP[ ,-c(3,6)]

colnames(eBLUP) = c("a", "se_a","d", "se_d")
eBLUP$g = eBLUP$a + eBLUP$d
head(eBLUP)

write.table(eBLUP, "results/BLUPS-AD-GBLUP.txt")

# Variance Components -----------------------------------------------------
vc.md = cbind.data.frame(
  summary(md)$varcomp)[ ,c(1,2)] ; vc.md

vc.md[ ,c(1:2)] = sapply(vc.md[ ,c(1:2)], round, digit =2)
write.table(vc.md, "results/varcomp_AD-GBLUP.txt")

# Herdabilidade -----------------------------------------------------------
h2a = asreml::vpredict(md,h2a~((V1)/( V1 + V2 + V3))); h2a
d2 = asreml::vpredict(md,d2~((V2)/( V1 + V2 + V3))); d2
H2 = asreml::vpredict(md,H2~((V1 + V2)/(V1 + V2 + V3))); H2

# Additive: Accuracy and Reliability ------------------------------

Va <- vc.md[1,1]; Va
eBLUP$PEV_a <- eBLUP$se_a^2
eBLUP$acc_a <- sqrt(1 - (eBLUP$PEV_a)/Va)
eBLUP$rel_a <- eBLUP$acc_a^2
mean_acc.a = mean(eBLUP$acc_a, na.rm = TRUE); mean_acc.a
mean_real.a = mean(eBLUP$rel_a, na.rm = TRUE); mean_real.a

#Accuracy and Reliability - Dominance
Vd <- vc.md[2,1]; Vd
eBLUP$PEV_d <- eBLUP$se_d^2
eBLUP$acc_d <- sqrt(1 - eBLUP$PEV_d/Vd)
eBLUP$rel_d <- eBLUP$acc_a^2
mean_acc.d = mean(eBLUP$acc_d, na.rm = TRUE); mean_acc.d
mean_real.d = mean(eBLUP$rel_d, na.rm = TRUE); mean_real.d

# Ganho -------------------------------------------------------------------
str(eBLUP)


sel_g = eBLUP[order(eBLUP$g, decreasing = T),][1:30, ] # 1% dos selecionados equivale a aproximadamente 10% dos individuos no experimento
sel_g$Rank = seq(1:nrow(sel_g))

gain.g = (mean(sel_g$g))/
  mean(pheno$DAP, na.rm=T) * 100 ; gain.g

# Encontra a familia, parcela, etc com base no pheno
sel_gen = pheno[pheno$GENA  %in% rownames(sel_g) , ]

head(sel_gen); dim(sel_gen)
ordering <- factor(sel_gen$GENA, levels = rownames(sel_g))
selected_gen <- sel_gen[order(ordering), , drop = FALSE]

write.table(selected_gen, "results/selected_gen_AD-GBLUP.txt")

mean_sel = mean(selected_gen$DAP,na.rm = TRUE); mean_sel
mean_zero = mean(as.numeric(pheno$DAP), na.rm=TRUE); mean_zero

herits = rbind(h2a, d2, H2); herits
write.table(herits, "results/heritabilities_AD-GBLUP.txt")

acc = rbind(mean_acc.a, mean_real.a, mean_acc.d, mean_real.d, gain.g, mean_sel, mean_zero); acc
write.table(acc, "results/accuracies_AD-GBLUP.txt")

# save image --------------------------------------------------------------

save.image("save_images/images_AD-GBLUP.RData")


# FAM MALE effect --------------------------------------------------------------

str(pheno)
View(pheno)
asreml.options(ai.sing=TRUE, threads=-1)
md_father <- asreml(
  fixed = DAP ~ REP,
  random = ~ vm(GENA, G.inv) +     # aditivo
    vm(GEND, D.inv) +     # dominante
    FATHER,                  # family
  na.action = na.method(x = "exclude", y = "exclude"), 
  workspace = "4gb",
  data = pheno,
  ai.sing=TRUE
)

summary(md_father)$varcomp



#BLUEs
BLUE = data.frame(summary(md_father, coef=TRUE)$coef.fixed); dim(BLUE)

#write.table(BLUE, "results/BLUES-AD-GBLUP.txt")

# Efeito de familia
BLUP = data.frame(summary(md_father, coef=TRUE)$coef.random); dim(BLUP)
pattern = "vm\\(GEND, D\\.inv\\)_|vm\\(GENA, G\\.inv\\)_"
indices_fam <- grep(pattern, rownames(BLUP))
length(indices_fam)
fam_eff_father <- BLUP[-indices_fam, ]
write.csv(fam_eff_father, "results/fam_eff_father.csv")

# # Efeito aditivo
# BLUP = data.frame(summary(md_father, coef=TRUE)$coef.random); dim(BLUP)
# pattern <- "vm\\(GENA"
# indices_a <- grep(pattern, rownames(BLUP))
# length(indices_a)
# a <- BLUP[indices_a, ] ; dim(a)
# rownames(a) = gsub("^vm\\(GENA, G\\.inv\\)_", "", rownames(a))
# 
# # Efeito dominancia
# BLUP = data.frame(summary(md_father, coef=TRUE)$coef.random); dim(BLUP)
# pattern <- "vm\\(GEND"
# indices_d <- grep(pattern, rownames(BLUP))
# length(indices_d)
# d <- BLUP[indices_d, ] ; dim(d)
# rownames(d) = gsub("^vm\\(GEND, D\\.inv\\)_", "", rownames(d))
# 
# # eBLUP ------------------------------------------------------------
# 
# eBLUP = a %>%
#   rownames_to_column("rowname") %>%
#   inner_join(d %>% rownames_to_column("rowname"), by = "rowname") %>%
#   column_to_rownames("rowname") ; dim(eBLUP)
# 
# eBLUP = eBLUP[ ,-c(3,6)]
# 
# colnames(eBLUP) = c("a", "se_a","d", "se_d")
# eBLUP$g = eBLUP$a + eBLUP$d
# head(eBLUP)
# 
# write.table(eBLUP, "results/BLUP-md-father.txt")

