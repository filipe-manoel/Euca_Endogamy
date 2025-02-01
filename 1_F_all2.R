rm(list=ls())
gc()

# Packages ----------------------------------------------------------------

library(tidyverse)

# Making PLINK files ------------------------------------------------------

# PED ---------------------------------------------------------------------

snp_data = read.csv("data_raw/snps_todos_ind.csv", sep = ";")
snp_data[1:5, 1:5]
dim(snp_data)

rownames(snp_data) = snp_data[ ,1]
snp_data = snp_data[ ,-1]

# Corrected pedigree
pedc= data.frame(readxl::read_excel("data_raw/pedigree_corrigido.xlsx", sheet = "pedC"))
head(pedc)
dim(pedc)


diff_inds = setdiff(as.character(rownames(snp_data)), as.character(pedc$IND))


FAM = rep(0, length(diff_inds))
IND = diff_inds
Father = rep(0,length(diff_inds))
Mother = rep(0,length(diff_inds))
SEX = rep(0,length(diff_inds))
FEN = rep(-9, length(diff_inds))

parents= data.frame(cbind(FAM, IND, Father, Mother, SEX, FEN))
pedc = rbind(parents, pedc)
dim(pedc)
dim(snp_data)

pedc= pedc[as.character(pedc$IND) %in% as.character(rownames(snp_data)), ]

# Remove duplicates based on the 'ID' column
pedc = pedc[!duplicated(pedc$IND), ]
dim(pedc)
pedc = pedc[order(pedc$IND), ]

snp_data = snp_data[order(rownames(snp_data)), ]

as.character(rownames(snp_data))==as.character(pedc$IND)


# Function to split alleles for each SNP column
split_alleles = function(df) {
  new_df = df  # Create a copy to add columns without altering original data
  for (snp in colnames(df)) {
    alleles = do.call(rbind, strsplit(df[[snp]], ""))  # Split each genotype
    new_df[[paste0(snp, "_A1")]] = alleles[, 1]        # First allele
    new_df[[paste0(snp, "_A2")]] = alleles[, 2]        # Second allele
    new_df[[snp]] = NULL                               # Remove original column if desired
  }
  return(new_df)
}

# Apply the function
snp_data_split= data.frame(split_alleles(snp_data))
snp_data_split[1:5, 1:15]

ped_file = cbind(pedc, snp_data_split)
ped_file[1:5, 1:15]
dim(ped_file)

write.table(ped_file, "data/data_plink/snp_all.ped", quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep = " ")


# Remove the outbred ------------------------------------------------------
# 
# genotypes = data.frame(readxl::read_excel("data_raw/pedigree_corrigido.xlsx", sheet = "inbread"))
# head(genotypes)
# 
# inbread_genotypes =  unique(genotypes$IND[genotypes$teste == "auto"])
# length(inbread_genotypes)
# 
# snp_inbread = ped_file[ped_file$IND %in% inbread_genotypes, ]
# dim(snp_inbread)
# 
# # Optional: Save to file if needed
# write.table(snp_inbread, file = "data/data_plink/snp_inbread.ped", quote = FALSE,
#             row.names = FALSE, col.names = FALSE, sep = " ")

# MAP ---------------------------------------------------------------------

map_file = read.csv("data_raw/Map_snps.csv", sep = ";")
colnames(map_file) = c("Chromossome", "SNP", "Gen_dist", "BP_position")
head(map_file)
dim(map_file)

snp_names = colnames(snp_data)

# all(as.character(map_file$SNP)) == all(as.character(snp_names))
# all
# snp_names

map_file = map_file[match(snp_names, map_file$SNP), ]
dim(map_file)


# Write to .map file
write.table(map_file, "data/data_plink/snp_all.map", quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep = " ")


# Convert .ped and .map to Binary Files in PLINK --------------------------

#Function
runPLINK <- function(PLINKoptions = "") system(paste("utilities/plink.exe", PLINKoptions))


# Run PLINK command in R
runPLINK("--file data/data_plink/snp_all --make-bed --out data/data_plink/snp_all_binary")


# Quality control with PLINK ---------------------------------------------------


## --geno 0.1: Excludes SNPs with more than 10% missing data.
## --maf 0.05: Excludes SNPs with a minor allele frequency less than 0.05.
## --hwe 0.001: Excludes SNPs that do not conform to Hardy-Weinberg equilibrium (HWE) at a p-value threshold of 0.001.


runPLINK("--bfile data/data_plink/snp_all_binary --geno 0.1 --maf 0.05 --make-bed --out data/data_plink/snp_all_qc")


# RzooRoH - Runs of Homozigosity -----------------------------------------------------------------

# Convert from ped format to Oxford GEN format

runPLINK("--bfile data/data_plink/snp_all_qc --recode oxford --out data/data_plink/data_all_ROH")

#1- Reading the data

library(RZooRoH)
library(data.table)

gen_file <- RZooRoH::zoodata(genofile = "data/data_plink/data_all_ROH.gen", zformat = "gp")
head(gen_file@genos[,1:6])

#nome dos individuos
names = read.table("data/data_plink/data_all_ROH.sample")
ind_names = names[3:nrow(names), 2]
length(ind_names)

#nome das familias
fam_names = names[3:nrow(names), 1]
length(fam_names)

#concatene ind e fam

fam_ind <- paste(fam_names, ind_names, sep = "_")

gen_file@sample_ids = fam_ind



#2- Reading the model

# Defining 10 classses of Haplotypes by descendence (HBD)
# mod1R <- RZooRoH::zoomodel(predefined=FALSE,K=2,krates=c(10,10),layers=FALSE)
# mod1R

mod10R <- zoomodel()
mod10R


# 3 - Analysis ----------------------------------------------------------------

#To obtain the locus specific HBD probabilities for all individuals

library(parallel)   # Base package for parallel computing

# Detect cores and set up cluster
n_cores = detectCores() - 1

inbreading_results <- zoorun(mod10R, gen_file, localhbd = TRUE, nT= n_cores) #nT paralelizaa função

inbreading_results@ids = fam_ind
inbreading_results@sampleids = fam_ind

# save.image("save_images/image_all_4_12.RData")
load("save_images/image_all_4_12.RData")


# Parameters --------------------------------------------------------------

#Quanto cada classe contribui para cada indíduo


# Inbreeding --------------------------------------------------------------


# To obtain an inbreeding coefficient including all HBD classes, 
# autozygosity must be summed over all HBD classes or can be obtained
# as 1 minus the non-HBD proportion (which is the 10th column).

x <- 1- inbreading_results@realized[,10]
hist(x,nc=20,main="",xlab="Inbreeding coefficient",col='coral2')

#ate o T=64
x <- t(apply(inbreading_results@realized[,3:6],1,cumsum))
hist(x[,4],nc=20,main="",xlab="Inbreeding coefficient (T = 64)",col='tomato')



# Summary -----------------------------------------------------------------

summary(realized(inbreading_results))

summary(cumhbd(inbreading_results,100))


allseg <- rohbd(zres = inbreading_results)
summary(allseg$length)

hist(allseg$length/1000000,xlab="Length of HBD segment (in cM)",main="",col='tomato',nc=100)


dtw <- zoodoris(zres = inbreading_results, minv = 5, maxv = 25, glen = 2450, method = "sharing")

matplot(dtw$win_st,cbind(dtw$sharing),
        xlab = 'Length of HBD segment in cM', ylab = 'Proportion of the genome',
        type ='p', pch=c(21),bg=c('royal blue'),
        col=c('black'))


zooplot_prophbd(list(Soay = inbreading_results), cols = 'tomato', style = 'boxplot')



zooplot_prophbd(list(Soay=inbreading_results),style='lines', cumulative = TRUE)

pop2 <- list(Soay=inbreading_results)
zooplot_individuals(pop2, cumulative = TRUE)

zooplot_individuals(pop2, cumulative = TRUE,
                    toplot = toplot)


#> Warning in zooplot_individuals(pop2, cumulative = TRUE):
#> All individuals are plotted; use toplot to select individuals



zooplot_partitioning(list(Wiltshire=inbreading_results), ylim = c(0,0.7), nonhbd = FALSE)


# teste -------------------------------------------------------------------

a = data.frame(inbreading_results@realized[ ,1:9])
dim(a)
a$soma = rowSums(a)
a$IND = inbreading_results@ids

a_order = a[order(a$soma, decreasing = T), ]

high_F_ind = a_order[1:20, ]

toplot = list(high_F_ind$IND)
class(toplot)


#svglite("outputs/HBD_HighF.svg")
zooplot_partitioning(list(HighF=inbreading_results),
                     toplot = toplot, plotids= TRUE,
                     nonhbd = FALSE,
                     vertical = TRUE)
#dev.off()


#svglite("outputs/F_cumulative_values.svg")
zooplot_prophbd(list(HighF=inbreading_results),style='lines', cumulative = TRUE)
#dev.off()


# Grafico cromossomo 5

y1 <- probhbd(zres = inbreading_results, zooin = gen_file, id = 18, chrom = 5,
              startPos = 0, endPos = 5)

x <-gen_file@bp[gen_file@chrbound[5,1]:gen_file@chrbound[5,2]]
x <- x[x >=0 & x<= 74630338]/1000000

plot(y1~x,type='b',ylim=c(0,1),ylab='HBD probability',col='red',
     xlab='Position on chr5 (Mb)')

y2 <- probhbd(zres = inbreading_results, zooin = gen_file, id = 2, chrom = 5,
              startPos = 0, endPos = 74630338)
par(new=TRUE)
plot(y1~x,type='b',ylim=c(0,1),ylab='',col='royalblue',xlab='',axes=FALSE)
par(new=TRUE)
plot(y2~x,type='b',ylim=c(0,1),ylab='',col='orange',xlab='',axes=FALSE)


# F - Inbreeding coefficient --------------------------------------------------------------


zooplot_individuals(inbreading_results, cumulative = TRUE, 
                    toplot = toplot)


F_roh = cumhbd(inbreading_results, T = NULL) #An array with the compute inbreeding coefficients for all the individuals in the analysis.
F_roh = data.frame(F_roh)
#F_roh$fator2 = stringr::str_extract(inbreading_results@ids, "9[23].*")
F_roh$IND =inbreading_results@sampleids
mean(F_roh$F_roh)
head(F_roh)
dim(F_roh)

# clean the data

F_roh$fator3 = stringr::str_extract(F_roh$IND, "_[92].*")
F_roh$fator4 = sapply(strsplit(F_roh$IND, "_"), `[`, 2)

# Criando a terceira coluna
F_roh = F_roh %>%
  mutate(
    fator5 = apply(
      cbind(fator3, fator4), 
      1, 
      function(row) {
        words = na.omit(row) # Remove os NAs
        words = words[nchar(words) > 3] # Filtra palavras com mais de 3 letras
        paste(words, collapse = " ") # Junta as palavras restantes
      }
    )
  )


F_roh$fator6 = sapply(strsplit(F_roh$fator5, " "), 
                      function(x) if (length(x) > 1) x[2] else x[1])

F_roh$INDIVIDUAL = gsub("_", "", F_roh$fator6)


F_roh$FAMILY = sapply(strsplit(F_roh$IND, "_"), `[`, 1)
F_roh$FAMILY = ifelse(F_roh$FAMILY==0, F_roh$INDIVIDUAL, F_roh$FAMILY)

F_roh$FAMILY = gsub("1-50-VT04", "VT04", F_roh$FAMILY)


F_roh = F_roh %>% select(FAMILY, INDIVIDUAL, F_roh)

write.table(F_roh, "outputs/F_roh_all.txt")


# Mean F_roh per family



# Regressions 3 years -----------------------------------------------------


# Phenotypic data 3 years
pheno = readxl::read_excel("data_raw/Dados_phen_3anos.xlsx")
str(pheno)
dim(pheno)

pheno_roh = pheno[as.character(pheno$IND)%in%as.character(F_roh$IND), ]
dim(pheno_roh)

F_roh = F_roh %>% rename(IND = INDIVIDUAL)
df1 = merge(pheno_roh, F_roh, by="IND")
names(df1)

df1 = df1 %>% select(IND, FAMILY, Mae, Pai, CLASSE, F_roh, dap, h_total, secao_ind, vol_ind, IMA)
df1[1:5] = lapply(df1[1:5], as.factor) #transform in factors
levels(df1$FAMILY)

write.table(df1, "outputs/pheno_data_F.txt")

# Family - Regression DAP in F by Family -------------------------------------------

# Calculate R-squared for each family
r2_data = df1 %>%
  group_by(FAMILY) %>%
  summarise(
    R2 = summary(lm(dap ~ F_roh))$r.squared
  )

mean_endogamy = df1 %>%
  group_by(FAMILY) %>%
  summarise(ME = mean(F_roh))

mean(mean_endogamy$ME)
#F= 0.2040474

# Reorder the Family factor based on F (descending order)
mean_endogamy$FAMILY = as.factor(mean_endogamy$FAMILY)
mean_endogamy_order = mean_endogamy[order(mean_endogamy$ME, decreasing = T), ]

levels(df1$FAMILY) == levels(mean_endogamy_order$FAMILY)

# Ensure df2$COLUMN has the desired order as a factor
mean_endogamy_order$FAMILY = factor(mean_endogamy_order$FAMILY, levels = unique(mean_endogamy_order$FAMILY))

write.table(mean_endogamy_order, "outputs/F_by_family.txt")

# Reorder df1 based on the factor levels of df2$COLUMN
df1_ordered = df1 %>%
  mutate(FAMILY = factor(FAMILY, levels = levels(mean_endogamy_order$FAMILY))) %>%
  arrange(FAMILY)


# Create the faceted plot by family
ggplot(df1_ordered, aes(x = F_roh, y = dap)) +
  geom_point(color = "#069C56", alpha = 0.7, size=3) +            # Scatter plot
  geom_smooth(method = "lm", color = "#FF681E", se = TRUE, 
              fill = "darkgray", alpha = 0.2) +         #Regression line    
  xlim(0, 0.6)+
  ylim(-5, 20)+
  facet_wrap(~ FAMILY) +                               # Facets sorted by R2
  geom_text(data = mean_endogamy, 
            aes(x = 0.5, y = -2, label = paste0("F = ", round(ME, 2))), 
            inherit.aes = FALSE, size = 2, color = "black") +  # Add Mean Endogamy
  geom_text(data = r2_data, 
            aes(x = 0.5, y = -4, label = paste0("R² = ", round(R2, 2))), 
            inherit.aes = FALSE, size = 2, color = "black") +  # Add R2
  
  labs(
    title = "Regression of F_ROH on DAP by Family",
    x = "F_ROH",
    y = "DAP"
  ) +
  theme_minimal()


# Male - Regression DAP in F by Male --------------------------------------

# Calculate R-squared for each family
r2_data = df1 %>%
  group_by(Pai) %>%
  summarise(
    R2 = summary(lm(dap ~ F_roh))$r.squared
  )

mean_endogamy = df1 %>%
  group_by(Pai) %>%
  summarise(ME = mean(F_roh))

mean(mean_endogamy$ME)
#F = 0.2366702

# Reorder the Family factor based on F (descending order)
mean_endogamy$Pai = as.factor(mean_endogamy$Pai)
mean_endogamy_order = mean_endogamy[order(mean_endogamy$ME, decreasing = T), ]

levels(df1$Pai) == levels(mean_endogamy_order$Pai)

# Ensure df2$COLUMN has the desired order as a factor
mean_endogamy_order$Pai = factor(mean_endogamy_order$Pai,
                                 levels = unique(mean_endogamy_order$Pai))

write.table(mean_endogamy_order, "outputs/Male_F_by_Male.txt")

# Reorder df1 based on the factor levels of df2$COLUMN
df1_ordered = df1 %>%
  mutate(Pai = factor(Pai, levels = levels(mean_endogamy_order$Pai))) %>%
  arrange(Pai)


# Create the faceted plot by family
ggplot(df1_ordered, aes(x = F_roh, y = dap)) +
  geom_point(color = "#069C56", alpha = 0.7, size=3) +            # Scatter plot
  geom_smooth(method = "lm", color = "#FF681E", se = TRUE, 
              fill = "darkgray", alpha = 0.2) +         #Regression line    
  xlim(0, 0.6)+
  ylim(-5, 20)+
  facet_wrap(~ Pai) +                               # Facets sorted by R2
  geom_text(data = mean_endogamy, 
            aes(x = 0.5, y = -1, label = paste0("F = ", round(ME, 2))), 
            inherit.aes = FALSE, size = 3, color = "black") +  # Add Mean Endogamy
  geom_text(data = r2_data, 
            aes(x = 0.5, y = -4, label = paste0("R² = ", round(R2, 2))), 
            inherit.aes = FALSE, size = 3, color = "black") +  # Add R2
  labs(
    title = "Regression of F_ROH on DAP by Pai",
    x = "F_ROH",
    y = "DAP"
  ) +
  theme_minimal()



# Individual - Regression DAP in F by Individual -------------------------------------------


# Calculate R-squared for each family
r2_ind = df1 %>%
  summarise(
    R2 = summary(lm(dap ~ F_roh))$r.squared
  )


mean_endogamy_ind = df1 %>%
  summarise(ME = mean(F_roh))



ggplot(df1_ordered, aes(x = F_roh, y = dap)) +
  geom_point(color = "#069C56", alpha = 0.7, size=3) +            # Scatter plot
  geom_smooth(method = "lm", color = "#FF681E", se = TRUE, 
              fill = "darkgray", alpha = 0.2) +         #Regression line    
  geom_text(data = mean_endogamy_ind, 
            aes(x = 0.6, y = 20, label = paste0("F = ", round(ME, 2))), 
            inherit.aes = FALSE, size = 4, color = "black") +  # Add Mean Endogamy
  geom_text(data = r2_ind, 
            aes(x = 0.6, y = 18, label = paste0("R² = ", round(R2, 2))), 
            inherit.aes = FALSE, size = 4, color = "black") +  # Add R2
  
  labs(
    title = "Regression of F_ROH on DAP",
    x = "F_ROH",
    y = "DAP"
  ) +
  theme_minimal()





# 3 IND per Family --------------------------------------------------------


sel_fam = droplevels(subset(mean_endogamy, ME > 0.20)) # select high endogamy families (F>0.45)


# Filter and drop families with F_roh <0.45
df2 = df1_ordered %>%
  filter(FAMILY %in% levels(sel_fam$FAMILY)) %>%
  mutate(FAMILY = droplevels(FAMILY))
dim(df2)
levels(df2$FAMILY)

df2 = df2[order(df2$F_roh, decreasing = TRUE), ]
df2 = df2[order(factor(df2$FAMILY,levels = mean_endogamy_order$FAMILY)), ]


# Select top 3 individuals with greatest F_roh values per Family
top_ind <- df2 %>%
  group_by(FAMILY) %>% # Group by Family
  top_n(3, F_roh) %>%                            # Select the top 3 by F_roh
  ungroup()                                      # Ungroup to return to original structure
dim(top_ind)

write.table(top_ind, "outputs/30_highF_inds.txt")


mean_F_sel = top_ind %>%
  group_by(FAMILY) %>%
  summarise(MeanEndogamy = mean(F_roh))

mean_F_sel = mean_F_sel[order(mean_F_sel$MeanEndogamy, decreasing = T), ]

mean(mean_F_sel$MeanEndogamy)
#0.5192871


# Using BLUP --------------------------------------------------------------

# Regression g in F - 3years-----------------------------------------------------

# Phenotypic data 3 years
# pheno = readxl::read_excel("data_raw/Dados_phen_3anos.xlsx")
# str(pheno)
# 
# pheno_roh = pheno[pheno$fator2%in%F_roh$fator2, ]
# dim(pheno_roh)
# 
# 
# 
# # Phenotypic data 3 years
# eBLUP = read.table("outputs/eBLUP_all.txt")
# str(eBLUP)
# 
# eBLUP = eBLUP[rownames(eBLUP)%in%F_roh$fator2, ]
# dim(eBLUP)
# head(eBLUP)
# eBLUP$fator2 = rownames(eBLUP)
# 
# df4 = merge(eBLUP, F_roh, by="fator2")
# pheno_roh = pheno_roh[pheno_roh$fator2%in%df4$fator2, ]
# df4 = merge(df4, pheno_roh, by="fator2")
# dim(df4)
# 
# df4 = df4 %>% rename(Family = fator1, IND = fator2)
# 
# 
# 
# # Calculate R-squared for each family
# r2_data = df4 %>%
#   group_by(Family) %>%
#   summarise(
#     R2 = summary(lm(a ~ F_roh))$r.squared
#   )
# 
# 
# mean_endogamy = df4 %>%
#   group_by(Family) %>%
#   summarise(MeanEndogamy = mean(F_roh))
# 
# 
# # Reorder the Family factor based on R^2 (descending order)
# df4$Family = factor(df4$Family, levels = mean_endogamy$Family[order(mean_endogamy$MeanEndogamy, decreasing = T)])
# 
# 
# # Create the faceted plot
# ggplot(df4, aes(x = F_roh, y = a)) +
#   geom_point(color = "#069C56", alpha = 0.7, size=3) +            # Scatter plot
#   geom_smooth(method = "lm", color = "#FF681E", se = TRUE, 
#               fill = "darkgray", alpha = 0.2) +         #Regression line    
#   #geom_vline(data = mean_endogamy, aes(xintercept = MeanEndogamy), 
#   #color = "darkgreen", linetype = "dashed", size = 0.8) +  # Mean Endogamy line
#   facet_wrap(~ Family) +                               # Facets sorted by R2
#   geom_text(data = mean_endogamy, 
#             aes(x = 0.4, y = 140, label = paste0("F = ", round(MeanEndogamy, 2))), 
#             inherit.aes = FALSE, size = 4, color = "black") +  # Add Mean Endogamy
#   geom_text(data = r2_data, 
#             aes(x = 0.4, y = 110, label = paste0("R² = ", round(R2, 2))), 
#             inherit.aes = FALSE, size = 4, color = "black") +  # Add R2
#   
#   labs(
#     title = "Regression of F_ROH on genotypic value (g) by Family",
#     x = "F_ROH",
#     y = "DAP"
#   ) +
#   theme_minimal()



# BLUP on FHBD (a, d e g) ----------------------------------------------------


eBLUP = read.table("results/BLUPS-AD-GBLUP.txt")

# Perform approximate row name matching
library(stringdist)
matches = stringdistmatrix(df1_ordered$IND, rownames(eBLUP), method = "jw")
closest_match = apply(matches, 1, which.min)


# Perform approximate matching between unique row names in df1 and df2
unique_names_df1 = unique(df1_ordered$IND)
matches = stringdistmatrix(unique_names_df1, rownames(eBLUP), method = "jw")
closest_match = apply(matches, 1, which.min)

# Create a mapping of matched row names
mapping = data.frame(
  Row1 = unique_names_df1,
  Row2 = rownames(eBLUP)[closest_match]
)

# Merge the matched values from df2 into df1
df1_ordered$Row = df1_ordered$IND  # Preserve original row names
df1_matched = merge(
  df1_ordered, mapping, by.x = "Row", by.y = "Row1", all.x = TRUE
)  # Add the corresponding row names from df2
df2_matched = eBLUP[mapping$Row2, ]  # Extract matching rows from df2
df2_matched$Row = rownames(df2_matched)

# Merge matched data from df2 into df1
final_df = merge(
  df1_matched, df2_matched, by = "Row", all.x = TRUE
)



# Reorder df1 based on the factor levels of df2$COLUMN
final_df = final_df %>%
  mutate(Pai = factor(Pai, levels = levels(mean_endogamy_order$Pai))) %>%
  arrange(Pai)


# F on a

# Calculate R-squared for each family
r2_data = final_df %>%
  group_by(Pai) %>%
  summarise(
    R2 = summary(lm(a ~ F_roh))$r.squared
  )

ggplot(final_df, aes(x = F_roh, y = a)) +
  geom_point(color = "#069C56", alpha = 0.7, size=3) +            # Scatter plot
  geom_smooth(method = "lm", color = "#FF681E", se = TRUE, 
              fill = "darkgray", alpha = 0.2) +         #Regression line    
  xlim(0, 0.6)+
  ylim(-5, 3)+
  facet_wrap(~ Pai) +                               # Facets sorted by R2
  geom_text(data = mean_endogamy, 
            aes(x = 0.5, y = -4, label = paste0("F = ", round(ME, 2))), 
            inherit.aes = FALSE, size = 3, color = "black") +  # Add Mean Endogamy
  geom_text(data = r2_data, 
            aes(x = 0.2, y = -4, label = paste0("R² = ", round(R2, 2))), 
            inherit.aes = FALSE, size = 3, color = "black") +  # Add R2
  labs(
    title = "Regression of F_ROH on a by Pai",
    x = "F_ROH",
    y = "a"
  ) +
  theme_minimal()



# F on d

# Calculate R-squared for each family
r2_data = final_df %>%
  group_by(Pai) %>%
  summarise(
    R2 = summary(lm(d ~ F_roh))$r.squared
  )

ggplot(final_df, aes(x = F_roh, y = d)) +
  geom_point(color = "#069C56", alpha = 0.7, size=3) +            # Scatter plot
  geom_smooth(method = "lm", color = "#FF681E", se = TRUE, 
              fill = "darkgray", alpha = 0.2) +         #Regression line    
  xlim(0, 0.6)+
  ylim(-5, 3)+
  facet_wrap(~ Pai) +                               # Facets sorted by R2
  geom_text(data = mean_endogamy, 
            aes(x = 0.5, y = -4, label = paste0("F = ", round(ME, 2))), 
            inherit.aes = FALSE, size = 3, color = "black") +  # Add Mean Endogamy
  geom_text(data = r2_data, 
            aes(x = 0.2, y = -4, label = paste0("R² = ", round(R2, 2))), 
            inherit.aes = FALSE, size = 3, color = "black") +  # Add R2
  labs(
    title = "Regression of F_ROH on d by Pai",
    x = "F_ROH",
    y = "d"
  ) +
  theme_minimal()

# F on g


# Calculate R-squared for each family
r2_data = final_df %>%
  group_by(Pai) %>%
  summarise(
    R2 = summary(lm(g ~ F_roh))$r.squared
  )


ggplot(final_df, aes(x = F_roh, y = g)) +
  geom_point(color = "#069C56", alpha = 0.7, size=3) +            # Scatter plot
  geom_smooth(method = "lm", color = "#FF681E", se = TRUE, 
              fill = "darkgray", alpha = 0.2) +         #Regression line    
  xlim(0, 0.6)+
  ylim(-7, 3)+
  facet_wrap(~ Pai) +                               # Facets sorted by R2
  geom_text(data = mean_endogamy, 
            aes(x = 0.5, y = -6.5, label = paste0("F = ", round(ME, 2))), 
            inherit.aes = FALSE, size = 3, color = "black") +  # Add Mean Endogamy
  geom_text(data = r2_data, 
            aes(x = 0.2, y = -6.5, label = paste0("R² = ", round(R2, 2))), 
            inherit.aes = FALSE, size = 3, color = "black") +  # Add R2
  labs(
    title = "Regression of F_ROH on g by Pai",
    x = "F_ROH",
    y = "a"
  ) +
  theme_minimal()

