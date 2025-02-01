
rm(list=ls())
gc()


# Library -----------------------------------------------------------------


# List of packages to load
packages <- c("tidyverse", "asreml", "lme4", "car", "snpReady", "AGHmatrix",
              "ASRgenomics", "redxl", "data.table", "gridExtra", "PupillometryR",
              "readxl")


#install.packages("StageWise")

# Use pacman to load and install packages if needed
pacman::p_load(char = packages)

# Load functions 
#source("ic_reml_asr.R")

# Phenotypic data ---------------------------------------------------------

pheno <- read_excel("data_raw/Dados_phen_3anos.xlsx")
head(pheno)
str(pheno)
dim(pheno)

pheno[1:19] = lapply(pheno[1:19], as.factor) #transform in factors
str(pheno)

# Gmatrix -----------------------------------------------------------------

# Obtain G.inv


# SNPs data  ----------------------------------------------------------

SNPs <- read.table(file = "data_raw/snps_todos_ind.csv", 
                   sep = ";",
                   quote = "\"",
                   header = T)
SNPs[1:5, 1:10]
dim(SNPs)

geno.mat = SNPs[ , -c(1:4)]
geno.mat[1:5, 1:10]

rownames(geno.mat)=SNPs$ID
colnames(geno.mat)=colnames(SNPs [ , -c(1:4)])
geno.mat = as.matrix(geno.mat)
geno.mat[geno.mat == 0] = NA


M <- snpReady::raw.data(
  data = geno.mat,
  frame = "wide",
  base = TRUE,
  sweep.sample = 0.95,
  call.rate = 0.90,
  maf = 0.01,
  imput = FALSE,
  imput.type = "wright",
  outfile = "012"# the allel with minor frequency receives 2
)

M1 <- M$M.clean |>
  round()


M1[1:5, 1:10]

G <- ASRgenomics::G.matrix(M = M1, method = "VanRaden")$G
attr(G, "rowNames") <- as.character(rownames(M1))
attr(G, "colNames") <- as.character(rownames(M1))

G[1:5, 1:5]
dim(G)

saveRDS(G, "data/matrices/G.rds")

# Get A for blending

# Get pedigree file from phenotype file.

ped_c = read_xlsx("data_raw/pedigree_corrigido.xlsx", sheet = "pedC", col_names = TRUE)
ped_c = ped_c[!duplicated(ped_c$IND), ]
dim(ped_c)

head(ped_c)
ped_c = ped_c %>% select(IND,Father, Mother)
head(ped_c)
ped_c= data.frame(ped_c)
ped_c[is.na(ped_c)] = 0
prep.A = nadiv::prepPed(ped_c)
A = as.matrix(nadiv::makeA(prep.A))
dim(A)
A[1:5, 1:5]

saveRDS(A, "data/matrices/A.rds")

A = A[order(rownames(A)), order(colnames(A))]
G = G[order(rownames(G)), order(colnames(G))]

#Subsetting the matrix to only individuals present in the matrix Ga (i.e., excluding the parents).
A <- A[rownames(A) %in% rownames(G),  colnames(A) %in% colnames(G)]  #droplevels in rows of missing genotypes
G <- G[rownames(G) %in% rownames(A),  colnames(G) %in% colnames(A)]  #droplevels in rows of missing genotypes

#Making sure that the matrices are in the same order before blending:
dim(ped_c)
dim(A)
dim(G)

all(rownames(A) == rownames(G)) #order and name of individuals

Aclean <- match.G2A(A = as.matrix(A), G = as.matrix(G), clean = TRUE, ord = TRUE, mism = TRUE)$Aclean
Aclean[1:5, 1:5]
dim(Aclean)

G.blend <- G.tuneup(G = G, A= Aclean, blend = TRUE, pblend = 0.02)$Gb
G.inv <- G.inverse(G = G.blend, sparseform = TRUE)$Ginv  # sparseform = TRUE

attr(G.inv,"INVERSE")<- TRUE
attr(G.inv, "rowNames")
head(G.inv)
saveRDS(G.inv, "data/matrices/G.inv.rds")


# H matrix ----------------------------------------------------------------

library(AGHmatrix)
#Computing Ac matrix (Munoz)
H_Munoz <- Hmatrix(A=Aclean, G=G, markers = M1,
                   ploidy=2, method="Munoz",
                   roundVar=2,
                   maf=0.05)
# ‘Munoz‘ shrinks the G matrix towards the A matrix scaling the molecular relatadness
#by each relationship classes; 

H_Munoz.inv = solve(H_Munoz) #Get inverse of matrix Ga_blended by solving.

## Getting ready for ASReml-r V4 package
#Transform to sparse matrix (i.e., sorting in a different way, only records were there are no zeros).
H_Munoz.inv <-as(H_Munoz.inv ,"sparseMatrix")

#Get lower diagonal of the matrix and get ready for use in ASReml (package MCMCglmm).
H_Munoz.inv  <-MCMCglmm::sm2asreml(H_Munoz.inv)

#Set "INVERSE" as an attribute. This is needed for ASReml-R V4.
attr(H_Munoz.inv, "INVERSE") <- TRUE
attr(H_Munoz.inv, "rowNames")
dim(H_Munoz.inv)

saveRDS(H_Munoz.inv, "data/matrices/H_Munoz.inv.rds")


# H matrix ----------------------------------------------------------------
#Computing H matrix (Martini)/ Legarra)2009)

Hmat_Legarra <- Hmatrix(A=Aclean, G=G, method="Martini",
                        ploidy=2,
                        maf=0.05,
                        tau = 1,
                        omega = 1)


# ‘Martini‘ is a modified version from Legarra et al. (2009) where
# combines A and G matrix using scaling factors. When method is equal ‘Martini‘ and ‘tau=1‘ and
# ‘omega=1‘ you have the same H matrix as in Legarra et al. (2009).

H_Legarra.inv = solve(Hmat_Legarra) #Get inverse of matrix Ga_blended by solving.

## Getting ready for ASReml-r V4 package
#Transform to sparse matrix (i.e., sorting in a different way, only records were there are no zeros).

H_Legarra.inv <-as(H_Legarra.inv,"sparseMatrix")

#Get lower diagonal of the matrix and get ready for use in ASReml (package MCMCglmm).

H_Legarra.inv <-MCMCglmm::sm2asreml(H_Legarra.inv)

#Set "INVERSE" as an attribute. This is needed for ASReml-R V4.
attr(H_Legarra.inv, "INVERSE") <- TRUE
attr(H_Legarra.inv, "rowNames")

saveRDS(H_Legarra.inv, "data/matrices/H_Legarra.inv.rds")


# Dominance Matrix --------------------------------------------------------

## Calculating the dominant genomic relationship matrix (Gd) and its inverse

#Obtain the matrix Gd (package AGHmatrix).
D = Gmatrix(SNPmatrix = as.matrix(M1), method = "Vitezica")
D[1:5, 1:5]
# It is not possible to obtain the pedigree-based dominant relationship matrix (D)
# for blending, since we are working with a half-sif families

attr(D, "rowNames") <- as.character(rownames(M1))
attr(Gd, "colNames") <- as.character(colnames(M))

D.inv <- G.inverse(G = D, sparseform = TRUE)$Ginv
head(D.inv)


attr(D.inv,"INVERSE")<-TRUE
attr(D.inv, "rowNames")
saveRDS(D.inv, "data/matrices/D.inv.rds")

