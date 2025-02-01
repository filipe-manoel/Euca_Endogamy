#Endogamia

rm(list = ls())
library(AlphaSimR)
library(ggplot2)
library(ggfortify) 
library(beepr)
library(openxlsx)
library(svglite)

setwd("Simulation")
getwd()
#------------------------------------------------------------------------------#
rep <- 5000
source("Pop_simulada.R")
source("ROH_AlphaSimR_POP.R")


# Exibe o maior valor de "area_sobreposicao" encontrado
print(paste("Maior area_sobreposicao:", maior_area))
print(grafico)

mediaP
VA 
VD 
VP

#------------------------------------------------------------------------------#
# Carregando o arquivo melhor_repeticao.RData
load("C:/Users/LENOVO/OneDrive/POS-DOC/Endogamia/Paper_ROH/codes/melhor_repeticao.RData") #89%

# Definir o número de ciclos
ciclos <- 100
source("C:/Users/LENOVO/OneDrive/POS-DOC/Endogamia/Paper_ROH/Simulation/Calc_D.R")
head(estatisticas)

#------------------------------------------------------------------------------#
library(tidyverse)
source("C:/Users/LENOVO/OneDrive/POS-DOC/Endogamia/Paper_ROH/Simulation/Gráficos.R")

#------------------------------------------------------------------------------#
library(openxlsx)

# Salvar o dataframe como um arquivo Excel
write.xlsx(estatisticas, file = "Simulation/results/resultados.xlsx")

save.image("Simulation/results/resultadosR.RData")


# Salvando os gráficos de boxplot
ggsave("Simulation/results/boxplot.svg", plot = box, width = 8, height = 6, dpi = 300)
ggsave("Simulation/results/boxplot2.svg", plot = boxD2, width = 8, height = 6, dpi = 300)


# Salvando os gráficos de RMSE
ggsave("Simulation/results/rmse_plot.svg", plot = rmse, width = 8, height = 6, dpi = 300)
ggsave("Simulation/results/rmse_plot2.svg", plot = rmse2, width = 8, height = 6, dpi = 300)

#------------------------------------------------------------------------------#
getwd()

beep(3)
