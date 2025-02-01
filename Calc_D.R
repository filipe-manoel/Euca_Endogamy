#Cálculo de D e F


calc_D <- function(combinacoes, ciclos) {
  
  # Definições iniciais
  corr = 3 # Correção relativa ao ROH (Ele não considera os últimos 3 indivíduos)
  
  # Criar uma lista para armazenar os resultados
  resultados_lista <- list()
  
  # Loop para calcular os resultados para cada combinação e para os ciclos
  for (ciclo in 1:ciclos) {
    for (i in 1:nrow(combinacoes)) {
      nAuto <- combinacoes$nAuto[i]
      nCruz <- combinacoes$nCruz[i]
      
      # Seleção dos indivíduos autorreprodutores
      id_auto_sel <- sample(auto_pop@id, nAuto, replace = FALSE)
      auto_pop_sel <- auto_pop[auto_pop@id %in% id_auto_sel]
      
      # Seleção dos indivíduos de correção
      corr_sel <- sample(hib_pop@id, corr, replace = FALSE)
      corr_pop <- hib_pop[hib_pop@id %in% corr_sel]
      
      # Seleção aleatória dos indivíduos para cruzamento
      id_sel <- sample(hib_pop@id, nCruz, replace = FALSE) 
      hib_pop_sel <- hib_pop[hib_pop@id %in% id_sel]
      
      # Unir populações
      populacao <- mergePops(list(auto_pop_sel, hib_pop_sel, corr_pop))
      populacao <- setPheno(populacao, H2)
      
      # Calcular variância genética
      varG(populacao)
      
      # Calcular o coeficiente de endogamia
      ROH <- ROH_POP(populacao)
      coef_endogamia <- ROH$Coef_Endogamia_Total
      F <- mean(coef_endogamia$Coef_Endogamia_Total)
      
      # Cálculos de D e DG
      S1 <- (1 / F) + 2
      D <- (100 / meanP(populacao)) * ((meanP(hib_pop_sel) - meanP(auto_pop_sel)) / ((1 / (2 - S1)) - (S1 / (2 - S1))))
      DG <- (100 / meanG(populacao)) * ((meanG(hib_pop_sel) - meanG(auto_pop_sel)) / ((1 / (2 - S1)) - (S1 / (2 - S1))))
      
      # Cálculos de D2 e DG2
      D2 <- (-meanP(auto_pop_sel) + meanP(hib_pop_sel)) / meanP(hib_pop_sel)*100
      DG2 <- (-meanG(auto_pop_sel) + meanG(hib_pop_sel)) / meanG(hib_pop_sel)*100
      
      # Adicionar os resultados à lista com o ciclo
      resultados_lista[[length(resultados_lista) + 1]] <- c(D, DG, D2, DG2, F, nAuto, nCruz, ciclo)
    }
  }
  
  # Converter a lista em um dataframe
  resultados <- do.call(rbind, resultados_lista)
  
  # Definir os nomes das colunas
  colnames(resultados) <- c("D", "DG", "D2", "DG2", "F", "nAuto", "nCruz", "Ciclo")
  
  # Retornar o dataframe com os resultados
  return(resultados)
}

#------------------------------------------------------------------------------#
# Gerar os valores de nAuto e nCruz
nAuto <- seq(29, 15, by = -1) # De 30 a 15
nCruz <- seq(1, 15, by = 1)  # De 0 a 15

# Criar o dataframe com as combinações
combinacoes <- data.frame(nAuto = nAuto, nCruz = nCruz)



# Chamar a função para obter os resultados
resultados <- calc_D(combinacoes, ciclos)



library(dplyr)
# Adicionar a coluna ao dataframe resultados
resultados <- resultados %>%
  as.data.frame() %>% # Caso resultados não seja um data.frame
  mutate(Relacao = nAuto / (nAuto + nCruz))

# Exibir as primeiras linhas do dataframe atualizado
head(resultados)



estatisticas <- resultados %>%
  as.data.frame() %>% 
  group_by(nAuto, nCruz) %>%
  summarise(
    D_mean = mean(D), D_min = min(D), D_max = max(D), D_var = var(D),
    DG_mean = mean(DG), DG_min = min(DG), DG_max = max(DG), DG_var = var(DG),
    D2_mean = mean(D2), D2_min = min(D2), D2_max = max(D2), D2_var = var(D2),
    DG2_mean = mean(DG2), DG2_min = min(DG2), DG2_max = max(DG2), DG2_var = var(DG2),
    F_mean = mean(F), F_min = min(F), F_max = max(F), F_var = var(F),
    Relacao = unique(nAuto / (nAuto + nCruz)) # Relacao constante por grupo
  ) %>%
  ungroup()

#------------------------------------------------------------------------------#

#MSE

MSE <- data.frame(nAuto = resultados$nAuto, nCruz = resultados$nCruz, 
                  MSE = ((resultados$D - estatisticas$DG_mean[which.min(estatisticas$DG_var)])^2) / nrow(resultados))

MSE <- MSE %>% 
  group_by(nAuto, nCruz) %>% 
  summarise(MSE_value = sum(MSE))


estatisticas$MSE = MSE$MSE_value

MSE2 <- data.frame(nAuto = resultados$nAuto, nCruz = resultados$nCruz, 
                  MSE2 = ((resultados$D2 - estatisticas$DG2_mean[which.min(estatisticas$DG2_var)])^2) / nrow(resultados))

MSE2 <- MSE2 %>% 
  group_by(nAuto, nCruz) %>% 
  summarise(MSE_value2 = sum(MSE2))


estatisticas$MSE2 = MSE2$MSE_value2

#RMSE

estatisticas$RMSE = sqrt(estatisticas$MSE)
estatisticas$RMSE2 = sqrt(estatisticas$MSE2)


