# Função para realizar cálculo do coeficiente de endogamia via ROH para uma população AlphaSimR
ROH_POP <- function(pop) {
  # Instalar e carregar o pacote RZooRoH
  if (!requireNamespace("RZooRoH", quietly = TRUE)) {
    install.packages("RZooRoH")
  }
  library(RZooRoH)
  
  # Extração dos genótipos e mapa genético
  snps <- pullSnpGeno(pop)
  genMap <- getGenMap()
  
  # Selecionar os IDs dos marcadores no dataframe `snps`
  marcadores_snps <- colnames(snps)
  
  # Filtrar o dataframe `genMap` para incluir apenas os marcadores existentes em `snps`
  genMap <- genMap[genMap$id %in% marcadores_snps, ]
  
  # Validar e ajustar a ordem dos marcadores em genMap
  genMap <- genMap[match(marcadores_snps, genMap$id), ]
  
  # Garantir consistência das dimensões
  if (nrow(genMap) != ncol(snps)) {
    stop("Inconsistência entre genMap e snps após a filtragem!")
  }
  
  # Verificar se os genótipos são numéricos e converter, se necessário
  if (!is.numeric(snps)) {
    snps <- apply(snps, 2, as.numeric)
  }
  
  # Preparar os dados para o RZooRoH
  dados_formatados <- data.frame(
    chr = genMap$chr,
    pos = genMap$pos,
    t(snps)  # Transpor os genótipos para o formato de colunas por indivíduos
  )
  
  # Salvar o arquivo formatado em um arquivo temporário
  temp_file <- tempfile(fileext = ".txt")
  write.table(dados_formatados, temp_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Ler os dados no RZooRoH
  zooin_obj <- zoodata(
    genofile = temp_file,
    zformat = "gt",
    chrcol = 1,
    poscol = 2
  )
  
  # Definir o modelo HMM com 5 classes HBD
  zmodel <- zoomodel()
  
  # Executar o modelo para identificar segmentos HBD e calcular autozigosidade
  resultados <- zoorun(zmodel, zooin_obj, fb = TRUE)
  
  # Coeficiente de Endogamia (Autozigosidade Realizada)
  autozygosity <- resultados@realized
  
  # Coeficiente de Endogamia Total (1 - Proporção não-HBD)
  coef_endogamia_total <- 1 - autozygosity[, ncol(autozygosity)]
  
  # Coeficiente de Endogamia Cumulativo (até T = 4)
  coef_endogamia_cumulativo <- t(apply(autozygosity[, 1:2], 1, cumsum))
  
  # Geração dos gráficos
  par(mfrow = c(1, 2)) # Configuração para dois gráficos lado a lado
  
  # Gráfico do coeficiente de endogamia total
  hist(
    coef_endogamia_total, 
    breaks = 20, 
    main = "Coeficiente de Endogamia Total",
    xlab = "Inbreeding Coefficient (Total)", 
    col = "coral2"
  )
  
  # Gráfico do coeficiente de endogamia cumulativo
  hist(
    coef_endogamia_cumulativo[, 2], 
    breaks = 20, 
    main = "Coeficiente de Endogamia (T = 4)",
    xlab = "Inbreeding Coefficient (T = 4)", 
    col = "tomato"
  )
  
  # Retornar os coeficientes como uma lista
  return(list(
    Coef_Endogamia_Total = data.frame(
      Individuo = resultados@sampleids,
      Coef_Endogamia_Total = coef_endogamia_total
    ),
    Coef_Endogamia_Cumulativo = data.frame(
      Individuo = resultados@sampleids,
      Coef_Endogamia_T4 = coef_endogamia_cumulativo[, 2]
    )
  ))
}

# Exemplo de uso:
# resultados_endogamia <- ROH_POP(familias_pop[["Familia_737"]])
# print(resultados_endogamia$Coef_Endogamia_Total)
# print(resultados_endogamia$Coef_Endogamia_Cumulativo)
