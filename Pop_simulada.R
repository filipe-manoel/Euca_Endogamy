#Obtenção da população simulada



# Variáveis para armazenar o maior valor de area_sobreposicao
maior_area <- 0  # Começa com um valor muito baixo
contador <- 0  # Contador de execuções bem-sucedidas

# Loop para rodar o código até atingir o número de repetições bem-sucedidas
while (contador < rep) {
  
  tryCatch({
    
    nInd = 500
    # Parâmetros Iniciais
    foundergenomes <- runMacs(nInd = nInd, 
                              nChr = 11, 
                              segSites = 1000,
                              species = "GENERIC",
                              split = 2,
                              ploidy = 2)
    SP <- SimParam$new(foundergenomes)
    
    SP$altAddTraitAD(nQtlPerChr = 600, 
                     mean = 1, 
                     varA = 3, 
                     varD = 1, 
                     inbrDepr = 10, 
                     limMeanDD = c(0.1, 0.3),
                     limVarDD = c(0.5, 1),
                     silent = FALSE)
    
    SP$addSnpChip(100)
    
    # População Base
    nPro <- 10
    H2 <- 0.40
    
    # Base population
    basepop1 <- newPop(foundergenomes, SP)
    basepop1 <- setPheno(basepop1, H2 = H2, reps = 1)
    
    #------------------------------------------------------------------------------#
    #Cruzamentos aleatórios com taxa de autofecundação (Simulando natureza)
    
    # Número de iterações
    nRand <- 100
    
    # Parâmetros gerais
    nSeeds <- 1         # Número de sementes por planta
    probSelf <- 0.8     # Probabilidade de autofecundação (20%)
    trait <- 1          # Traço a ser usado para seleção
    use <- "pheno"      # Critério de seleção (fenótipos)
    H2 <- 0.40           # Herdabilidade para setPheno
    
    rand_sel <- function(pop, nRand, nInd, nSeeds, probSelf, H2, SP) {
      for (i in 1:nRand) {
        # Realiza a seleção
        pop <- selectOP(pop, nInd = nInd, nSeeds = nSeeds, probSelf = probSelf, simParam = SP)
        
        # Atualiza os fenótipos
        pop <- setPheno(pop, H2 = H2, reps = 1)
        
        print(meanG(pop))
        
        # Opcional: Mostra progresso
        cat("Iteração", i, "concluída\n")
      }
      
      # Retorna a população final
      return(pop)
    }
    
    nCycles = 1
    
    
    basepop1 = rand_sel(basepop1,nRand, nInd, nSeeds, probSelf, H2, SP)
    basepop1 = setPheno(basepop1, H2)
    basepop1 = selectInd(basepop1, 10)
    
    #------------------------------------------------------------------------------#
    #Seleção recorrente
    
    nInd = 10
    # Função para executar um ciclo de seleção recorrente
    rec_selection <- function(pop, nCycles, H2, nPro, nInd) {
      populations <- list() # Lista para armazenar as populações de cada ciclo
      
      for (cycle in 1:nCycles) {
        cat("Iniciando ciclo:", cycle, "\n")
        
        # Armazena a população atual antes da seleção
        populations[[paste0("Cycle_", cycle)]] <- pop
        
        # Plano de cruzamento
        crossplan <- t(as.data.frame(combn(pop@id, 2)))
        pop <- makeCross(pop, crossplan, nProgeny = nPro)
        
        # Fenotipagem
        pop <- setPheno(pop, H2)
        print(meanG(pop))
        
        # Seleção
        pop <- selectInd(pop, nInd, use = "pheno")
        #pop <- self(pop, 1)
        #pop <- self(pop, 1)
      }
      
      # Armazena a população final após o último ciclo
      populations[[paste0("Cycle_", nCycles + 1)]] <- pop
      
      return(populations)
    }
    
    # Executar o programa de seleção recorrente
    all_pops <- rec_selection(basepop1, nCycles, H2, nPro, nInd)
    
    # Resultado final
    cat("Seleção recorrente concluída!\n")
    
    # Acessar populações
    for (cycle in names(all_pops)) {
      cat(cycle, ": População com", nInd(all_pops[[cycle]]), "indivíduos\n")
    }
    
    #------------------------------------------------------------------------------#
    #População 2
    # Parâmetros Iniciais para a População 2
    
    nInd2 <- nInd
    
    
    foundergenomes2 <- runMacs(nInd = nInd2, 
                               nChr = 11, 
                               segSites = 1000,
                               species = "GENERIC",
                               split = 2,
                               ploidy = 2)
    SP2 <- SimParam$new(foundergenomes2)
    
    SP2$altAddTraitAD(nQtlPerChr = 600, 
                      mean = 1, 
                      varA = 3, 
                      varD = 1,  
                      inbrDepr = 10, 
                      limMeanDD = c(0.1, 0.3),
                      limVarDD = c(0.5, 1),
                      silent = FALSE)
    
    SP2$addSnpChip(100)
    
    # População Base 2
    nPro2 <- 10
    H2_2 <- 0.40
    
    basepop2 <- newPop(foundergenomes2, SP2)
    basepop2 <- setPheno(basepop2, H2 = H2_2, reps = 1)
    
    basepop2 = rand_sel(basepop2,nRand, nInd, nSeeds, probSelf, H2, SP)
    basepop2 = setPheno(basepop2, H2)
    basepop2 = selectInd(basepop2, 10)
    
    nCycles2 <- 1
    
    # Executar o programa de seleção recorrente para a População 2
    all_pops2 <- rec_selection(basepop2, nCycles2, H2_2, nPro2, nInd2)
    
    # Resultado final para a População 2
    cat("Seleção recorrente concluída para a População 2!\n")
    
    # Acessar populações da População 2
    for (cycle in names(all_pops2)) {
      cat(cycle, ": População com", nInd(all_pops2[[cycle]]), "indivíduos\n")
    }
    
    #------------------------------------------------------------------------------#
    nPro_hybrid = 30
    H2 = 0.4
    
    crossplan <- expand.grid(all_pops[[(nCycles+1)]]@id, all_pops2[[(nCycles2+1)]]@id) # Todos cruzamentos entre pop1 e pop2
    crossplan <- as.matrix(crossplan)
    pop = makeCross2(all_pops[[(nCycles+1)]], all_pops2[[(nCycles2+1)]], crossplan, nProgeny = nPro_hybrid)
    
    pop = setPheno(pop,H2)
    pop = selectInd(pop, 20, use = "pheno")
    
    #------------------------------------------------------------------------------#
    #Pop base para seleção recorrente
    
    ciclos = 1:10
    for (i in ciclos){
      crossplan <- t(as.data.frame(combn(pop@id, 2)))
      pop <- makeCross(pop, crossplan, nProgeny = 10)
      pop <- setPheno(pop, H2)
      pop = selectInd(pop, 5)
      print(meanG(pop))
      
    }
    crossplan <- t(as.data.frame(combn(pop@id, 2)))
    pop <- makeCross(pop, crossplan, nProgeny = 10)
    pop <- setPheno(pop, H2)
    pop = selectWithinFam(pop, 1)
    
    pop <- self(pop,nProgeny = 30)
    pop <- setPheno(pop, H2)
    pop = selectWithinFam(pop, 1)
    pop = selectFam(pop, 6)
    
    #ROH_pop = ROH_POP(pop)
    #varG(pop)
    #mean(ROH_pop$Coef_Endogamia_Total$Coef_Endogamia_Total)
    
    #Autofecundações
    parental = selectInd(pop, 1)
    auto_pop = self(parental, nProgeny = 100)
    auto_pop = setPheno(auto_pop,H2)
    
    #ROH_auto = ROH_POP(auto_pop)
    #mean(ROH_auto$Coef_Endogamia_Total$Coef_Endogamia_Total)
    
    
    #Híbridos
    
    crossplan <- expand.grid(parental@id, pop@id) 
    crossplan <- as.matrix(crossplan)
    crossplan = crossplan[-1,]
    hib_pop = makeCross2(parental,pop, crossplan, nProgeny = 10)
    hib_pop = setPheno(hib_pop, H2)
    
    
    #------------------------------------------------------------------------------#
    
    finalpop = populacao <- mergePops(list(auto_pop, hib_pop))
    
    print(finalpop)
    
    
    # Verifica se finalpop@pheno é um vetor ou matriz
    if (is.vector(finalpop@pheno)) {
      pheno_data <- finalpop@pheno
    } else {
      pheno_data <- finalpop@pheno[, 1]  # Usando a primeira coluna, se for matriz
    }
    
    # Calculando a densidade da distribuição de finalpop@pheno
    densidade_pheno <- density(pheno_data)
    
    # Definindo os parâmetros da distribuição normal
    media <- 12.5
    variancia <- 20.45
    desvio_padrao <- sqrt(variancia)
    
    # Criando uma sequência de valores no eixo x para a distribuição normal
    x <- seq(2.20, 20.80, length.out = 10000)
    
    # Calculando a função densidade da distribuição normal
    y <- dnorm(x, mean = media, sd = desvio_padrao)
    
    
    # Plotando a densidade de finalpop@pheno com limite no eixo Y
    plot(densidade_pheno, 
         main = "Distribution for the Simulated and Real Data",
         xlab = "Fenotype", 
         ylab = "Density",
         col = "steelblue", 
         lwd = 2,
         ylim = c(0, 0.15),
         xlim = c(0, 25)
    )  # Define o limite do eixo Y
    
    # Adicionando a curva da distribuição normal teórica
    lines(x, y, col = "darkgreen", lwd = 2, lty = 2)
    
    # Adicionando uma linha vertical para marcar a média de finalpop@pheno
    abline(v = mean(pheno_data), col = "red", lwd = 2, lty = 2)
    
    # Adicionando uma linha vertical para marcar a média da distribuição normal teórica
    abline(v = media, col = "darkgreen", lwd = 2, lty = 2)
    
    # Adicionando legenda
    legend("topright", legend = c("Simulation", "Real Data"),
           col = c("steelblue", "darkgreen"), lwd = 2, lty = c(1, 2))
    
    # Salvando o gráfico como um objeto
    grafico_densidade <- recordPlot()
    
    # Para visualizar o gráfico salvo posteriormente
    replayPlot(grafico_densidade)
    
    
    #------------------------------------------------------------------------------#
    
    # Verifica se finalpop@pheno é um vetor ou matriz
    if (is.vector(finalpop@pheno)) {
      pheno_data <- finalpop@pheno
    } else {
      pheno_data <- finalpop@pheno[, 1]  # Usando a primeira coluna, se for matriz
    }
    
    # Calculando a densidade da distribuição observada (finalpop@pheno)
    densidade_pheno <- density(pheno_data)
    
    # Definindo os parâmetros da distribuição normal teórica
    media <- 12.5
    variancia <- 20.46
    desvio_padrao <- sqrt(variancia)
    
    # Criando uma sequência de valores para o eixo x
    x <- densidade_pheno$x
    
    # Calculando a densidade da distribuição normal teórica nos mesmos pontos de x
    y_theoretical <- dnorm(x, mean = media, sd = desvio_padrao)
    
    # Calculando a área de sobreposição
    area_sobreposicao <- sum(pmin(densidade_pheno$y, y_theoretical)) * (x[2] - x[1])
    
    # Exibindo a área de sobreposição
    print(paste("Ciclo Atual:", contador))
    
    print(paste("Atual area_sobreposicao:", area_sobreposicao))
    
    
    if (area_sobreposicao > maior_area) {
      maior_area <- area_sobreposicao
      
      # Salva os objetos atuais em um arquivo .RData
      save.image(file = paste0("melhor_repeticao.RData"))
      grafico = grafico_densidade 
      mediaP = meanP(finalpop)
      VA = varA(finalpop)
      VD = varD(finalpop)
      VP = varP(finalpop)
    }
    
    print(paste("Maior area_sobreposicao:", maior_area))
    
    # Incrementa o contador apenas se a execução for bem-sucedida
    contador <- contador + 1
    
  }, error = function(e) {
    # Mensagem de erro (opcional)
    message("Erro capturado: ", conditionMessage(e))
    # Continua o loop sem incrementar o contador
  })
}


