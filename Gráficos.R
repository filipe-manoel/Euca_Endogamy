#Gráficos

# Encontrando o valor de DG_mean com a menor variância

dg_mean_value <- estatisticas$DG_mean[which.min(estatisticas$DG_var)]

# Criando o boxplot
box <- ggplot(resultados, aes( interaction(nAuto, nCruz), y = D)) +
  geom_boxplot() +
  geom_hline(yintercept = dg_mean_value, color = "red", linetype = "dashed", size = 1) +  # Linha vermelha
  labs(
    title = "Estimated inbreeding depression",
    x = "",
    y = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotacionando os rótulos do eixo X para melhor leitura
    plot.title = element_text(hjust = 0.5)  # Centralizando o título
  )

print(box)

# Encontrando o valor de DG2_mean com a menor variância
dg2_mean_value <- estatisticas$DG2_mean[which.min(estatisticas$DG2_var)]

# Criando o boxplot
boxD2 <- ggplot(resultados, aes( interaction(nAuto, nCruz), y = D2)) +
  geom_boxplot() +
  geom_hline(yintercept = dg2_mean_value, color = "red", linetype = "dashed", size = 1) +  # Linha vermelha
  labs(
    title = "Estimated inbreeding depression D2",
    x = "",
    y = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotacionando os rótulos do eixo X para melhor leitura
    plot.title = element_text(hjust = 0.5)  # Centralizando o título
  )

print(boxD2)

#------------------------------------------------------------------------------#

library(ggplot2)
library(forcats)



rmse <- ggplot(estatisticas, aes(x =interaction(nAuto, nCruz), y = RMSE)) +
  geom_bar(stat = "identity", fill = "white", color = "black", size = 0.8) +  # Aumentando a espessura do contorno
  labs(title = "RMSE", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotacionando os rótulos do eixo X para melhor leitura
    plot.title = element_text(hjust = 0.5)  # Centralizando o título
  )
# Exibindo o gráfico
print(rmse)

rmse2 <- ggplot(estatisticas, aes(x =interaction(nAuto, nCruz), y = RMSE2)) +
  geom_bar(stat = "identity", fill = "white", color = "black", size = 0.8) +  # Aumentando a espessura do contorno
  labs(title = "RMSE2", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotacionando os rótulos do eixo X para melhor leitura
    plot.title = element_text(hjust = 0.5)  # Centralizando o título
  )
# Exibindo o gráfico
print(rmse2)
