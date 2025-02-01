# Figuras

remove(list = ls()) ;gc()

library(readxl); library(tidyverse)
# Composiçao de cada familia e endogamia

composicao = readxl::read_excel("data/composicao_especie.xlsx")
str(composicao)

composicao[2:9] = lapply(composicao[2:9], function(x) round(as.numeric(x), digits = 4))
str(composicao)

comp_male = composicao[ , c(1:6, 9)]
str(comp_male)

comp_male= comp_male[order(comp_male$F_Male, decreasing = TRUE), ]

comp_male = comp_male[!is.na(comp_male$F_Male),]
dim(comp_male)
comp_male = comp_male %>% rename("Unknown" = "Desconhecido")

male_long <- comp_male %>% 
  pivot_longer(
    cols = `E. urophylla`:`Unknown`, 
    names_to = "specie",
    values_to = "composition"
  )


male_long$ID_GEN = factor(male_long$ID_GEN, levels = comp_male$ID_GEN)





male_long  %>%
  ggplot(aes(x = factor(ID_GEN), y = composition, fill = fct_reorder(specie, composition))) +
  geom_col(position = "stack") +
  geom_text(aes(label = round(F_Male, 2), y = 1.05),
            color = "black", 
            size = 4) +
  scale_fill_manual(values = c("#0D353F","#9D683B","gray","#CBB57B", "#017979"),
                    name = "Species")+ 
  scale_y_continuous(breaks = seq(0, 1))+
  theme_bw()+
  ylab("Expected Proportion of Genetic Contribution (%)")+
  xlab("Selected Families")+
  theme(legend.text = element_text(face = c("italic")),
        axis.text.y = element_text(hjust = 1,color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, 
                                   color =  "black") 
  )


