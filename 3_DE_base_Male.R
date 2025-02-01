rm(list = ls()); gc()

library(tidyverse)
library(readxl)

# Phenotypic data ---------------------------------------------------------

phenoF = read.table("outputs/pheno_data_F.txt")
phenoF[1:5] = lapply(phenoF[1:5], as.factor)
str(phenoF)

ind_per_family = aggregate(IND ~ Pai, data = phenoF, FUN = length)
names(ind_per_family)[2] = "NUM_IND"


# Means per Family --------------------------------------------------------

# Count individuals per family and CROSS class
count_ind = aggregate(IND ~ Pai + CLASSE, data = phenoF, FUN = length)
names(count_ind)[3] = "COUNT"

# Calculate the mean TRAIT value for each FAMILY and CROSS class
mean_per_family = aggregate(dap ~ Pai + CLASSE, data = phenoF, FUN = function(x) mean(x, na.rm = TRUE))


# Reshape the data to wide format for easy subtraction
trait_means_wide = mean_per_family %>%
  pivot_wider(names_from = CLASSE, values_from = dap, values_fill = 0)

trait_means_wide[trait_means_wide==0]=NA

# Compute the subtractions
trait_means_wide = trait_means_wide %>%
  mutate(
    CRUZ_menos_AUTO = CRUZ - AUTO
  )

trait_means_wide = merge(trait_means_wide, ind_per_family, by = c("Pai"))
trait_means_wide$MEAN_ALL = rowMeans(trait_means_wide[c(2,3)], na.rm = TRUE)


# s based on F ------------------------------------------------------------

dim(phenoF)
str(phenoF)

mean_endogamy = phenoF %>%
  group_by(Pai) %>%
  summarise(ME = mean(F_roh))

mean(mean_endogamy$ME)
#0.2366702

# Reorder the Family factor based on F (descending order)
mean_endogamy$Pai = as.factor(mean_endogamy$Pai)
mean_endogamy_order = mean_endogamy[order(mean_endogamy$ME, decreasing = T), ]

#levels(df1$Pai) == levels(mean_endogamy_order$Pai)

# Ensure df2$COLUMN has the desired order as a factor
mean_endogamy_order$Pai = factor(mean_endogamy_order$Pai, 
                                 levels = unique(mean_endogamy_order$Pai))

write.table(mean_endogamy_order, "outputs/F_by_Male.txt")


F_by_fam = mean_endogamy_order
str(F_by_fam)

means = trait_means_wide %>% 
  full_join(F_by_fam, by="Pai")

# s = 2F/(1+F)
means$s = (2*means$ME)/(1+means$ME)

write.table(means, "outputs/Male_phenotipic_means_classes.txt")

# DE = (100/media)((M-S)/(1/(2-s)-(s/(2-s))))
##########################################################
means$DE_CRUZ = (100/means$MEAN_ALL)* (means$CRUZ_menos_AUTO/((1/(2-means$s))-(means$s/(2-means$s))))

means = means %>% rename(F_HBD = ME)

write.table(means, "outputs/DE_per_Male.txt")

mean(means$DE_CRUZ, na.rm = TRUE)
#91.76274


mean(means$AUTO, na.rm = TRUE)
#9.732886
mean(means$CRUZ, na.rm = TRUE)
#13.87962


ID = ((13.87962-9.732886)/13.87962)*100
#ID  #29.87642

0.2987642*12.65
#3.779367

mean(means$F_HBD, na.rm = TRUE)
# mean FHBD = 0.2366702


means2 = means[!is.na(means$CRUZ_menos_AUTO), ]
#######################################################


mean(means2$AUTO, na.rm = TRUE)
#9.463966
mean(means2$CRUZ, na.rm = TRUE)
#13.4356

# Regression DE on F ------------------------------------------------------


means = means[!is.na(means$CRUZ_menos_AUTO), ] ; dim(means)

# Calculate R-squared for each family
r2_data = means %>%
  summarise(
    R2 = summary(lm(DE_CRUZ ~ s))$r.squared
  )


# Create the faceted plot
p4 = ggplot(means, aes(x = s, y = DE_CRUZ)) +
  geom_point(color =  "#2c6e49", alpha = 0.7, size=3) +            # Scatter plot
  geom_smooth(method = "lm", color = "#D68c45", se = TRUE, 
              fill = "darkgray", alpha = 0.2) +         #Regression line    
  geom_text(data = r2_data, 
            aes(x = 0.5, y = 150, label = paste0("R² = ", round(R2, 2))), 
            inherit.aes = FALSE, size = 4, color = "black") +  # Add R2
  labs(
    title = "Regression of Endogamy Depression of Parentals on s",
    x = "s",
    y = "Endogamy Depression"
  ) +
  theme_minimal() ;p4



# Calculate R-squared for each family
r2_data = means %>%
  summarise(
    R2 = summary(lm(DE_CRUZ ~ F_HBD))$r.squared
  )

# Create the faceted plot
p5 = ggplot(means, aes(x = F_HBD, y = DE_CRUZ)) +
  geom_point(color = "#2c6e49", alpha = 0.7, size=3) +            # Scatter plot
  geom_smooth(method = "lm", color = "#D68c45", se = TRUE, 
              fill = "darkgray", alpha = 0.2) +         #Regression line    
  geom_text(data = r2_data, 
            aes(x = 0.3, y = 200, label = paste0("R² = ", round(R2, 2))), 
            inherit.aes = FALSE, size = 5, color = "black") +  # Add R2
  ggrepel::geom_text_repel(data = means,
    aes(label = Pai),
    size = 3,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.4, "lines")) +
  labs(
    #title = "Family based Regression of Endogamy Depression of Outbread on F_HBD",
    x = "F",
    y = "Endogamy Depression"
  ) +
  theme_minimal() +
  theme(legend.text = element_text(face = c("italic")),
        axis.text.y = element_text(hjust = 1,color="black", size = 10),
        axis.text.x = element_text(angle = 0, hjust = 1, 
                                   color =  "black", size = 10)) ;p5


library(patchwork)
p4+p5


# Grafico de Barras -------------------------------------------------------

means_all = read.table("outputs/DE_per_Male.txt")

names(means_all)
means_all_de = means_all %>% select(Pai, CRUZ, AUTO, F_HBD)
means_all_de = means_all_de[means_all_de$Pai!=0, ]
means_all_de = means_all_de[order(means_all_de$F_HBD, decreasing=TRUE), ]
dim(means_all_de)

means_all_de_long = means_all_de %>%
  pivot_longer(cols = c(CRUZ, AUTO))

# Carregar ggplot2 (caso necessário)
library(ggplot2)

means_all_de_long$name = factor(means_all_de_long$name, levels = c("CRUZ", "AUTO"))

means_all_de_long$Pai = factor(means_all_de_long$Pai, levels = means_all_de$Pai)

# Criar o gráfico de barras
p6 = ggplot(means_all_de_long, aes(x = Pai, y = value, fill = name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = round(F_HBD, 2), y = 19),
            color = "black", 
            size = 3, angle = 45) +
  scale_y_continuous(breaks = seq(0, 18, by=3))+
  theme_bw()+
  labs(
    #title = "Gráfico de Barras para Indivíduos e Características",
    x = "Parents",
    y = "DBH (cm)"
  ) +
  scale_fill_manual(values = c("#2c6e49", "#D68c45")) +  # Cores personalizadas
  theme_minimal() +
  theme(legend.text = element_text(face = c("italic")),
        axis.text.y = element_text(hjust = 1,color="black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, 
                                   color =  "black", size = 10) 
  );p6


library(patchwork)
p6 + p5


# Criar o gráfico de barras da DE baseado no fenotipo

means_all_red = means_all[!is.na(means_all$CRUZ_menos_AUTO), ]
means_all_red = means_all_red[order(means_all_red$F_HBD, decreasing = TRUE), ]
means_all_red$Pai= factor(means_all_red$Pai, levels = means_all_red$Pai)

p7 = ggplot(means_all_red, aes(x = Pai, y = CRUZ_menos_AUTO, fill = "Pai")) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = round(F_HBD, 2), y = 8),
            color = "black", 
            size = 3, angle = 45) +
  scale_y_continuous(breaks = seq(0, 8, by=2))+
  theme_bw()+
  labs(
    #title = "Gráfico de Barras para Indivíduos e Características",
    x = "Parents",
    y = "Reduction in DBH (cm)"
  ) +
  scale_fill_manual(values = c("#D68c45")) +  # Cores personalizadas
  theme_minimal() +
  theme(legend.text = element_text(face = c("italic")),
        axis.text.y = element_text(hjust = 1,color="black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, 
                                   color =  "black", size = 10) 
  );p7



# Create the faceted plot

# Calculate R-squared for each family
r2_data = means_all_red %>%
  summarise(
    R2 = summary(lm(CRUZ_menos_AUTO ~ F_HBD))$r.squared
  )

p8 = ggplot(means_all_red, aes(x = F_HBD, y = CRUZ_menos_AUTO)) +
  geom_point(color = "#2c6e49", alpha = 0.7, size=3) +            # Scatter plot
  geom_smooth(method = "lm", color = "#D68c45", se = TRUE, 
              fill = "darkgray", alpha = 0.2) +         #Regression line    
  geom_text(data = r2_data, 
            aes(x = 0.35, y = 10, label = paste0("R² = ", round(R2, 2))), 
            inherit.aes = FALSE, size = 5, color = "black") +  # Add R2
  ggrepel::geom_text_repel(data = means,
                           aes(label = Pai),
                           size = 3,
                           box.padding = unit(0.4, "lines"),
                           point.padding = unit(0.4, "lines")) +
  labs(
    #title = "Family based Regression of Endogamy Depression of Outbread on F_HBD",
    x = "F",
    y = "Reduction in DBH (cm)"
  ) +
  theme_minimal() +
  theme(legend.text = element_text(face = c("italic")),
        axis.text.y = element_text(hjust = 1,color="black", size = 10),
        axis.text.x = element_text(angle = 0, hjust = 1, 
                                   color =  "black", size = 10)) ;p8



library(patchwork)
(p6 + p5) / (p7 + p8)






