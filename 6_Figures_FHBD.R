# save.image("save_images/image_all_4_12.RData")
load("save_images/image_all_4_12.RData")

# Inbreeding --------------------------------------------------------------

x <- 1- inbreading_results@realized[,10]
hist(x,nc=20,main="",xlab="Inbreeding coefficient",col='coral2')

#ate o T=64
x <- t(apply(inbreading_results@realized[,3:6],1,cumsum))
hist(x[,4],nc=20,main="",xlab="Inbreeding coefficient (T = 64)",col='tomato')

# Fig1 - 20 mais endogamicos  -------------------------------------------------------------------
a = data.frame(inbreading_results@realized[ ,1:9])
dim(a)
a$soma = rowSums(a)
a$IND = inbreading_results@ids

a_order = a[order(a$soma, decreasing = T), ]
high_F_ind = a_order[1:20, ]

toplot = list(high_F_ind$IND)

library(svglite)
svglite("results/20_HighF.svg")
p1 = zooplot_partitioning(list(HighF=inbreading_results),
                     toplot = toplot, plotids= TRUE,
                     nonhbd = FALSE,
                     vertical = TRUE); p1
dev.off()


# Fig2 - 20 maiores g  -------------------------------------------------------------------

a = data.frame(inbreading_results@realized[ ,1:9])
dim(a)
a$soma = rowSums(a)
a$IND = inbreading_results@ids
a_order = a[order(a$soma, decreasing = T), ]

BLUPS = read.table("results/BLUPS-AD-GBLUP.txt")
BLUPS = BLUPS[order(BLUPS$g, decreasing = TRUE), ]

high_g_ind = BLUPS[1:20, ]
toplot2 = list(rownames(high_g_ind))

# Extrair os nomes
toplot2_df = data.frame(IND = unlist(toplot2))

# Extrair as últimas 6 letras de cada coluna
a$last_6 <- substr(a$IND, nchar(a$IND) - 5, nchar(a$IND))
toplot2_df$last_6 <- substr(toplot2_df$IND, nchar(toplot2_df$IND) - 5, nchar(toplot2_df$IND))

# Encontrar as linhas onde as últimas 6 letras correspondem
matching_rows <- merge(a, toplot2_df, by = "last_6")
matching_rows = matching_rows[1:20, ]
dim(matching_rows)

length( unique(matching_rows) %in% unique(inbreading_results@ids))
toplot2 = list(matching_rows$IND.x)

svglite("results/20_High_G_FHBD.svg")
p2 = zooplot_partitioning(list(HighF=inbreading_results),
                     toplot = toplot2, plotids= TRUE,
                     nonhbd = FALSE,
                     vertical = TRUE); p2
dev.off()


