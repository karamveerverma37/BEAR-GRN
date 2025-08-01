setwd("/gpfs/Home/kmk7420/BEAR-GRN/BEAR-GRN/INFERRED.GRNS")
#create a empy list to add scRNA and scATAC data
K562 <- list()
RNA <- list()
RNA$CellOracle <- read.csv("K562/CellOracle/K562_human_filtered_out_E7.5_rep1_final_GRN.csv")
RNA$FigR <- read.csv("K562/FigR/K562_filtered_network.csv")

#add RNA data
K562$RNA <- RNA
K562
#save.image(file="Inferred_GRNs.RData")
