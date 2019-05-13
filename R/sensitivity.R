# This script solves the question:
# Are there genes where (a lot of or certain) interactions are 
# only annotated for the main gene identifier, but not annotated for any 
# of its reported proteoforms, while there are proteoforms reported?  

library(ggplot2)
library(cowplot)
library(plyr)

hitsAllProteoforms <- read.csv("../resources/sensitivity/hitsByAllProteoformsForProteinsWithMultipleProteoforms.csv", header = T, fileEncoding="UTF-8-BOM")
 
# --------------- Plot the density of reactions and pathways mapped by unmodified proteoforms of proteins that have multiple proteoforms

hitsBasicProteoforms <- read.csv("../resources/sensitivity/hitsByBasicProteoformsForProteinsWithMultipleProteoforms.csv", header = T, fileEncoding="UTF-8-BOM")

hits <- merge(hitsAllProteoforms, hitsBasicProteoforms, by = "protein")
colnames(hits) <- c("protein", "accession_reactions", "accession_pathways", "basicProteoforms_reactions", "basicProteoforms_pathways")
hits$reactions_share <- (hits$basicProteoforms_reactions*100)/hits$accession_reactions
hits$pathways_share <- (hits$basicProteoforms_pathways*100)/hits$accession_pathways

p_basic_pathways <- ggplot(hits) +
    geom_density(aes(x=pathways_share),color="darkblue", fill="lightblue") +
    geom_vline(aes(xintercept=mean(pathways_share)), color="blue", linetype="dashed", size=1)
p_basic_pathways

p_basic_reactions <- ggplot(hits) +
    geom_density(aes(x=reactions_share),color="darkgreen", fill="lightgreen") +
    geom_vline(aes(xintercept=mean(reactions_share)), color="green", linetype="dashed", size=1)
p_basic_reactions

# --------------- Plot the density of reactions and pathways mapped by modified proteoforms of proteins that have multiple proteoforms

hitsModifiedProteoforms <- read.csv("../resources/sensitivity/hitsByModifiedProteoformsForProteinsWithMultipleProteoforms.csv", header = T, fileEncoding="UTF-8-BOM")

hits2 <- merge(hitsAllProteoforms, hitsModifiedProteoforms, by = "protein")
colnames(hits2) <- c("protein", "accession_reactions", "accession_pathways", "modifiedProteoforms_reactions", "modifiedProteoforms_pathways")
hits2$reactions_share <- (hits2$modifiedProteoforms_reactions*100)/hits2$accession_reactions
hits2$pathways_share <- (hits2$modifiedProteoforms_pathways*100)/hits2$accession_pathways

p_modified_pahtways <- ggplot(hits2) +
    geom_density(aes(x=pathways_share),color="darkblue", fill="lightblue") +
    geom_vline(aes(xintercept=mean(pathways_share)), color="blue", linetype="dashed", size=1)
p_modified_pahtways
p_modified_reactions <- ggplot(hits2) +
    geom_density(aes(x=reactions_share),color="darkgreen", fill="lightgreen") +
    geom_vline(aes(xintercept=mean(reactions_share)), color="green", linetype="dashed", size=1)
p_modified_reactions


# -----------------------

p <- plot_grid(p_basic_reactions, p_basic_pathways, p_modified_reactions, p_modified_pahtways, labels = c("A", "B", "C", "D"), ncol = 2)
save_plot("plot2by2.png", p,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
)

#---------------- Plot the density of number of proteoforms

proteoformsPerProtein <- read.csv("../resources/sensitivity/numberOfProteoformsPerProtein.csv")

p <- ggplot(proteoformsPerProtein) +
    geom_density(aes(x=num_proteoforms),color="darkblue", fill="maroon") +
    geom_vline(aes(xintercept=mean(num_proteoforms)), color="red", linetype="dashed", size=1) +
    xlim(0, 10) + ggtitle("Number of proteoforms for each protein")
p

proteoformsPerProtein <- read.csv("../resources/sensitivity/numberOfProteoformsPerProteinWithMultipleProteoforms.csv")

p <- ggplot(proteoformsPerProtein) +
    geom_density(aes(x=num_proteoforms),color="darkblue", fill="maroon") +
    geom_vline(aes(xintercept=mean(num_proteoforms)), color="red", linetype="dashed", size=1) +
    xlim(0, 10) + ggtitle("Number of proteoforms for each protein with multiple proteoforms")
p

# --------------- Comparison of proteins with 2, 3, 4 and 5 proteoforms, 
# by share of reactions and pathways mapped with modified proteoforms from the total possible
# reactions and pathways mapped with the accession

hitsProteins2ProteoformsByAllProteoforms <- read.csv("../resources/sensitivity/hitsByAllProteoformsForProteinsWith2Proteoforms.csv", header = T, fileEncoding="UTF-8-BOM")
hitsProteins2ProteoformsByModifiedProteoforms <- read.csv("../resources/sensitivity/hitsByModifiedProteoformsForProteinsWith2Proteoforms.csv", header = T, fileEncoding="UTF-8-BOM")
hitsProteins3ProteoformsByAllProteoforms <- read.csv("../resources/sensitivity/hitsByAllProteoformsForProteinsWith3Proteoforms.csv", header = T, fileEncoding="UTF-8-BOM")
hitsProteins3ProteoformsByModifiedProteoforms <- read.csv("../resources/sensitivity/hitsByModifiedProteoformsForProteinsWith3Proteoforms.csv", header = T, fileEncoding="UTF-8-BOM")
hitsProteins4ProteoformsByAllProteoforms <- read.csv("../resources/sensitivity/hitsByAllProteoformsForProteinsWith4Proteoforms.csv", header = T, fileEncoding="UTF-8-BOM")
hitsProteins4ProteoformsByModifiedProteoforms <- read.csv("../resources/sensitivity/hitsByModifiedProteoformsForProteinsWith4Proteoforms.csv", header = T, fileEncoding="UTF-8-BOM")
hitsProteins5ProteoformsByAllProteoforms <- read.csv("../resources/sensitivity/hitsByAllProteoformsForProteinsWith5Proteoforms.csv", header = T, fileEncoding="UTF-8-BOM")
hitsProteins5ProteoformsByModifiedProteoforms <- read.csv("../resources/sensitivity/hitsByModifiedProteoformsForProteinsWith5Proteoforms.csv", header = T, fileEncoding="UTF-8-BOM")

hitsProteins2Proteoforms <- merge(hitsProteins2ProteoformsByAllProteoforms, hitsProteins2ProteoformsByModifiedProteoforms, by = "protein")
hitsProteins3Proteoforms <- merge(hitsProteins3ProteoformsByAllProteoforms, hitsProteins3ProteoformsByModifiedProteoforms, by = "protein")
hitsProteins4Proteoforms <- merge(hitsProteins4ProteoformsByAllProteoforms, hitsProteins4ProteoformsByModifiedProteoforms, by = "protein")
hitsProteins5Proteoforms <- merge(hitsProteins5ProteoformsByAllProteoforms, hitsProteins5ProteoformsByModifiedProteoforms, by = "protein")

hitsProteins2Proteoforms$num_proteoforms <- 2
hitsProteins3Proteoforms$num_proteoforms <- 3
hitsProteins4Proteoforms$num_proteoforms <- 4
hitsProteins5Proteoforms$num_proteoforms <- 5

hitsProteinsMultipleProteoforms <- rbind(hitsProteins2Proteoforms, hitsProteins3Proteoforms, hitsProteins4Proteoforms, hitsProteins5Proteoforms)
hitsProteinsMultipleProteoforms$num_proteoforms <- as.factor(hitsProteinsMultipleProteoforms$num_proteoforms)

colnames(hitsProteinsMultipleProteoforms) <- c("protein", "accession_reactions", "accession_pathways", "modifiedProteoforms_reactions", "modifiedProteoforms_pathways", "num_proteoforms")
hitsProteinsMultipleProteoforms$reactions_share <- (hitsProteinsMultipleProteoforms$modifiedProteoforms_reactions*100)/hitsProteinsMultipleProteoforms$accession_reactions
hitsProteinsMultipleProteoforms$pathways_share <- (hitsProteinsMultipleProteoforms$modifiedProteoforms_pathways*100)/hitsProteinsMultipleProteoforms$accession_pathways

mu <- ddply(hitsProteinsMultipleProteoforms, "num_proteoforms", summarise, grp.mean=mean(reactions_share))
head(mu)

p_reactions <- ggplot(hitsProteinsMultipleProteoforms, aes(x=reactions_share, fill=num_proteoforms)) +
    geom_density(alpha=0.4) +
    geom_vline(data=mu, aes(xintercept=grp.mean, color=num_proteoforms), linetype="dashed")
p_reactions

mu <- ddply(hitsProteinsMultipleProteoforms, "num_proteoforms", summarise, grp.mean=mean(pathways_share))
head(mu)

p_pathways <- ggplot(hitsProteinsMultipleProteoforms, aes(x=pathways_share, fill=num_proteoforms)) +
    geom_density(alpha=0.4) +
    geom_vline(data=mu, aes(xintercept=grp.mean, color=num_proteoforms), linetype="dashed")
p_pathways

grid_comparison <- plot_grid(p_reactions, p_pathways, labels = c("A", "B"), nrow = 2)
grid_comparison
save_plot("../docs/mappingMultipleProteoforms.png", grid_comparison,
          ncol = 1, nrow = 2,
          base_aspect_ratio = 1.1
)
