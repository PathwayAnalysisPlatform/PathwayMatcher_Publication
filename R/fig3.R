# 
# This script plots the interactors of a gene in a gene- and proteoform-centric paradigm.
#
startTimeAll <- proc.time()


# Libraries

library(ggplot2)
library(ggrepel)
library(scico)
library(graphlayouts)
library(dplyr)
library(purrr)
library(igraph)


# Parameters

set.seed(240419)

## Input files

proteoformEdgesFile <- "resources/network_1.8.1/proteoformInternalEdges.tsv.gz"
accessionsEdgesFile <- "resources/network_1.8.1/proteinInternalEdges.tsv.gz"

## Plot theme

theme_set(theme_bw(base_size = 11))


## Colors

palette <- 'cork'
accessionColor <- scico(n = 1, begin = 0.15, end = 0.15, palette = palette)
proteoformColor <- scico(n = 1, begin = 0.85, end = 0.85, palette = palette)


## Modification names

modNameDFs <- data.frame(
    accession = c(
        "00046", 
        "00047", 
        "00064", 
        "00076",
        "00084", 
        "00078", 
        "00085", 
        "01148",
        "01149"
    ),
    fullName = c(
        "O-phospho-L-serine", 
        "O-phospho-L-threonine", 
        "N6-acetyl-L-lysine", 
        "symmetric dimethyl-L-arginine",
        "N6,N6-dimethyl-L-lysine",
        "omega-N-methyl-L-arginine",
        "N6-methyl-L-lysine", 
        "ubiquitinylated lysine",
        "sumoylated lysine"
    ),
    shortName = c(
        "pS", 
        "pT", 
        "aceK", 
        "dimethR", 
        "dimethK", 
        "methR", 
        "methK", 
        "ubiK",
        "sumoK"
    ),
    stringsAsFactors = F
)




# Functions

#' Returns the accession corresponding to a proteoform.
#' 
#' @param proteoform the proteoform identifier
#' 
#' @return the protein accession of the proteoform
getAccession <- function(proteoform) {
    
    accession <- substr(proteoform, start = 1, stop = regexpr(proteoform, pattern = ";")-1)
    
    indexDash <- regexpr(accession, pattern = "-")
    
    if (indexDash > 1) {
        
        accession <- substr(accession, start = 1, stop = indexDash - 1)
        
    }
    
    return(accession)
    
}

#' Returns the name corresponding to a proteoform.
#' 
#' @param proteoform the proteoform identifier
#' 
#' @return the name of the proteoform
getLabel <- function(proteoform) {
    
    indexSemicolon <- regexpr(proteoform, pattern = ";")
    
    accession <- substr(proteoform, start = 1, stop = indexSemicolon - 1)
    indexDash <- regexpr(accession, pattern = "-")
    
    if (indexDash > 1) {
        
     stop("Isoform not implemented.")
    
    }
    
    if (indexSemicolon == nchar(proteoform)) {
        return("Canonical")
    }
    
    modString <- substr(proteoform, start = indexSemicolon + 1, stop = nchar(proteoform))
    
    mods <- strsplit(modString, split = ",")[[1]]
    modNames <- character(length(mods))
    
    for (i in 1:length(mods)) {
        
        mod <- mods[i]
        indexColon <- regexpr(mod, pattern = ":")
        modAccession <- substr(mod, start = 1, stop = indexColon - 1)
        
        if (! modAccession %in% modNameDFs$accession) {
            stop(paste0("Name for modification ", modAccession, " not found."))
        }
        
        modName <- modNameDFs$shortName[modNameDFs$accession == modAccession]
        modSite <- substr(mod, start = indexColon + 1, stop = nchar(mod))
        
        if (modSite != "null") {
            
            modName <- paste0(modName, modSite)
            
        }
        
        modNames[i] <- modName
        
    }
    
    proteoformName <- paste(modNames, collapse = " ")
    
    return(proteoformName)
    
}

#' Returns the modifications found in a proteoform.
#' 
#' @param proteoform the proteoform identifier
#' 
#' @return the modifications found in a proteoform
getModifications <- function(proteoform) {
    
    modifications <- substr(proteoform, start = regexpr(proteoform, pattern = ";")+1, stop = 100000)
    split1 <- strsplit(x = modifications, ",")
    modifications <- substr(split1[[1]], start = 1, stop = regexpr(split1[[1]], pattern = ":")-1)
    
    return(unique(modifications))
    
}


# Main script
 
## Load data

print(paste(Sys.time(), " Loading data", sep = ""))

edgesProteoforms <- read.table(proteoformEdgesFile, header = T, sep = "\t", quote = "", comment.char = "", stringsAsFactors = F)
edgesAccessions <- read.table(accessionsEdgesFile, header = T, sep = "\t", quote = "", comment.char = "", stringsAsFactors = F)


## Format

edgesProteoforms <- edgesProteoforms[, c("id1", "id2")]
edgesAccessions <- edgesAccessions[, c("id1", "id2")]

accession1 <- sapply(X = edgesProteoforms$id1, FUN = getAccession)
accession2 <- sapply(X = edgesProteoforms$id2, FUN = getAccession)


## Proteoform and protein occurrence

allProteoforms <- unique(c(edgesProteoforms$id1, edgesProteoforms$id2))
allProteoformsAccessions <- sapply(X = allProteoforms, FUN = getAccession)

proteoformsPerAccession <- as.data.frame(table(allProteoformsAccessions))


## Select edges interacting with protein of interest

targetAccession <- "P04637"

proteoformsExample <- unique(c(edgesProteoforms$id1[accession1 == targetAccession], edgesProteoforms$id2[accession2 == targetAccession]))
proteoformsExampleModifications <- c()

for (proteoform in proteoformsExample) {
    
    proteoformsExampleModifications <- c(proteoformsExampleModifications, getModifications(proteoform))
    
}

proteoformsExampleModifications <- unique(proteoformsExampleModifications)

edgesProteoformsExample <- edgesProteoforms[accession1 == targetAccession | accession2 == targetAccession, ]
edgesAccessionsExample <- edgesAccessions[edgesAccessions$id1 == targetAccession | edgesAccessions$id2 == targetAccession, ]


## Make graph

graphProteoforms <- graph_from_data_frame(edgesProteoformsExample)
graphAccessions <- graph_from_data_frame(edgesAccessionsExample)

graphProteoforms <- igraph::simplify(graphProteoforms, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")
graphAccessions <- igraph::simplify(graphAccessions, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")


## Plot graphs

layout <- layout_igraph_stress(graphAccessions)

verticesDF <- data.frame(id = V(graphAccessions)$name, x = layout[, 1], y = layout[, 2], stringsAsFactors = F)
verticesDF %>%
    mutate(
        target = ifelse(id == targetAccession, targetAccession, "Other")
    ) -> verticesDF
x <- verticesDF$x
names(x) <- verticesDF$id
y <- verticesDF$y
names(y) <- verticesDF$id

verticesDF %>% 
    filter(
        target == targetAccession
    ) -> targetVertices

edgesList <- get.edgelist(graphAccessions)
edgesDF <- data.frame(name1 = edgesList[, 1], name2 = edgesList[, 2], stringsAsFactors = F)
edgesDF %>%
    mutate(
        x1 = x[name1],
        y1 = y[name1],
        x2 = x[name2],
        y2 = y[name2]
    ) -> edgesDF

graphPlot <- ggplot() + 
    geom_segment(data = edgesDF, aes(x = x1, y = y1, xend = x2, yend = y2), col = accessionColor, alpha = 0.2, size = 0.5) + 
    geom_point(data = verticesDF, aes(x = x, y = y, col = target, size = target, alpha = target)) + 
    geom_label(data = targetVertices, aes(x = x, y = y, label = "italic(TP53)"), col = "red", parse = T) + 
    scale_color_manual(values = c(accessionColor, "red")) + 
    scale_size_manual(values = c(1, 2)) + 
    scale_alpha_manual(values = c(0.5, 1)) +
    ggtitle("Gene-centric") + 
    theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "white"),
        plot.title = element_text(hjust = 0.5)
    )


png(paste0("docs/fig_3A.png"), height = 12, width = 12, units = "cm", res = 600)
plot(graphPlot)
dummy <- dev.off()

layout <- layout_igraph_backbone(graphProteoforms)

verticesDF <- data.frame(id = V(graphProteoforms)$name, x = layout[, 1], y = layout[, 2], stringsAsFactors = F)
verticesDF %>%
    mutate(
        accession = map(id, getAccession),
        target = ifelse(accession == targetAccession, targetAccession, "Other")
    ) -> verticesDF
x <- verticesDF$x
names(x) <- verticesDF$id
y <- verticesDF$y
names(y) <- verticesDF$id

verticesDF %>% 
    filter(
        target == targetAccession
    ) %>%
    mutate(
        label = map(id, getLabel),
        number = row_number()
    ) -> targetVertices

edgesList <- get.edgelist(graphProteoforms)
edgesDF <- data.frame(name1 = edgesList[, 1], name2 = edgesList[, 2], stringsAsFactors = F)
edgesDF %>%
    mutate(
        x1 = x[name1],
        y1 = y[name1],
        x2 = x[name2],
        y2 = y[name2],
        accession1 = map(name1, getAccession),
        accession2 = map(name2, getAccession),
        target = ifelse(accession1 == targetAccession & accession2 == targetAccession, targetAccession, "Other")
    ) %>%
    arrange(target) -> edgesDF

graphPlot <- ggplot() + 
    geom_segment(data = edgesDF, aes(x = x1, y = y1, xend = x2, yend = y2, col = target, size = target), alpha = 0.2) + 
    geom_point(data = verticesDF, aes(x = x, y = y, col = target, size = target, alpha = target)) + 
    geom_label_repel(data = targetVertices, aes(x = x, y = y, label = number), col = "red") + 
    scale_color_manual(values = c(proteoformColor, "red")) + 
    scale_size_manual(values = c(1, 2)) + 
    scale_alpha_manual(values = c(0.5, 1)) +
    ggtitle("Proteoform-centric") + 
    theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "white"),
        plot.title = element_text(hjust = 0.5)
    )


png(paste0("docs/fig_3B.png"), height = 12, width = 12, units = "cm", res = 600)
plot(graphPlot)
dummy <- dev.off()

