# 
# This script extracts and plots the degree distribution difference between proteins and proteoforms.
#
startTimeAll <- proc.time()


# Libraries

library(ggplot2)
library(ggrepel)
library(igraph)
library(scico)
library(gtable)
library(grid)
library(dplyr)
library(purrr)


# Parameters

## Network files

proteoformEdgesFile <- "resources/network_1.8.1/proteoformInternalEdges.tsv.gz"
accessionsEdgesFile <- "resources/network_1.8.1/proteinInternalEdges.tsv.gz"

## Colors

palette <- 'cork'
accessionColor <- scico(n = 1, begin = 0.15, end = 0.15, palette = palette)
proteoformColor <- scico(n = 1, begin = 0.85, end = 0.85, palette = palette)


# Functions

#' Returns the accession without isoform number.
#' 
#' @param accession the accession
#' @return the accession without isoform number
removeIsoform <- function(accession) {
  
  indexDash <- regexpr(accession, pattern = "-")
  
  if (indexDash > 1) {
    
    accession <- substr(accession, start = 1, stop = indexDash - 1)
    
  }
  
  return(accession)
  
}

#' Returns the accession corresponding to a proteoform.
#' 
#' @param proteoform the proteoform identifier
#' @return the protein accession of the proteoform
getAccession <- function(proteoform) {
  
  accession <- substr(proteoform, start = 1, stop = regexpr(proteoform, pattern = ";")-1)
  
  accession <- removeIsoform(accession)
  
  return(accession)
  
}


#' Returns the accession corresponding to a proteoform.
#' 
#' @param proteoform the proteoform identifier
#' @return the protein accession of the proteoform
getGenes <- function(accession) {
  
  i <- which(accessionsMapping$accession == accession)
  
  return(c(accessionsMapping$id[i]))
  
}


#' Returns the degree of the protein with the given accession.
#' 
#' @param accession the protein accession
#' @param accessionNames the accession names corresponding to the degrees
#' @param accessionDegrees the degrees of the protein accessions
#' @return the degree of the protein accession
getDegreeFromAccession <- function(accession, accessionNames, accessionDegrees) {
  
  i <- which(accessionNames == accession)
  
  if (length(i) == 0) {
    
    return(0);
    
  }
  
  return(accessionDegrees[i])
  
}


#' Returns the degree of the protein with the given accession when mapping back to the genes
#' 
#' @param accession the protein accession
#' @param geneNames the gene names corresponding to the degrees
#' @param geneDegrees the degrees of the gene names
#' @return the degree of the protein accession
getDegreeFromGene <- function(accession, geneNames, geneDegrees) {
  
  degree <- 0
  
  geneMapping <- getGenes(accession)
  
  for (geneName in geneMapping) {
    
    i <- which(geneNames == geneName)
    
    if (length(i) > 0) {
      
      degree <- degree + geneDegrees[i]
      
    }
  }
  
  
  return(degree)
  
}


#' Returns a data frame containing log2 binned degree probabilities for the given degrees.
#' 
#' @param degrees the degrees
#' @return a data frame containing log2 binned degree probabilities for the given degrees
getDegreDF <- function(degrees) {
  
  bin <- 0
  degreesBinned <- c()
  ps <- c()
  
  maxDegree <- max(degrees)
  
  while(T) {
    
    bi <- 2 ^ bin
    biPlusOne <- 2 ^ (bin+1)
    degree <- mean(bi:(biPlusOne-1))
    
    pi <- sum(degrees >= bi & degrees < biPlusOne) / (biPlusOne - bi)
    
    if (pi > 0) {
      
      ps <- c(ps, pi)
      degreesBinned <- c(degreesBinned, degree)
      
    }
    
    if (biPlusOne > maxDegree) {
      
      break()
      
    }
    
    bin <- bin + 1
    
  }
  
  ps <- ps / length(degrees)
  
  degreeDF <- data.frame(degree = degreesBinned, p = ps)
  
  return(degreeDF)
  
}

#' Dodge points in the form of a sina plot.
#' 
#' @param value the value
#' 
#' @return the difference to add in x
getDodge <- function(ratio) {
  
  ld <- length(ratioDensity)
  
  if (sum(ratioDensity$x == ratio) == 0) {
    
    i1 <- max(which(ratioDensity$x < ratio))
    i2 <- min(which(ratioDensity$x > ratio))
    
    x1 <- ratioDensity$x[i1]
    x2 <- ratioDensity$x[i2]
    y1 <- ratioDensity$y[i1]
    y2 <- ratioDensity$y[i2]
    
    y <- y1 + ((y2 - y1) * (ratio - x1) / (x2 - x1))
    
    y <- y * runif(min = -1, max = 1, n = 1)
    
    return(y)
    
  } else if (sum(ratioDensity$x == ratio) == 1) {
    
    y <- ratioDensity$y[ratioDensity$x == ratio]
    
    y <- y * runif(min = -1, max = 1, n = 1)
    
    return(y)
    
  } else {
    stop("Multiple density x values matching ratio ", ratio)
  }
}


# Main script
 
## Load data

print(paste(Sys.time(), " Loading data", sep = ""))

edgesProteoforms <- read.table(proteoformEdgesFile, header = T, sep = "\t", quote = "", comment.char = "", stringsAsFactors = F)
edgesAccessions <- read.table(accessionsEdgesFile, header = T, sep = "\t", quote = "", comment.char = "", stringsAsFactors = F)


## Format

edgesProteoforms <- edgesProteoforms[, c("id1", "id2")]
edgesAccessions <- edgesAccessions[, c("id1", "id2")] %>%
  mutate(
    id1 = map(id1, removeIsoform),
    id2 = map(id2, removeIsoform)
  )


## Make graph

graphProteoforms <- graph_from_data_frame(edgesProteoforms)
graphAccessions <- graph_from_data_frame(edgesAccessions)

graphProteoforms <- simplify(graphProteoforms, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")
graphAccessions <- simplify(graphAccessions, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")


## Extract proteoforms

allProteoforms <- V(graphProteoforms)$name
lengths <- sapply(allProteoforms, FUN = nchar, USE.NAMES = F)
separatorI <- sapply(allProteoforms, FUN = regexpr, pattern = ';', USE.NAMES = F)
proteoforms <- allProteoforms[lengths > separatorI]

proteoform0 <- edgesProteoforms[! edgesProteoforms$id1 %in% proteoforms & ! edgesProteoforms$id2 %in% proteoforms, ]
proteoform1 <- edgesProteoforms[edgesProteoforms$id1 %in% proteoforms & ! edgesProteoforms$id2 %in% proteoforms
                                | ! edgesProteoforms$id1 %in% proteoforms & edgesProteoforms$id2 %in% proteoforms, ]
proteoform2 <- edgesProteoforms[edgesProteoforms$id1 %in% proteoforms & edgesProteoforms$id2 %in% proteoforms, ]

proteoform0Graph <- graph_from_data_frame(proteoform0)
proteoform1Graph <- graph_from_data_frame(proteoform1)
proteoform2Graph <- graph_from_data_frame(proteoform2)

proteoform0Graph <- simplify(proteoform0Graph, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")
proteoform1Graph <- simplify(proteoform1Graph, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")
proteoform2Graph <- simplify(proteoform2Graph, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")


## Get proteoform degree using other matchings

degreeProteoforms <- data.frame(accession = V(graphProteoforms)$name, degree = degree(graphProteoforms), type = "Proteoform", stringsAsFactors = F)

degreeProteoforms %>%
  mutate(
    modified = nchar(accession) > regexpr(accession, pattern = ';')
  ) %>%
  filter(
    modified & degree > 0
  ) %>%
  select(
    -modified
  ) -> degreeProteoforms

data.frame(
  accession = unique(sapply(X = degreeProteoforms$accession, FUN = getAccession, USE.NAMES = F)),
  type = "Accession",
  stringsAsFactors = F
) %>%
  rowwise() %>%
  mutate(
    degree = getDegreeFromAccession(accession, accessionNames = V(graphAccessions)$name, accessionDegrees = degree(graphAccessions))
  ) %>%
  filter(degree > 0) -> degreeAccessions

rbind(degreeProteoforms, degreeAccessions) %>%
  mutate(
    degreeLog = log10(degree)
  ) -> degreeDF

degreeProteoforms %>%
  select(-type) %>%
  rename(
    proteoform = accession,
    degreeProteoform = degree
  ) %>%
  mutate(
    accession = map_chr(proteoform, getAccession)
  ) -> matchingDF

matchingDF %>% left_join(
  degreeAccessions %>%
    select(-type) %>%
    rename(
      degreeAccession = degree
    ),
   by = "accession"
) %>%
  filter(degreeAccession > 0) %>%
  mutate(
    ratio = log10(degreeProteoform / degreeAccession)
  ) -> matchingDF


## Plot degree distributions and ratios

degreeDF$typeFactor <- factor(degreeDF$type)
degreeDF %>%
  group_by(typeFactor) %>%
  summarise(
    medianDegree = median(degreeLog)
  ) %>%
  mutate(
    x = 0.5,
    xend = ifelse(typeFactor == "Accession", 1, 2)
  ) -> medianDF

breaks <- c(0:3, medianDF$medianDegree)

degreePlot <- ggplot() + theme_bw(base_size = 11) + 
  geom_violin(data = degreeDF, aes(x = as.numeric(typeFactor), y = degreeLog, col = typeFactor, fill = typeFactor), alpha = 0.5) + 
  geom_segment(data = medianDF, aes(x = x, xend = xend, y = medianDegree, yend = medianDegree, col = typeFactor), linetype = "dashed") +
  geom_boxplot(data = degreeDF, aes(x = as.numeric(typeFactor), y = degreeLog, col = typeFactor), alpha = 0.5, width = 1/3) + 
  scale_color_manual(values = c(accessionColor, proteoformColor)) + 
  scale_fill_manual(values = c(accessionColor, proteoformColor)) + 
  scale_x_continuous(name = "", breaks = c(1, 2), labels = c("Gene", "Proteoform"), expand = c(0, 0), limits = c(0.5, 2.5)) +
  scale_y_continuous(name = "Degree", breaks = breaks, labels = c(10^(breaks)), expand = c(0, 0), limits = c(0, 1.05 * max(degreeDF$degreeLog))) + 
  theme(legend.position = "none",
        axis.text.y = element_text(color = c(rep("black", 4), accessionColor, proteoformColor)),
        axis.ticks.y = element_line(color = c(rep("black", 4), accessionColor, proteoformColor)),
        panel.border = element_rect(color = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = c(rep("grey90", 4), NA, NA)),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(color = "black", size = 0.25))


png("docs/fig_5C.png", height = 12, width = 8, units = "cm", res = 600)
degreePlot
dummy <- dev.off()

ratioDensity <- density(matchingDF$ratio)
matchingDF %>% 
  mutate(
    dx = map_dbl(ratio, getDodge)
  ) %>%
  arrange(ratio) -> matchingDF

medianRatio <- median(matchingDF$ratio)
breaks <- c(-2:1, round(medianRatio, digits = 2))
breakLabels <- c(10^(-2:1), round(10^medianRatio, digits = 2))

ratioPlot <- ggplot(data = matchingDF) + theme_bw(base_size = 11) + 
  geom_point(aes(x = 1 - 0.6 * dx, y = ratio, col = ratio), size = 0.5, alpha = 0.8) +
  geom_segment(aes(x = 0.5, xend = 1, y = medianRatio, yend = medianRatio), linetype = "dashed") +
  geom_boxplot(aes(x = 1, y = ratio), col = "black", alpha = 0.5, width = 1/3, outlier.shape = NA) +
  scale_color_scico(begin = 0, end = (1 + -max(matchingDF$ratio) / min(matchingDF$ratio))/2, palette = "cork") +
  scale_x_continuous(name = "", breaks = c(1), labels = "Ratio", expand = c(0, 0), limits = c(0.5, 1.5)) +
  scale_y_continuous(name = "Ratio", breaks = breaks, labels = breakLabels) + 
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = c(rep("grey90", 4), NA)),
    panel.border = element_rect(color = "white"),
    axis.line = element_line(color = "black", size = 0.25)
  )

# axis.text.y = element_text(color = c(rep("black", 2), NA, rep("black", 2))),

png("docs/fig_5D.png", height = 12, width = 4, units = "cm", res = 600)
ratioPlot
dummy <- dev.off()


## Plot degree p

degreeAccessions <- getDegreDF(matchingDF$degreeAccession)
degreeProteoforms <- getDegreDF(matchingDF$degreeProteoform)

degreeAll <- c(degreeAccessions$degree, degreeProteoforms$degree)
pAll <- c(degreeAccessions$p, degreeProteoforms$p)
matching <- c(rep("Gene", nrow(degreeAccessions)), rep("Proteoform", nrow(degreeProteoforms)))

plotDF <- data.frame(degree = degreeAll, p = pAll, matching, stringsAsFactors = F)
plotDF$degree <- log10(plotDF$degree)
plotDF$p <- log10(plotDF$p)
plotDF$matching <- factor(plotDF$matching, levels = c("Gene", "Proteoform"))

plot <- ggplot() + theme_bw(base_size = 11)
plot <- plot + geom_point(data = plotDF, aes(x = degree, y = p, shape = matching, col = matching), alpha = 0.8, size = 2)

plot <- plot + scale_color_manual(name = "Matching", values = c(accessionColor, proteoformColor))

plot <- plot + xlab("Degree [log10]")
plot <- plot + ylab("p [log10]")

plot <- plot + theme(legend.position = "none")


## Plot degree comparison

matchingDF %>% 
  mutate(
    degreeAccessionLog = log10(degreeAccession),
    degreeProteoformLog = log10(degreeProteoform)
  ) -> matchingDF

maxDegree <- max(matchingDF$degreeAccessionLog, matchingDF$degreeProteoformLog)

plot <- ggplot() + theme_bw(base_size = 11) + 
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", size = 0.5) + 
  geom_point(data = matchingDF, aes(x = degreeAccessionLog, y = degreeProteoformLog, col = ratio), alpha = 0.8) +
  scale_color_scico(begin = 0, end = (1 + -max(matchingDF$ratio) / min(matchingDF$ratio))/2, palette = "cork") +
  scale_x_continuous(name = "Degree Gene", limits = c(0, maxDegree), breaks = 0:3, labels = 10^(0:3)) + 
  scale_y_continuous(name = "Degree Proteoform", limits = c(0, maxDegree), breaks = 0:3, labels = 10^(0:3)) + 
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "white"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.25)
  )

png("docs/fig_5E.png", height = 12, width = 12, units = "cm", res = 600)
plot
dummy <- dev.off()


