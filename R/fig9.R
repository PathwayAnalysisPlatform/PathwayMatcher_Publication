# 
# Plot percentage of proteoforms matched with the phosphoproteoforms created from an experimental dataset.
# (Ochoa D, Jarnuczak AF, Gehre M, Soucheray M, Kleefeldt AA, Vieitez C, et al. The functional landscape of the human phosphoproteome. bioRxiv. 2019:541656. doi:10.1101/541656, Supp table 2)
# Each proteoform is matched against proteoforms in the database that share the same accession.
# From all possible proteoforms, we calculate the percentage of them that got matched using each matching.
#
start_time <- Sys.time()

# Libraries

library(ggplot2)
library(cowplot)
library(plyr)
source("R/plotPercentages.R")

# Parameters
percentagesFile <- "resources/sensitivity/percentagesFilePhosphoproteoforms.csv"

# Main script

print(paste("Loading data from", percentagesFile))
percentages <- read.csv(percentagesFile, header = T)

print(paste("Plotting percentages"))
percentagesPlot <- plotPercentages(percentages)
percentagesPlot

png("docs/fig11.png", height = 12, width = 12, units = "cm", res = 600)
plot(percentagesPlot)
dummy <- dev.off()

print("Finished")
print(Sys.time() - start_time)
