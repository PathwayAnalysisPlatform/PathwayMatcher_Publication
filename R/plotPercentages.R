## Libraries

library(scico)

## Plot theme

theme_set(theme_bw(base_size = 11))
dotSize <- 2

## Colors

palette <- 'cork'
accessionColor <- scico(n = 1, begin = 0.15, end = 0.15, palette = palette)
proteoformColor <- scico(n = 1, begin = 0.85, end = 0.85, palette = palette)

## Plot functions

plotPercentages <- function(percentages){
    percentages$MatchType <- as.factor(percentages$MatchType)
    percentages$MatchType <- factor(percentages$MatchType, levels=unique(percentages$MatchType[order(percentages$Percentage)]), ordered=TRUE)
    
    p <- ggplot(percentages, aes(x=MatchType, y=Percentage)) +
        geom_point(col = proteoformColor, alpha = 0.5, size = dotSize) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "white"),
            axis.line = element_line(color = "black", size = 0.25)
            ) +
        xlab("Proteoform Matching Type") +
        ylab("Percentage of proteoforms matched") + ylim(0.0, 100.0)
    p
}

plotPercentagesSeparated <- function(percentages){
    summarised <- ddply(percentages, .(MatchType), summarize,  Percentage=mean(Percentage))
    percentages$MatchType <- factor(percentages$MatchType, levels = summarised$MatchType[order(-summarised$Percentage)])
    
    p <- ggplot(percentages, aes(x=MatchType, y=Percentage, color=Category))
    if(nrow(percentages) > 14){
        p <- p + geom_boxplot(position = position_dodge(0.0), width = 1.2, size = 0.5) +
            geom_jitter(width = 0.1, alpha = 0.5, size = dotSize)
    } else {
        p <- p + geom_point(alpha = 0.5, size = 2)
    }
    p <- p + scale_color_manual(values = c(accessionColor, proteoformColor)) +
        scale_fill_manual(values = c(accessionColor, proteoformColor)) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_line(color = c(NA, rep("grey95", 3))),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "white"),
            axis.line = element_line(color = "black", size = 0.25)
              ) +
        xlab("Proteoform Matching Type") +
        ylab("Percentage of proteoforms matched") + ylim(0.0, 100.0)
    p
}