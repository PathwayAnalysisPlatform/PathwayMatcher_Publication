## Plot theme

theme_set(theme_bw(base_size = 11))

## Plot functions

plotPercentages <- function(percentages){
    percentages$MatchType <- as.factor(percentages$MatchType)
    percentages$MatchType <- factor(percentages$MatchType, levels=unique(percentages$MatchType[order(percentages$Percentage)]), ordered=TRUE)
    
    p <- ggplot(percentages, aes(x=MatchType, y=Percentage)) +
        geom_point() +
        theme(axis.text = element_text(angle = 45, hjust = 1)) +
        xlab("Proteoform Matching Type") +
        ylab("Percentage of proteoforms matched") + ylim(0.0, 100.0)
    p
}

plotPercentagesSeparated <- function(percentages){
    summarised <- ddply(percentages, .(MatchType), summarize,  Percentage=mean(Percentage))
    percentages$MatchType <- factor(percentages$MatchType, levels = summarised$MatchType[order(-summarised$Percentage)])
    
    p <- ggplot(percentages, aes(x=MatchType, y=Percentage, color=Category))
    if(nrow(percentages) > 14){
        p <- p + geom_boxplot(aes(fill = Category), position = position_dodge(0.0), width = 1.2, size = 0.5)
    }
    p <- p + geom_point() +
        theme(axis.text = element_text(angle = 45, hjust = 1)) +
        xlab("Proteoform Matching Type") +
        ylab("Percentage of proteoforms matched") + ylim(0.0, 100.0)
    p
}