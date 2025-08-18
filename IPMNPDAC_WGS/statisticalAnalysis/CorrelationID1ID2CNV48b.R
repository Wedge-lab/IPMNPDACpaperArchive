
library("ggpubr")

datapath = '/mnt/b/1_3July23VersionDataPlotCode/SupplementaryInformation/FigureS2/'
id83cn48b = read.csv(paste0(datapath,'41id83cn48bforcorr.csv'))

p1 <- ggscatter(id83cn48b, x = "ID1", y = "CNV48B", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, 
          #cor.method = "pearson",
          cor.method = 'spearman',
          cor.coef.size = 6,  
          xlab = "ID1 activity", ylab = "CN48B activity",
          color = "blue", 
          size = 2, # Points color, shape and size
          ylim=c(0,161),
          font.x = c(14, "plain"),
          font.y = c(14, "plain"),
          font.xtickslab = 14,
          font.ytickslab = 14,
          add.params = list(color = "red", fill = "lightgray")) +
          theme(plot.margin = margin(10, 10, 10, 50)) # add space to the right

p2 <- ggscatter(id83cn48b, x = "ID2", y = "CNV48B", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, 
          #cor.method = "pearson",
          cor.method = 'spearman',
          cor.coef.size = 6, 
          xlab = "ID2 activity", ylab = "CN48B activity",
          color = "blue", 
          size = 2, # Points color, shape and size
          ylim=c(0,161),
          font.x = c(14, "plain"),
          font.y = c(14, "plain"),
          font.xtickslab = 14,
          font.ytickslab = 14,
          add.params = list(color = "red", fill = "lightgray"))+
          theme(plot.margin = margin(10, 10, 10, 50))  # add space to the left
blank <- ggplot() + theme_void()


#combined <- ggarrange(p1, blank, p2,
                      #ncol = 3,
                     # widths = c(2, 0.2, 2),
                      #align = "h")

options(repr.plot.width = 14, repr.plot.height = 5)

#ggarrange(p1, blank, p2,ncol = 3, nrow = 1)
ggarrange(p1, p2,ncol = 2, nrow = 1)
#ggexport(combined, filename = "combined_plot_fixed.pdf", width = 8, height = 6)
