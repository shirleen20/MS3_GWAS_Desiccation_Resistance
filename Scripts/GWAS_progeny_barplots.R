#____________________________________26/11/2025___________________________________________


# We create Barplots for GWAS progenies but without considering mortality at 8hours

library("tidyverse")
library("ggplot2")
library("addinslist")
library("ggThemeAssist")
library("Rmisc")
library("psych")
library("ggpubr")
library("gridExtra")
library("cowplot")


rm(list=ls())

#________________Import desiccation resistance data for progenies used for GWAS ___________________

DF <- read_csv("Data/GWAS_progenies.csv")   

#_______First calculate the median time to mortality and se values for parental lines __________

#. Expand data

DF01 <- DF %>% 
  group_by(Bioassay, Generation, Line, Replicate, Hours) %>% 
  dplyr::filter(Hours != "8") %>% 
  slice(rep(row_number(), times = Mortality)) %>% 
  ungroup() 

DF_summary <- DF01 %>% 
  dplyr::group_by(Line) %>% 
  dplyr::summarise(Mean=mean(Hours),Median=median(Hours),n=length(Hours),sd=sd(Hours),se=sd(Hours)/sqrt(length(Hours))) %>% 
  ungroup()  


#______________________________

# Make bar plot for QTL1F1



# Make bar plot for QTL1F4

DF_QTLF4 <- DF %>% dplyr::filter(Line =="QTL1F4") 


Red <- c("16", "20", "47")
Gray100 <- c("8", "20", "23","26","29","32","35","38","41","44","50","53","56","59","62")

DF_QTLF4_df <- DF_QTLF4 %>% 
  dplyr::mutate(Zone = ifelse(Hours %in% Red, "Red", "Gray100"))


p1 <- DF_QTLF4_df %>% 
  ggplot() +
  #stat_summary(aes(x = Hours, y = Mortality, fill = Zone),fun.y = sum, geom = "col", width=2, colour = "black")+
  stat_summary(aes(x = Hours, y = Mortality, fill = Zone),  #color by Zone
    fun = function(x) {
      s <- sum(x, na.rm = TRUE)
      if (s == 0) NA_real_ else s },# remove flat bars
    geom = "col",width = 2,colour = "black") +
  ggtitle("F4 (n = 594, median = 32, se = 0.3)")+ 
  theme_bw()+
  scale_fill_manual(values = c(Red = "red", Gray100 ="gray100")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9),axis.text.x = element_blank(),axis.title = element_blank())+ #Minimal theme + blue fill color
  ylim(0,120) + scale_x_continuous(breaks = seq(0,70, by=10)) +
  theme(panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) 

#______________________________

# Make bar plot for QTL1F5

DF_QTLF5 <- DF %>% dplyr::filter (Line =="QTL1F5") 

# add another column to DF_QTLF4 that has the following colour categories

Red <- c("16","20","48.5")
Gray100 <- c("8","20","23","26","29","32","35","38","41","44")

DF_QTLF5_df <- DF_QTLF5 %>% 
  dplyr::mutate(Zone = ifelse(Hours %in% Red, "Red", "Gray100"))

p2 <- DF_QTLF5_df %>% 
  ggplot() +
  #stat_summary(aes(x = Hours, y = Mortality, fill = Zone),fun.y = sum, geom = "col", width=2, colour = "black")+
  stat_summary(aes(x = Hours, y = Mortality, fill = Zone),  #color by Zone
               fun = function(x) {
                 s <- sum(x, na.rm = TRUE)
                 if (s == 0) NA_real_ else s },# remove flat bars
               geom = "col",width = 2,colour = "black") +
  ggtitle("F5 (n = 594, median = 35, se = 0.3)")+ 
  theme_bw()+
  scale_fill_manual(values = c(Red = "red", Gray100 ="gray100")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9),axis.text.x = element_blank(),axis.title = element_blank())+ #Minimal theme + blue fill color
  ylim(0,80)+
  scale_x_continuous(breaks = seq(0,70,by=10))+
  theme(panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) 


#______________________________

# Make bar plot for QTL1F6

DF_QTLF6 <- DF %>% dplyr::filter (Line =="QTL1F6") 

# add another column to DF_QTLF4 that has the following colour categories

Red <- c("16","20","47.5")
Gray100 <- c("8","20","23","26","29","32","35","38","41","44","47")

DF_QTLF6_df <- DF_QTLF6 %>% 
  dplyr::mutate(Zone = ifelse(Hours %in% Red, "Red", "Gray100"))

p3 <- DF_QTLF6_df %>% 
  ggplot() +
  #stat_summary(aes(x = Hours, y = Mortality, fill = Zone),fun.y = sum, geom = "col", width=2, colour = "black")+
  stat_summary(aes(x = Hours, y = Mortality, fill = Zone),  #color by Zone
               fun = function(x) {
                 s <- sum(x, na.rm = TRUE)
                 if (s == 0) NA_real_ else s },# remove flat bars
               geom = "col",width = 2,colour = "black") +
  ggtitle("F6 (n = 588, median = 35, se =  0.4)")+ #(n=598 with 8 hours included, Median=35)
  theme_bw()+
  scale_fill_manual(values = c(Red = "red", Gray100 ="gray100")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9),axis.title = element_blank())+ #Minimal theme + blue fill color
  ylim(0,80) +
  scale_x_continuous(breaks = seq(0,70,by=10)) +
  theme(panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) 


#___________________Now combine all the plots (p1.1-7.1) from above_______________________________

# Create a figure by combining the different plots
Figure <- ggarrange(p1, p2, p3, ncol=1, nrow =3)

# Annotate the figure by adding a common labels
Figure_annotated <- annotate_figure(Figure,
                                    bottom = text_grob("Hours of desiccation", color = "black",hjust =0.6, face = "bold", size = 10),
                                    left = text_grob("Number of flies", color = "black",rot = 90,face = "bold", size = 10))

ggsave(plot = Figure_annotated, width = 5, height = 4.5, units = "in", dpi = 300,filename = "Figures/Phenotypic distribution in GWAS progenies.jpg")

#_______________________END____________________________
