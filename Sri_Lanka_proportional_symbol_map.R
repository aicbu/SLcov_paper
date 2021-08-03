####################################################################################
### 
# R Code for global tree visualization of proportional symbol map of SARS-CoV-2
# sampling sites in Sri Lanka.
# code for reproducible research (prepared for CODECHECK codecheck.org.uk)
#
# This code is part of, and published under the same terms and conditions as, the following publication 
# Chandima Jeewandara, Deshni Jayathilaka, Diyanath Ranasinghe, Nienyun Sharon Hsu, Dinuka Ariyaratne, 
# Tibutius Thanesh Jayadas, Deshan Madushanka, Benjamin B. Lindsey, Laksiri Gomes, Matthew D. Parker, 
# Ananda Wijewickrama, Malika Karunaratne, Graham S. Ogg, Thushan I. de Silva, Gathsaurie Neelika Malavige.
# (2021). Genomic and epidemiological analysis of SARS-CoV-2 viruses in Sri Lanka. 
# Frontiers in Microbiology, section Infectious Diseases. 
#
#
# Code authored by Diyanath Ranasinghe.
# was first written for R(v4.0.1), R Studio(v1.3.1073)
# Last tested on  R(v4.0.1), R Studio(v1.3.1073)
# 
# Recommended first use is to manually execute each Stage incrementally in RStudio.
#
###
####################################################################################
### STAGE 1 - load libraries and import data
####################################################################################

cat("\014") ## clear RStudio console
rm(list=ls()) ## clear the Global Environment

### Loading required packages
library(ggplot2)
library(maps)
library(ggrepel)
library(cowplot)  
library(dplyr)
library(readxl)

### Loading Sri Lankan map from global data 
LK <- map_data("world") %>% filter(region=="Sri Lanka")

### Import B.1.411 sampling counts
dat_sl <- read_xlsx('sl_sampling_data.xlsx', sheet = "B.1.411")

### Import B.1.1.7 sampling counts
dat_uk <- read_xlsx('sl_sampling_data.xlsx', sheet = "B.1.1.7")

####################################################################################
### STAGE 2 - Plotting separate proportional symbol maps for B.1.411 and B.1.1.7
####################################################################################

### Plot B.1.411 sampling counts 
map_1 <- ggplot() + geom_polygon(data = LK, aes(x=long, y = lat, group = group), 
  fill="gray", alpha=0.8)+ coord_quickmap() + 
  geom_point(data=dat_sl, aes(x=long, y=lat, color = Location, size = Tests),
  alpha=.8) + scale_size_continuous(range=c(1,14)) + guides(color=F)+ 
  theme(panel.border = element_rect(colour = "black", fill=F))+
  geom_text_repel(data=dat_sl, aes(x=long, y=lat, label=loc), size=3.2, 
  box.padding = .45, segment.colour = NA, max.overlaps =50)+
  ggtitle("Geographical distribution of B.1.411 sequences in Sri Lanka")+
  theme(legend.text = element_text(size=10),
        legend.direction = "horizontal",
        legend.background = element_rect(color='black', fill = "lightgray"),
        legend.key = element_rect(fill = "lightgray" ),
        panel.border = element_rect(colour = "black", fill=F),legend.position = "na")+
  scale_color_manual(values=c("#cb4bba",
                              "#7cd849",
                              "#cecb46",
                              "#8758d9",
                              "#6870be",
                              "#66da89",
                              "#cc4472",
                              "#c683ba",
                              "#4d8958",
                              "#d74831",
                              "#5dd2c1",
                              "#d98d3a",
                              "#65a9d8",
                              "#c87563",
                              "#c1cd8a"))+
  labs(size = "No of sequences ")+
  annotate(geom = "rect", xmin=79.8, xmax=79.98, ymin=6.8, ymax=6.98, 
           color = "#4d4d4d", fill=NA, linetype = 3)

### Plot B.1.17 sampling counts 
map_2 <- ggplot() +geom_polygon(data = LK, aes(x=long, y = lat, group = group), fill="gray", 
                       alpha=0.7)+ coord_quickmap() + labs(size = "No of sequences ")+
  geom_point(data=dat_uk, aes(x=long, y=lat, color = Location, size = Tests),
             alpha=.8, show.legend = T) + scale_size_continuous(range=c(1,14)) +
  theme(panel.border = element_rect(colour = "black", fill=F))+
  geom_text_repel(data=dat_uk, aes(x=long, y=lat, label=loc), size=3.2, 
                  box.padding = .45, segment.colour = NA, max.overlaps = 50)+
  guides(color=F)+ theme(legend.position = "na")+
  ggtitle("Geographical distribution of B.1.1.7 sequences in Sri Lanka")+
  scale_color_manual(values=c('#8758d9', '#6870be', '#cc4472','#659431', 
                              '#4d8958', '#5dd2c1', '#65a9d8'))+
  annotate(geom = "rect", xmin=79.8, xmax=79.98, ymin=6.8, ymax=6.98, 
             color = "#4d4d4d", fill=NA, linetype = 3)

####################################################################################
### STAGE 3 - Plotting zoomed in panel for suburbs of Colombo
####################################################################################

### Create data layers for B.1.4.11(sl) and B.1.1.7(uk) in suburbs of Colombo

col_data_sl <- dat_sl %>% filter(between(lat,6.8,7)) %>% filter(between(long,79.8,79.98))
col_data_uk <- dat_uk %>% filter(between(lat,6.8,7)) %>% filter(between(long,79.8,79.98))

### Plotting zoomed in map
zoom_map  <- ggplot() +
  geom_polygon(data = LK, aes(x=long, y = lat, group = group), fill="gray", 
               alpha=0.6)+ coord_sf(xlim = c(79.8, 79.98), ylim = c(6.8, 6.98), expand = T)+
  geom_point(data=col_data_sl, aes(x=long, y=lat, colour = 'B.1.411', size= Tests),alpha=0.9)+
  scale_size_continuous(range=c(1,14)) +
  geom_point(data=col_data_uk, aes(x=long, y=lat, colour='B.1.1.7', size= Tests),alpha=0.9)+
  labs(color='Lineage', size='No of Sequences')+
  geom_text_repel(data=col_data_sl, aes(x=long, y=lat, label=name), size=3.5, 
                  box.padding = 0.6, segment.colour = NA) +
  scale_color_manual(values=c('#d9004b', '#60dd93'))+
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        legend.background = element_rect(color='black', fill = "lightgray"),
        legend.key = element_rect(fill = "lightgray", ),
        panel.border = element_rect(colour = "black", fill=F),
        legend.box = "Horizontal", legend.position = "bottom") 

####################################################################################
### STAGE 4 - Combining all plots and saving the output
####################################################################################
        
### Combine maps
map <- plot_grid(map_1,zoom_map, map_2, nrow = 1, ncol = 3 ) 

### Save output
ggsave("SL373_prop_map.pdf",map , width=16.5, height=10, limitsize = T)

####################################################################################
# end
####################################################################################
 
  