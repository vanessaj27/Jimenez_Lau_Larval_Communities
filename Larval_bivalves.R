####### TN401 Bivalve Data Figures #######


library(vegan)
library(ggplot2)
library(ggforce)
library(mvabund)
library(ape)
library(pairwiseAdonis)
library(devtools)
library(dplyr)
library(GGally)
library(tidyverse)
library(ggsci)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# New # 09/17/2025
setwd("/Users/Nessa/Desktop")
data <- read.csv('TN401_Larval_Sample_Data_vj_0925.csv')
str(data)

# Remove unnecessary columns for analysis 
data <- data[,-c(13:19, 23:30)]

# Remove Mata Tolu samples
data <- data[data$RecoveryDive != 'J21416',]

# Remove Groups not included in analysis
data <- subset(data, Genus.Group != 'Egg')
data <- subset(data, Genus.Group != 'Holoplankton')
data <- subset(data, Genus.Group != 'Shell')
data <- subset(data, Genus.Group != 'Pteropod')

# Make Sites factors
data$Site <- factor(data$Site, levels = c('Kilo Moana', 'Tow Cam', 'Tahi Moana', 'ABE', 'Tui Malila', 'Mariner'))
levels(data$Site)

# Make Genus a factor
data$Genus.Group <- as.factor(data$Genus.Group)
levels(data$Genus.Group)

data$Larval.Type <- as.factor(data$Larval.Type)
levels(data$Larval.Type)

# Make Meters above bottom a factor
data$Meters.Above.Bottom <- as.factor (data$Meters.Above.Bottom)
levels(data$Meters.Above.Bottom)

# Subset Bivalve Data 
bivalve_data <- subset(data, Larval.Type == 'Bivlave' )

# Subset Mclane/Syprid data from Larval Trap Data 
bivalve_mcsy <- subset(bivalve_data, Meters.Above.Bottom != 'na')
bivalve_lt <- subset(bivalve_data, Meters.Above.Bottom == 'na')


# Mclane/Syprid Bivalve Stacked Plots 
bivalve_mcsy$larvae.30m3 <- as.numeric (bivalve_mcsy$larvae.30m3)

ggplot(bivalve_mcsy, aes(x = Site, y = larvae.30m3, fill = Genus.Group)) +
  geom_bar(position='stack', width = 0.75, stat = 'identity') +
  facet_grid(Meters.Above.Bottom~., scales = 'free') +
  labs(x = '',
       y = expression(' Larvae per 30' ~m^3)) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45,  hjust=1, color = 'black'),
        legend.text = element_text(face = 'italic'),
        legend.title = element_blank()) +
  scale_x_discrete(name = "",
                   limits = c("Kilo Moana", 'Tow Cam', 'Tahi Moana', 'ABE', 'Tui Malila', 'Mariner')) +
  theme_bw() 

pal <- c('#00cec9','#fab1a0', '#74b9ff','#ff7675', '#a29bfe','#ffeaa7', '#00b894','#e17055', '#6c5ce7', "#B3DE69",   '#fdcb6e' )

m_8mab <- subset(bivalve_mcsy, Meters.Above.Bottom == '8-15 mab')
m_8 <- ggplot(m_8mab, aes(x = Site, y = larvae.30m3, fill = Genus.Group)) +
  geom_bar(position='stack', stat = 'identity') +
  facet_grid(Meters.Above.Bottom~., scales = 'free') +
  labs(x = '',
       y = '',
       fill = 'Bivalve Taxa') +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text( color = 'black'),
        legend.text = element_text(face = 'italic')) +
  scale_fill_manual(limits = c("Bathymodiolus"  ,  "Galeommatoide"    ,   "Lucinidae"  , "Myida" ,  "Mytilidae"     ,   "Pinnidae"  , "Vesicomyid" ,  "Other"  ), values = pal) 

  #scale_x_discrete(name = "", limits = c("KM", 'TC', 'TM', 'ABE', 'TuM', 'Mar'))
m_8

m_2mab <- subset(bivalve_mcsy, Meters.Above.Bottom == '2 mab')
m_2 <- ggplot(m_2mab, aes(x = Site, y = larvae.30m3, fill = Genus.Group)) +
  geom_bar(position='stack', stat = 'identity') +
  facet_grid(Meters.Above.Bottom~., scales = 'free') +
  labs(x = '',
       y = '',
       fill = 'Bivalve Taxa') +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(color = 'black'),
        legend.text = element_text(face = 'italic')) +
  theme_classic() +
  scale_fill_manual(limits = c("Bathymodiolus"  ,  "Galeommatoide"    ,   "Lucinidae"  , "Myida" ,  "Mytilidae"     ,   "Pinnidae"  , "Vesicomyid" ,  "Other"   ), values = pal) +
  scale_x_discrete(name = "",
                   limits = c("Kilo Moana", 'Tow Cam', 'Tahi Moana', 'ABE', 'Tui Malila', 'Mariner'))

m_2
#scale_x_discrete(name = "",
                   #limits = c("KM", 'TC', 'TM', 'ABE', 'TuM', 'Mar'))

m_45mab <- subset(bivalve_mcsy, Meters.Above.Bottom == '26-45 mab')
m_45 <- ggplot(m_45mab, aes(x = Site, y = larvae.30m3, fill = Genus.Group)) +
  geom_bar(position='stack', stat = 'identity') +
  facet_grid(Meters.Above.Bottom~., scales = 'free') +
  labs(x = '',
       y = '',
       fill = 'Bivalve Taxa') +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(color = 'black'),
        legend.text = element_text(face = 'italic')) +
  scale_fill_manual(limits = c("Bathymodiolus"  ,  "Galeommatoide"    ,   "Lucinidae"  , "Myida" ,  "Mytilidae"     ,   "Pinnidae"  , "Vesicomyid" ,  "Other"   ), values = pal) +
  scale_x_discrete(name = "",
                   limits = c("Kilo Moana", 'Tow Cam', 'Tahi Moana', 'ABE', 'Tui Malila', 'Mariner'))
m_45

library (ggpubr)
fig <- ggarrange(m_45 + theme(axis.text.x = element_blank(),
                              axis.title.y  = element_blank()) +
                   facet_wrap(.~'26-45 mab'),
                 m_8 +
                   theme(axis.text.x = element_blank(),
                         axis.title.y  = element_blank()) +
                   facet_wrap(.~'8-15 mab'), 
                 m_2 + theme(axis.title.y  = element_blank()) +
                   facet_wrap(.~'2 mab'),
                 ncol = 1, nrow = 3, common.legend = TRUE, legend = 'none')

fig <- annotate_figure(fig,
                       left = text_grob(expression('Larvae per 30' ~m^3), rot = 90))


########## Bivalve Relative Abundance Figures ###########
#bmol.w <- b.larva[,-c(1,3:11)]
bmol.w <- bivalve_mcsy %>%
  aggregate(larvae.30m3 ~ Site + Meters.Above.Bottom + Genus.Group, sum)
bmol.w <- reshape(data = bmol.w,
                  timevar = 'Genus.Group',
                  idvar = c('Site', 'Meters.Above.Bottom'),
                  direction = 'wide')
bmol.w[is.na(bmol.w)] <- 0
site_totals <- rowSums(bmol.w[,3:10])
rel_a <- bmol.w[,3:10] / site_totals
rel_a$Site <- bmol.w$Site
rel_a$Meters.Above.Bottom <- bmol.w$Meters.Above.Bottom

bmol.l <- gather(rel_a, Genus.Group, rel_a, larvae.30m3.Bathymodiolus:larvae.30m3.Vesicomyid, factor_key= TRUE)
#levels(bmol.l$ID) <- c("Bathymodiolus"  ,  "Galeommatoide"    ,   "Lucinidae"  , "Myida" ,  "Mytilidae"     ,   "Pinnidae"  , "Vesicomyid", )
mr_8mab <- subset(bmol.l, Meters.Above.Bottom == '8-15 mab')
mr_8 <- ggplot(mr_8mab, aes(x = Site, y = rel_a, fill = Genus.Group)) +
  geom_bar(position='stack', stat = 'identity') +
  facet_grid(Meters.Above.Bottom~., scales = 'free') +
  labs(x = '',
       y = '',
       fill = 'Bivalve Taxa') +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text( color = 'black'),
        legend.text = element_text(face = 'italic')) +
  scale_fill_manual(limits = c("Bathymodiolus"  ,  "Galeommatoide"    ,   "Lucinidae"  , "Myida" ,  "Mytilidae"     ,   "Pinnidae"  , "Vesicomyid", 'Other' ), values = pal) +
  scale_x_discrete(name = "",
                   limits = c("KM", 'TC', 'TM', 'ABE', 'TuM', 'Mar')) 


mr_2mab <- subset(bmol.l, Meters.Above.Bottom == '2 mab')
mr_2 <- ggplot(mr_2mab, aes(x = Site, y = rel_a, fill = Genus.Group)) +
  geom_bar(position='stack', stat = 'identity') +
  facet_grid(Meters.Above.Bottom~., scales = 'free') +
  labs(x = '',
       y = '',
       fill = 'Bivalve Taxa') +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(color = 'black'),
        legend.text = element_text(face = 'italic')) +
  theme_classic() +
  scale_fill_manual(limits = c("Bathymodiolus"  ,  "Galeommatoide"    ,   "Lucinidae"  , "Myida" ,  "Mytilidae"     ,   "Pinnidae"  , "Vesicomyid"), values = pal) +
  scale_x_discrete(name = "",
                   limits = c("KM", 'TC', 'TM', 'ABE', 'TuM', 'Mar'))

mr_45mab <- subset(bmol.l, Meters.Above.Bottom == '26-45 mab')
mr_45 <- ggplot(mr_45mab, aes(x = Site, y = rel_a, fill = Genus.Group)) +
  geom_bar(position='stack', stat = 'identity') +
  facet_grid(Meters.Above.Bottom~., scales = 'free') +
  labs(x = '',
       y = '',
       fill = 'Bivalve Taxa') +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(color = 'black'),
        legend.text = element_text(face = 'italic')) +
# Fix limits
  scale_fill_manual(limits = c("Bathymodiolus"  ,  "Galeommatoide"    ,   "Lucinidae"  , "Myida" ,  "Mytilidae"     ,   "Pinnidae"  , "Vesicomyid" ), values = pal)+
  scale_x_discrete(name = "",
                   limits = c("KM", 'TC', 'TM', 'ABE', 'TuM', 'Mar'))

fig2 <- ggarrange(mr_45 + theme(axis.text.x = element_blank(),
                                axis.title.y  = element_blank()) +
                    facet_wrap(.~'26-45 mab'),
                  mr_8 +
                    theme(axis.text.x = element_blank(),
                          axis.title.y  = element_blank()) +
                    facet_wrap(.~'8-15 mab'), 
                  mr_2 + theme(axis.title.y  = element_blank()) +
                    facet_wrap(.~'2 mab'),
                  ncol = 1, nrow = 3, common.legend = TRUE, legend = 'none')

fig2 <- annotate_figure(fig2,
                        left = text_grob(expression('Relative Abundance'), rot = 90))

######### Final Plot Plate with Bivalve Abundance/30m3 and Relativ Abundances ########
legend <- get_legend(mr_45)
legends <- ggarrange(legend, nrow=1, ncol= 1)

ggarrange(fig, NULL, fig2,legends,  ncol = 4, nrow =1, widths=c(1,0.09,1,0.7), common.legend = TRUE, labels = c('A', '', 'B', ""))



########## Bivalve Larval Trap Figures ###########

#### Figure out ind.day calculation again!!! 

bivalve_lt<- read.csv('Larval_Traps.csv')

bivalve_lt$Site <- factor(bivalve_lt$Site, levels = c('Kilo Moana', 'Tow Cam', 'Tahi Moana', 'ABE', 'Tui Malila', 'Mariner'))
levels(bivalve_lt$Site) <- c('KM', 'TC', 'TM', 'ABE', 'TuM', 'Mar')


#blt$ID <- factor(blt$Identity, levels = c( "Bathymodiolus"  ,  "Galeommatoide"    ,   "Lucinidae"  , "Myida" ,  "Mytilidae"     ,   "Pinnidae"  , "Vesicomyid" ,  "Unknown"   ))

levels(bivalve_lt$Site)

blt.a <- ggplot(blt, aes(x = Site2, y = ind.day, fill = ID)) +
  geom_bar(position='stack', stat = 'identity') +
  labs(x = '',
       y = expression ("Larvae / day",
                       fill = 'Bivalve Taxa')) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text( color = 'black'),
        legend.text = element_text(face = 'italic')) +
  scale_fill_manual(limits = c( "Bathymodiolus"  ,  "Galeommatoide"    ,   "Lucinidae"  , "Myida" ,  "Mytilidae"   ,   "Pinnidae"  , "Vesicomyid" ,  "Unknown"  ), values = pal) +
  scale_x_discrete(limits = c('KM', 'TC', 'TM', 'ABE', 'TuM', 'Mar'))

gglt <- ggplot(glt, aes(x = Site2, y = ind.day, fill = ID)) +
  geom_bar(position='stack', stat = 'identity', show.legend = FALSE) +
  labs(x = '',
       y = expression ( 'Larval per Day',
                        fill = 'Gastropod Taxa')) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text( color = 'black'),
        legend.text = element_text(face = 'italic')) +
  scale_fill_manual(values = lt_colors)+
  scale_x_discrete(limits = c('KM', 'TC', 'TM', 'ABE', 'TuM', 'Mar'))


########### END FOR NOW ########## 9/17/25








