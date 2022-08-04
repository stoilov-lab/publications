library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)


setwd('~/Projects/Msi1/Paper figures/Figures_inkescape/Draft_1/Code/WesterbBlot-and-PCR/')
wb<-read.csv(file = 'Gnat_WB.csv', head=T)
wb

wb<-wb %>% group_by(Group, Label) %>% mutate(scaledIntensity=Value/median(Value)) %>% 
  select(-Value) %>% pivot_wider(values_from=scaledIntensity, names_from="Label") %>%
  mutate(ratio=HA/T7, Sample="Protein") %>% ungroup() %>% select(-c(HA, T7, Group))

meanData<-wb %>%  group_by(Name, Sample) %>% summarise(`Gnat1 Ratio (HA/T7)`=mean(ratio), stdErr=sd(ratio)/sqrt(n()-1))

pcr<-read.csv(file="GNAT1_GFP_PCBP2_MSI1_2022-07-22_Walter.csv", header = T,quote="\'", skip = 19) %>% 
  filter(Content=="Unkn", Sample!="Untransfected") %>% mutate(Cq=as.numeric(Cq))

pcr<-pcr %>% select(-c(Content, Fluor)) %>% 
  pivot_wider(values_from = "Cq", names_from = "Target") %>%
  mutate(ratio=2^(T7-HA), 
         Name=ifelse(Sample=="GFP",Sample,str_to_title(Sample)), Sample="mRNA") %>% 
  select(-c(HA, T7))



meanData<-pcr %>% group_by(Name, Sample) %>% 
  summarise(`Gnat1 Ratio (HA/T7)`=mean(ratio), stdErr=sd(ratio)/sqrt(n()-1)) %>%
  rbind(.,meanData)


meanData<-meanData %>% mutate(Name=str_wrap(Name, width=10)) %>% 
  mutate(Name=factor(Name, levels=c("GFP","Pcbp2","Msi1")), 
         Sample=factor(Sample,levels=c("Protein","mRNA")))

allData<-bind_rows(wb,pcr) %>% mutate(Sample=factor(Sample,levels=c("Protein","mRNA")))

meanData %>%
  ggplot(aes(x=Name, y=`Gnat1 Ratio (HA/T7)`, fill=Sample)) +
  geom_col(position=position_dodge2(padding=0.05), color="black", size=1)+
  geom_errorbar(aes(ymin=`Gnat1 Ratio (HA/T7)`-stdErr, ymax=`Gnat1 Ratio (HA/T7)`+stdErr), 
                position=position_dodge2(padding=0.55), size=1.5)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.05)) +
  scale_fill_grey(start=0.5, end=0.9)+
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 1),
        axis.ticks = element_line(color = 'black', size = 1),
        axis.text = element_text(family="Arial",face="bold",size=14, color="black"),
        axis.title.y = element_text(family="Arial",face="bold",size=16, color="black"),
        axis.title.x=element_blank(),
        legend.text = element_text(family="Arial",face="bold",size=15, color="black"),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.95), 
        legend.box = "horizontal")


meanData %>% 
  ggplot(aes(x=Name, y=`Gnat1 Ratio (HA/T7)`)) +
  geom_col(position=position_dodge2(padding=0.05), fill="lightgrey",color="black", size=1)+
  geom_errorbar(aes(ymin=`Gnat1 Ratio (HA/T7)`-stdErr, ymax=`Gnat1 Ratio (HA/T7)`+stdErr), 
                position=position_dodge2(padding=0.55), size=1.5,
                width=0.5)+
  geom_point(data=allData, 
             aes(x=Name, y=ratio), 
             position=position_dodge2(width=0.7),
             size=3)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.25)) +
  theme_classic() +
  facet_grid(~Sample) +
  theme(axis.line = element_line(color = 'black', size = 1),
        axis.ticks = element_line(color = 'black', size = 1),
        axis.text = element_text(family="Arial",face="bold",size=14, color="black"),
        axis.title.y = element_text(family="Arial",face="bold",size=16, color="black"),
        axis.title.x=element_blank(),
        legend.text = element_text(family="Arial",face="bold",size=15, color="black"),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.95), 
        legend.box = "horizontal",
        strip.text.x= element_text(family="Arial",face="bold",size=15, color="black"))



pcrAov<-pcr %>% aov(ratio ~ Name, data=.)
summary(pcrAov)
TukeyHSD(pcrAov)

wbAov<-wb %>% aov(ratio ~ Name, data=.)
summary(wbAov)
TukeyHSD(wbAov)



