code library(readr)
library(dplyr)
library(ggplot2)

library(ggsankey)
library(plotly)
library(ggalluvial)


cn_data <- read_tsv("../in_data/tetraploid_passages_cn_changes.tsv")%>% filter(Clone_Name!='Clone 17b')


major_changes_plotting <- cn_data %>% filter(Passage_50_Major %in% c(1,2,3,4,5),Passage_4_Major<6) %>%group_by(Clone_Name,Passage_50_Major,Passage_4_Major)%>% summarise(Total_Width=sum(Width)) %>% group_by(Clone_Name,Passage_50_Major)%>% mutate(Proportion=Total_Width/sum(Total_Width))

major_changes_plotting %>% ggplot(aes(x=as.character(Passage_4_Major),y=Proportion))+geom_col()+facet_grid(rows=vars(Passage_50_Major),cols = vars(Clone_Name))+theme_bw()+labs(x='Major Copy Number at Passage 4')

minor_changes_plotting <- cn_data %>% filter(Passage_50_Minor %in% c(0,1,2)) %>%group_by(Clone_Name,Passage_50_Minor,Passage_4_Minor)%>% summarise(Total_Width=sum(Width)) %>% group_by(Clone_Name,Passage_50_Minor)%>% mutate(Proportion=Total_Width/sum(Total_Width))

minor_changes_plotting %>% ggplot(aes(x=as.character(Passage_4_Minor),y=Proportion))+geom_col()+facet_grid(rows=vars(Passage_50_Minor),cols = vars(Clone_Name))+theme_bw()+labs(x='Minor Copy Number at Passage 4')

major_changes_plotting_across_samples <- major_changes_plotting %>% group_by(Passage_4_Major,Passage_50_Major)%>% summarise(Total_Proportion=sum(Proportion))%>% filter(Passage_50_Major%in%c(3,4),Passage_4_Major>0 & Passage_4_Major<5)%>% mutate(Parsimonious=(Passage_4_Major==Passage_50_Major)|(Passage_4_Major==2&Passage_50_Major==3))


ggplot(data = major_changes_plotting_across_samples,
       aes(axis1 = as.character(Passage_4_Major), axis2=as.character(Passage_50_Major), y = Total_Proportion)) +
  geom_alluvium(aes(fill = Parsimonious)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +

  theme_void()
ggsave('../plots/major_cn_flow.pdf')
