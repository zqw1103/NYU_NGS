######task3
library(tidyverse)
pop.tbl_df <- read_delim(file = "Demo1pop.info",delim = " ", col_names = FALSE)
pop.tbl_df %>% View()
qopt.tbl_df <- read_delim(file = "Demo1NGSadmix_nowhite.qopt",delim = " ", col_names = F)
qopt.tbl_df %>% View()
qopt.tbl_df <- bind_cols(pop.tbl_df,qopt.tbl_df)
qopt.tbl_df %>% View()
names(qopt.tbl_df) <- c("pop","sample","g1","g2","g3")
qopt.tbl_df %>% View()


qopt.tbl_df.long <- qopt.tbl_df %>% 
  pivot_longer(cols = g1:g3, names_to = 'group', values_to = 'fraction')
qopt.tbl_df.long %>% View()

qopt.tbl_df.long %>%
  ggplot(aes(x=sample,y=fraction,fill=group)) + geom_col(color = "gray", size = 0.1) +
  facet_grid(~ pop, scales = "free", space = "free") +   
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x = element_text(angle = 90))

#####task4
pop.tbl_df <- read_delim(file = "Demo2pop.info",delim = " ", col_names = FALSE)
pop.tbl_df %>% View()
###k3
qopt.k3.tbl_df <- read_delim(file = "Demo2NGSadmix.k3_nowhite.qopt",delim = " ", col_names = F)
qopt.k3.tbl_df %>% View()
###k4
qopt.k4.tbl_df <- read_delim(file = "Demo2NGSadmix.k4_nowhite.qopt",delim = " ", col_names = F)
qopt.k4.tbl_df %>% View()
###k5
qopt.k5.tbl_df <- read_delim(file = "Demo2NGSadmix.k5_nowhite.qopt",delim = " ", col_names = F)
qopt.k5.tbl_df %>% View()
###bind_cols
?bind_cols()
qopt.k3.tbl_df <- bind_cols(pop.tbl_df, qopt.k3.tbl_df)
qopt.k4.tbl_df <- bind_cols(pop.tbl_df, qopt.k4.tbl_df)
qopt.k5.tbl_df <- bind_cols(pop.tbl_df, qopt.k5.tbl_df)

qopt.tbl_df <- bind_cols(pop.tbl_df, qopt.k3.tbl_df, qopt.k4.tbl_df, qopt.k5.tbl_df)
qopt.tbl_df %>% View()

##rename cols
names(qopt.k3.tbl_df) <- c("pop","sample","g1","g2","g3")
qopt.k3.tbl_df %>% View()
names(qopt.k4.tbl_df) <- c("pop","sample","g1","g2","g3","g4")
qopt.k4.tbl_df %>% View()
names(qopt.k5.tbl_df) <- c("pop","sample","g1","g2","g3","g4","g5")
qopt.k5.tbl_df %>% View()

###
names(qopt.tbl_df) <- c("pop","sample","3_g1","3_g2","3_g3","4_g1","4_g2","4_g3","4_g4","5_g1","5_g2","5_g3","5_g4","5_g5" )
qopt.tbl_df %>% View()

###pivot_longer
qopt.k3.tbl_df.long <- qopt.k3.tbl_df %>% 
  pivot_longer(cols = -c(pop, sample), names_to = 'group', values_to = 'fraction')
qopt.k3.tbl_df.long %>% View()

qopt.k4.tbl_df.long <- qopt.k4.tbl_df %>%
  pivot_longer(cols = -c(pop, sample), names_to = 'group', values_to = 'fraction')
qopt.k4.tbl_df.long %>% View()

qopt.k5.tbl_df.long <- qopt.k5.tbl_df %>%
  pivot_longer(cols = -c(pop, sample), names_to = 'group', values_to = 'fraction')
qopt.k5.tbl_df.long %>% View()

###
qopt.tbl_df.long <- qopt.tbl_df %>% 
  pivot_longer(cols = -c(pop, sample), names_to = c('k','group'), names_sep = '_', values_to = 'fraction')
qopt.tbl_df.long %>% View()


###plot
qopt.tbl_df.long %>% filter(k == '3') %>%
  ggplot(aes(x=sample,y=fraction,fill=group)) + geom_col(color = "gray", size = 0.1) +
  facet_grid(~ pop, scales = "free", space = "free") +   
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x = element_text(angle = 90))

qopt.tbl_df.long %>% filter(k == '4') %>%
  ggplot(aes(x=sample,y=fraction,fill=group)) + geom_col(color = "gray", size = 0.1) +
  facet_grid(~ pop, scales = "free", space = "free") +   
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x = element_text(angle = 90))

qopt.tbl_df.long %>% filter(k == '5') %>%
  ggplot(aes(x=sample,y=fraction,fill=group)) + geom_col(color = "gray", size = 0.1) +
  facet_grid(~ pop, scales = "free", space = "free") +   
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x = element_text(angle = 90))

##all k in one plot
qopt.tbl_df.long %>% 
  ggplot(aes(x=sample,y=fraction,fill=group)) + geom_col(color = "gray", size = 0.1) +
  facet_grid(k ~ pop, scales = "free", space = "free") +   
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x = element_text(angle = 90))
