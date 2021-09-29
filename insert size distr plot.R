IS <- read.table('IS.txt', header = F)
names(IS) <- c('insert_size', 'count')
library(ggplot2)
library(tidyr)
IS %>% ggplot(aes(x= insert_size, y= count)) +
  geom_col() +
  xlim(0,2000) +
  scale_y_continuous(labels = scales::scientific) +
  geom_vline(xintercept = c(90,559)) +
  geom_rect(aes(xmin = 559, xmax = Inf, ymin = 0, ymax = Inf, fill = 'deletion'), color =NA, alpha = 0.05) +
  geom_rect(aes(xmin = -Inf, xmax = 90, ymin = 0, ymax = Inf, fill = 'insertion'), color =NA, alpha = 0.05) +
  labs(title = "Insert size distribution", x = 'insert size (bp)') +
  theme_bw() +
  theme(plot.title = element_text(color = 'red', size = 12, hjust = 0.5, face = "bold" )) +
  scale_fill_manual("Anomalous distance", values = c('seagreen2', 'blue'), guide = guide_legend(override.aes = list(alpha = 1)))

typeof(IS)
