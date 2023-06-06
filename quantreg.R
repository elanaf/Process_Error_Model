# Quantile regression

library(quantreg)
load("r_values_table.RData")
load("main_dfs.RData")

#For this analysis, we want to find the quantiles of phrag growth by growth rate
## As a note, these calculations were completed with the original r values, which might change

#First, we need to make a table that includes both of these values

dat <- greenhouse %>% 
  filter(Date_Cleaned == "2022-05-16",
         Phrag_Presence == "W") %>% 
  select(Species, Density, Phrag_Presence, Block, Cover.Phrag) %>% 
  left_join(data_read, by = c("Species", "Density", "Phrag_Presence"))
#View(dat)

#Now make the calculations to see the relationship between the two
#Choose the quantiles by changing the tau

rqfit <- rq(Cover.Phrag ~ value, tau = c(.05, .25, .5, .75, .95),data = dat)
summary(rqfit, se=) #need to check about which method of se is best

#Now graph them
#base r graph
color <- c("#ffcccc", "#ff9999", "#ff6666", "#ff3333", "#ff0000")
library(gridExtra)
jpeg("quantile_regression.jpeg", height = 500, width = 500)
plot(Cover.Phrag ~ value, data = dat,
     xlab = "Growth rate (r)",
     ylab = substitute(paste("Proportional ", italic("Phragmites "), "cover"))
)
for (j in 1:ncol(rqfit$coefficients)) {
  abline(coef(rqfit)[, j], col = color[j])
}
legend(x = "topright", legend = c(0.05, 0.25, 0.5, 0.75, 0.95), 
       col = color, lty = 1, title = "Quantiles")
dev.off()

#ggplot graph 
dat %>% 
  ggplot(aes(x = value, y = Cover.Phrag))+
  geom_point() +
  geom_quantile(quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95))
