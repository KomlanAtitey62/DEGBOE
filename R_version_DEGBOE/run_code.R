setwd("/Users/atiteyk2/Documents")
getwd()

source("DEGBOE_atitey.R")
library(ggplot2)
##############
##############
load("dataset/LungMutationData.rdata")
mutation.dataset <- data.frame(LungMutationData)
dataset <- mutation.dataset[82,] # 5039 # 82

n = 100
N = 200
H = 0.58
cancer.komlan <- sir_atitey(n,M,H,dataset)

############## 
############## plot
obs.vec <- cancer.komlan$obs
sir.est <- cancer.komlan$sir

sir.data  <- data.frame(cbind(obs.vec, sir.est))

colors <- c("obs.vec"="blue", "sir.est"="red")
p <- ggplot(sir.data , aes(x = 1:100)) +
  geom_line(aes(y = obs.vec, color = "obs.vec"), size = 1) +
  geom_line(aes(y = sir.est, color = "sir.est"), size = 1) +
  labs(x = "Time",
       y = "X",
       color = "Legend") +
  scale_color_manual(values = colors)
p + theme_classic() #+ theme_bw()
