setwd("/Users/atiteyk2/Documents/matlab_manuscript")
getwd()

############### Compare Signal to the Model

########## EGFR H05
##########
EGFR_data <- t(EGFR_data)
EGFR_H05_data <- t(EGFR_H05_data)
all_EGFR_H05  <- data.frame(cbind(EGFR_data, EGFR_H05_data))

colors <- c("EGFR_data"="blue", "EGFR_H05_data"="red")
p <- ggplot(all_EGFR_H05, aes(x = 1:100)) +
  geom_line(aes(y = EGFR_data, color = "EGFR_data"), size = 1.5) +
  geom_line(aes(y = EGFR_H05_data, color = "EGFR_H05_data"), size = 1.5) +
  labs(x = "Time",
       y = "X",
       color = "Legend") +
  scale_color_manual(values = colors)
p + theme_classic() #+ theme_bw()

########## EGFR H07
##########
EGFR_data <- t(EGFR_data)
EGFR_H07_data <- t(EGFR_H07_data)
all_EGFR_H07  <- data.frame(cbind(EGFR_data, EGFR_H07_data))
std_EGFR_H07 <- data.frame(all_EGFR_H07*0.09)

colors <- c("EGFR_data"="blue", "EGFR_H07_data"="red")
p <- ggplot(all_EGFR_H07, aes(x = 1:100)) +
  geom_line(aes(y = EGFR_data, color = "EGFR_data"), size = 1.5) +
  geom_line(aes(y = EGFR_H07_data, color = "EGFR_H07_data"), size = 1.5) +
  labs(x = "Time",
       y = "X",
       color = "Legend") +
  scale_color_manual(values = colors)
p + theme_classic() #+ theme_bw()

########## EGFR H09
##########
EGFR_data <- t(EGFR_data)
EGFR_H09_data <- t(EGFR_H09_data)
all_EGFR_H09  <- data.frame(cbind(EGFR_data, EGFR_H09_data))
std_EGFR_H09 <- data.frame(all_EGFR_H09*0.15)

colors <- c("EGFR_data"="blue", "EGFR_H09_data"="red")
p <- ggplot(all_EGFR_H09, aes(x = 1:100)) +
  geom_line(aes(y = EGFR_data, color = "EGFR_data"), size = 1.5) +
  geom_line(aes(y = EGFR_H09_data, color = "EGFR_H09_data"), size = 1.5) +
  labs(x = "Time",
       y = "X",
       color = "Legend") +
  scale_color_manual(values = colors)
p + theme_classic() #+ theme_bw()

