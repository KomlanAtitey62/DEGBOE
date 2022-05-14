setwd("/Users/atiteyk2/Documents/matlab_manuscript")
getwd()

############### convert Output data from the DEGBOE implementation in matlab to .rdata file

library("R.matlab")

######### EGFR
#########
EGFR <- readMat("lung_data/EGFR_Signal.mat")
EGFR_H05 <- readMat("lung_data/EGFR_Model_H05.mat")
EGFR_H07 <- readMat("lung_data/EGFR_Model_H07.mat")
EGFR_H09 <- readMat("EGFR_Model_H09.mat")
err_EGFR_H05 <- readMat("lung_data/err_EGFR_H05.mat")
err_EGFR_H07 <- readMat("lung_data/err_EGFR_H07.mat")
err_EGFR_H09 <- readMat("lung_data/err_EGFR_H09.mat")
#
EGFR_data <- EGFR$obs.data 
EGFR_H05_data <- EGFR_H05$komlan
EGFR_H07_data <- EGFR_H07$komlan
EGFR_H09_data <- EGFR_H09$komlan
err_EGFR_H05_data <- err_EGFR_H05$err
err_EGFR_H07_data <- err_EGFR_H07$err
err_EGFR_H09_data <- err_EGFR_H09$err

######### KRAS
#########
KRAS <- readMat("lung_data/KRAS_Signal.mat")
KRAS_H05 <- readMat("lung_data/KRAS_Model_H05.mat")
KRAS_H07 <- readMat("lung_data/KRAS_Model_07.mat")
KRAS_H09 <- readMat("lung_data/KRAS_Model_H09.mat")
err_KRAS_H09 <- readMat("lung_data/err_KRAS_H05.mat")
err_KRAS_H07 <- readMat("lung_data/err_KRAS_H07.mat")
err_KRAS_H05 <- readMat("lung_data/err_KRAS_H09.mat")
#
KRAS_data <- KRAS$obs.data 
KRAS_H05_data <- KRAS_H05$komlan
KRAS_H07_data <- KRAS_H07$komlan
KRAS_H09_data <- KRAS_H09$komlan
err_KRAS_H05_data <- err_KRAS_H05$err
err_KRAS_H07_data <- err_KRAS_H07$err
err_KRAS_H09_data <- err_KRAS_H09$err

######### TP53
#########
TP53 <- readMat("lung_data/TP53_Signal.mat")
TP53_H05 <- readMat("lung_data/TP53_Model_H05.mat")
TP53_H07 <- readMat("lung_data/TP53_Model_H07.mat")
TP53_H09 <- readMat("lung_data/TP53_Model_H09.mat")
err_TP53_H05 <- readMat("lung_data/err_TP53_H05.mat")
err_TP53_H07 <- readMat("lung_data/err_TP53_H07.mat")
err_TP53_H09 <- readMat("lung_data/err_TP53_H09.mat")
#
TP53_data <- TP53$obs.data 
TP53_H05_data <- TP53_H05$komlan
TP53_H07_data <- TP53_H07$komlan
TP53_H09_data <- TP53_H09$komlan
err_TP53_H05_data <- err_TP53_H05$err
err_TP53_H07_data <- err_TP53_H07$err
err_TP53_H09_data <- err_TP53_H09$err













