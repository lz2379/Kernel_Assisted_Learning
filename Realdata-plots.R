


######################################################
####Figure 1: plot variables
######################################################
#load data
Data = read.csv("data-processed.csv")

par(mfrow = c(2, 4))
par(mar = c(3, 2.5, 0.5, 2))
par(cex = 1.3)
par(mgp = c(1.5, 0.5, 0))
par(family = "serif")

#Plot continuous variables
hist(Data$Age*10, main = "", xlab = "Age")
hist(Data$Weight, xlab = "Weight", main = "")
hist(Data$Height, xlab = "Height", main = "")

#Plot discrete variables
par(mar = c(3, 1, 0.5, 3))
pie(table(Data$Gender), 
    labels = c("Female", "Male"), 
    radius = 1, xlab = "Gender", 
    col  =  gray(seq(0.7,  1,  length  =  2)))
pie(table(Data$Black + 2*Data$Asian), 
    labels = c("Other", "Black", "Asian"), 
    radius = 1, 
    xlab = "          Race", 
    col  =  gray(seq(0.4,  1,  length  =  3)))
pie(table(Data$Enzyme + 2*Data$Amiodarone), 
    labels = c("None", "Enzyme", "Amiodarone", "Both"), 
    radius = 1, 
    xlab = "             Medicine", 
    col  =  gray(seq(0.4,  1,  length  =  3)), 
    init.angle  = 335)
pie(table(Data$VKORC1.AA + 2*Data$VKORC1.AG), 
    labels = c("G/G", "A/A", "A/G"), 
    xlab = "                 VKORC1", 
    radius = 1, 
    col  =  gray(seq(0.4,  1,  length  =  3)), 
    init.angle  = 330)
pie(table(Data$CYP2C9.12 + 2*Data$CYP2C9.13 + 4*Data$CYP2C9.other), 
    labels = c("*1/*1", "*1/*2", "*1/*3", "other"), 
    radius = 1, 
    xlab = "  CYP2C9", 
    col  =  gray(seq(0.4,  1,  length  =  4)), 
    init.angle  = 340)


#######different format
par(mfrow = c(4, 2))
par(mar = c(2, 2.5, 1.5, 2))
par(cex = 1)
par(mgp = c(1.5, 0.5, 0))
par(family = "serif")

#Plot continuous variables
hist(Data$Age*10, main = "Age", xlab = "")
hist(Data$Weight, xlab = "", main = "Weight")
hist(Data$Height, xlab = "", main = "Height")

#Plot discrete variables
par(mar = c(0.5, 1, 2, 3))
pie(table(Data$Gender), 
    labels = c("Female", "Male"), 
    radius = 1, 
    xlab = "", 
    main = "Gender", 
    col  =  gray(seq(0.7,  1,  length  =  2)))
pie(table(Data$Black + 2*Data$Asian), 
    labels = c("Other", "Black", "Asian"), 
    radius = 1, 
    main = "         Race", 
    col  =  gray(seq(0.4,  1,  length  =  3)))
pie(table(Data$Enzyme + 2*Data$Amiodarone), 
    labels = c("None", "Enzyme", "Amiodarone", "Both"), 
    radius = 1, 
    main = "      Medicine", 
    col  =  gray(seq(0.4,  1,  length  =  3)), 
    init.angle  = 335)
pie(table(Data$VKORC1.AA + 2*Data$VKORC1.AG), 
    labels = c("G/G", "A/A", "A/G"), 
    main = "               VKORC1", 
    radius = 1, 
    col  =  gray(seq(0.4,  1,  length  =  3)), 
    init.angle  = 330)
pie(table(Data$CYP2C9.12 + 2*Data$CYP2C9.13 + 4*Data$CYP2C9.other), 
    labels = c("*1/*1", "*1/*2", "*1/*3", "other"), 
    radius = 1, 
    main = "    CYP2C9", 
    col  =  gray(seq(0.4,  1,  length  =  4)), 
    init.angle  = 340)


######################################################
####plot Figure 2: suggested doses
######################################################

par(mfrow = c(1, 1))
par(cex = 1)
par(mar = c(3.5, 2.5, 1, 1))
par(cex = 1.7)
par(mgp = c(1.5, 0.5, 0))
par(family = "serif")
location = "results/realdata/"

#####(Figure 2a)##########  
lol_dose_hat <- unlist(read.csv(paste(location, "lol-real-dose-Y2.csv",  sep = "")))
kol_dose_hat <- unlist(read.csv(paste(location, "kol-real-dose-Y2.csv",  sep = "")))
KAL_dose_hat <- unlist(read.csv(paste(location, "KAL-real-dose-Y2.csv", sep = "")))
discreteQ_dose_hat <- unlist(read.csv(paste(location, "discreteQ-real-dose-Y2.csv", sep = "")))

#ylim_u = max(c(density(lol_dose_hat)$y, density(kol_dose_hat)$y, density(KAL_dose_hat)$y))
graypalette = gray.colors(1, start = 0.5, end = 0.5, gamma = 2.2)

plot(density(lol_dose_hat*(max(Data$Dose)-min(Data$Dose)) + min(Data$Dose)), lwd = 2,  lty = 2, 
     col = graypalette[1], ylim = c(0, 0.20), main = "", 
     xlim = c(min(Data$Dose), max(Data$Dose)))
mtext("Dose", side = 1, line = 2.5, cex = 1.7)
lines(density(kol_dose_hat*(max(Data$Dose)-min(Data$Dose)) + min(Data$Dose)), lwd = 2, col = graypalette[1], lty = 1)
lines(density(KAL_dose_hat*(max(Data$Dose)-min(Data$Dose)) + min(Data$Dose)), lwd = 2, col = 1, lty = 1)
lines(density(Data$Dose), lty = 2, lwd = 2)
legend("topright", legend = c("Original Dose", "KAL", "LOL", "KOL"), lty = c(2, 1, 2, 1), col = c(1, 1, graypalette[1], graypalette[1]))

######(Figure 2b)##########
hist(discreteQ_dose_hat*(max(Data$Dose)-min(Data$Dose)) + min(Data$Dose), breaks = 10, 
     main = "", xlab = "", xlim = c(0, 100))
mtext("Dose", side = 1, line = 2.5, cex = 1.7)

######################################################
####plot figure 3: value fuction
######################################################

KAL_real_value <- unlist(read.csv(paste(location, "KAL-real-value-Y2.csv", sep = "")))
lol_value_hat <- unlist(read.csv(paste(location,  "lol-real-value-Y2.csv", sep = "")))
kol_value_hat <- unlist(read.csv(paste(location,  "kol-real-value-Y2.csv", sep = "")))
discreteQ_value<- unlist(read.csv(paste(location,  "discreteQ-real-value-Y2.csv", sep = "")))

par(mar = c(3.5, 3, 0.5, 0.5))
par(mfrow = c(1, 1))
par(cex = 1.7)

graypalette = gray.colors(1, start = 0.5, end = 0.5, gamma = 2.2)

plot(density(discreteQ_value), col = 1, lty = 2, main = "", ylim = c(0, 38), xlim = c(-0.5, -0.05), lwd = 2)
mtext("Estimated value function", side = 1, line = 2.5,  cex = 1.7)

lines(density(KAL_real_value), col = 1, lty = 1, lwd = 2)
lines(density(lol_value_hat, na.rm = TRUE), col = graypalette[1], lty = 2, lwd = 2)
lines(density(kol_value_hat, na.rm = TRUE), col = graypalette[1], lty = 1, lwd = 2)
legend("topleft", lty = c(1, 2, 1, 2), col = c(1, graypalette[1], graypalette[1], 1), legend = c("KAL", "LOL", "KOL",  "Discretized Q"))
