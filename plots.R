


######################################################
####Figure 1: plot variables
######################################################
Data <- read.csv("data-processed.csv")
par(mfrow=c(2,4))
par(mar=c(3,2.5,0.5,2))
par(cex=1.3)
par(mgp=c(1.5,0.5,0))
par(family="serif")
hist(Data$Age*10,main="",xlab="Age")

hist(Data$Wegiht,xlab="Weight",main="")
hist(Data$Height,xlab="Height",main="")
par(mar=c(3,1,0.5,3))
pie(table(Data$Gender),labels=c("Female","Male"),radius=1,xlab="Gender",col = gray(seq(0.7, 1, length = 2)))
pie(table(Data$Black+2*Data$Asian),labels=c("Other","Black","Asian"),radius=1,xlab="          Race",col = gray(seq(0.4, 1, length = 3)))
pie(table(Data$Enzyme+2*Data$Amiodarone),labels=c("None","Enzyme","Amiodarone","Both"),radius=1,xlab="             Medicine",col = gray(seq(0.4, 1, length = 3)),init.angle =335)
pie(table(Data$VKORC1.AA+2*Data$VKORC1.AG),labels=c("G/G","A/A","A/G"),xlab="                 VKORC1",radius=1,col = gray(seq(0.4, 1, length = 3)),init.angle =330)
pie(table(Data$CYP2C9.12+2*Data$CYP2C9.13+4*Data$CYP2C9.other),
    labels=c("*1/*1","*1/*2","*1/*3","other"),radius=1,xlab="  CYP2C9",col = gray(seq(0.4, 1, length = 4)),init.angle =340)


######################################################
####plot Figure 2: suggested doses
######################################################

par(mfrow=c(1,1))
par(cex=1)
par(mar=c(3.5,2.5,1,1))
par(cex=1.7)
#####Figure 2a  
lol_dose_hat <- unlist(read.csv("results/lol-real-dose.csv", sep=""))
kol_dose_hat <- unlist(read.csv("results/kol-real-dose.csv", sep=""))
KAL_dose_hat <- unlist(read.csv("results/KAL-real-dose.csv",sep=""))
discreteQ_dose_hat <- unlist(read.csv("results/discreteQ-real-dose.csv",sep=""))

plot(density(lol_dose_hat*(max(Data$Dose)-min(Data$Dose))+min(Data$Dose)),lwd=2,col=3,ylim=c(0,0.16),main="",xlim=c(min(Data$Dose),max(Data$Dose)))
mtext("Dose",side=1,line=2.5,cex=1.7)
lines(density(kol_dose_hat*(max(Data$Dose)-min(Data$Dose))+min(Data$Dose)),lwd=2,col=4)
lines(density(KAL_dose_hat*(max(Data$Dose)-min(Data$Dose))+min(Data$Dose)),lwd=2,col=2)
lines(density(Data$Dose),lty=1,lwd=2)
legend("topright",legend=c("Original Dose","KAL","LOL","KOL"),lty=c(1,1,1,1),col=c(1,2,3,4))

######Figure 2b
hist(discreteQ_dose_hat*(max(Data$Dose)-min(Data$Dose))+min(Data$Dose),breaks=10,
     main="",xlab="",xlim=c(0,100))
mtext("Dose",side=1,line=2.5,cex=1.7)

######################################################
####plot figure 3: value fuction
######################################################
KAL_real_value <- unlist(read.csv("results/KAL-real-value.csv"))
lol_value_hat <- unlist(read.csv("results/lol-real-value.csv"))
kol_value_hat <- unlist(read.csv("results/kol-real-value.csv"))
discreteQ_value<- unlist(read.csv("results/discreteQ-real-value.csv"))
par(mar=c(3.5,3,0.5,0.5))
par(mfrow=c(1,1))
par(cex=1.7)
plot(density(discreteQ_value),col=6,main="",ylim=c(0,38),lwd=2)
mtext("Estimated value function",side=1,line=2.5, cex=1.7)

lines(density(KAL_real_value),col=2,lwd=2)
lines(density(lol_value_hat,na.rm=TRUE),col=3,lwd=2)
lines(density(kol_value_hat,na.rm=TRUE),col=4,lwd=2)
legend("topleft",lty=c(1,1,1,1),col=c(2,3,4,6),legend=c("KAL","LOL","KOL", "Discretized Q"))
