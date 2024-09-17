set.seed(123)
library("stochvol")
library(readxl)

#lee los datos
realdatadtf <- read_excel("E:/Universidad/PI 3/Monedas.xlsx") #cambiar a donde este la carpeta
realnames <- names(realdatadtf)[3:16] #obtiene los nombres de las manedas usadas
boxplotnames <- names(realdatadtf)[2:16]
boxplotnames <- sub("/USD","",boxplotnames)

realdata <- data.matrix(realdatadtf[,2:16])
days <- dim(realdata)[1]
currency <- dim(realdata)[2]

#saca los log retornos
logretdata <- matrix(nrow=days-1,ncol=currency)
for (i in 1:currency) {
  logretdata[,i] <- logret(realdata[,i])
}

#obtiene los parametros de h_t y guarda su media
h_t2<- matrix(nrow=days-1,ncol=currency)
for (i in 1:currency) {
  res <- svsample(logretdata[,i],draws = 18000,burnin = 3000, priormu = c(-10, 1), priorphi = c(20, 1.1), priorsigma = .1)
  h_t2[,i] <-res[["summary"]][["latent"]][,6]
}
i=1

#Crea las distribuciones por dia
x <- seq(-0.0005, 0.0005, by = 0.0000001)
norm <- matrix(nrow=length(x),ncol=currency*(days-1))
for (i in 1:currency) {
  for (j in 1:(days-1)) {
    norm[,(i-1)*(days-1)+j] <-t(dnorm(x,0,h_t2[j,i]))
  }
}

i=1
j=1
dim(norm)

#calcula el KL
library(LaplacesDemon)
KLlist <- matrix(nrow=1,ncol=currency-1)


for (i in 2:currency) {
  auxsum=0
  for (j in 1:(days-1)) {
    KLdist <- KLD(norm[,j],norm[,(i-1)*(days-1)+j])
    auxsum=auxsum+KLdist$sum.KLD.px.py
  }
  KLlist[,i-1] <- auxsum
}
class(KLdist)

names <- names(realdata)[2:currency]

#DataFrame con el KL a cada moneda
result <-as.data.frame(KLlist)
colnames(result) <- realnames
result <- as.data.frame(result)
best=min(result)

logretdatadt=as.data.frame(logretdata)
colnames(logretdatadt) <- boxplotnames
format(logretdatadt, scientific = FALSE)


library(writexl)
write_xlsx(logretdatadt,"E:/Universidad/PI 3/table logret.xlsx")
#paso a mano
plot(realdatadtf$Date[1:367],logret(data.matrix(realdatadtf$`EUR/USD`)),type="l",col="red",main = "Volatility EUR vs JPY",xlab="Date",ylab = "LogReturn")
lines(realdatadtf$Date[1:367],logret(data.matrix(realdatadtf$`JPY/USD`)),col="green")



#Box plots
logretdata2=logretdata[,-6]
boxplotnames2=boxplotnames[-6]

par(cex.main=2)
par(cex.axis=1.25)
par(cex.lab=1.7)
boxplot(logretdata2,names = boxplotnames2,xlab="Currency",ylab="Log-Returns",main = "Total")
boxplot(logretdata2[1:61,],names = boxplotnames2,xlab="Currency",ylab="Log-Returns",main = "Boxplot first trimester")
boxplot(logretdata2[62:122,],names = boxplotnames2,xlab="Currency",ylab="Log-Returns",main = "Boxplot second trimester")
boxplot(logretdata2[123:183,],names = boxplotnames2,xlab="Currency",ylab="Log-Returns",main = "Boxplot third trimester")
boxplot(logretdata2[184:244,],names = boxplotnames2,xlab="Currency",ylab="Log-Returns",main = "Boxplot fouth trimester")
boxplot(logretdata2[245:305,],names = boxplotnames2,xlab="Currency",ylab="Log-Returns",main = "Boxplot fifth trimester")
boxplot(logretdata2[306:366,],names = boxplotnames2,xlab="Currency",ylab="Log-Returns",main = "Boxplot sixth trimester")

#Trimestral tabla
KLlist_One <- matrix(nrow=1,ncol=currency-1)
for (i in 2:currency) {
  auxsum=0
  for (j in 1:61) {
    KLdist_One <- KLD(norm[,j],norm[,(i-1)*(days-1)+j])
    auxsum=auxsum+KLdist_One$sum.KLD.px.py
  }
  KLlist_One[,i-1] <- auxsum
}

KLlist_Two <- matrix(nrow=1,ncol=currency-1)
for (i in 2:currency) {
  auxsum=0
  for (j in 62:122) {
    KLdist_Two <- KLD(norm[,j],norm[,(i-1)*(days-1)+j])
    auxsum=auxsum+KLdist_Two$sum.KLD.px.py
  }
  KLlist_Two[,i-1] <- auxsum
}

KLlist_Three <- matrix(nrow=1,ncol=currency-1)
for (i in 2:currency) {
  auxsum=0
  for (j in 123:183) {
    KLdist_Three <- KLD(norm[,j],norm[,(i-1)*(days-1)+j])
    auxsum=auxsum+KLdist_Three$sum.KLD.px.py
  }
  KLlist_Three[,i-1] <- auxsum
}

KLlist_Four <- matrix(nrow=1,ncol=currency-1)
for (i in 2:currency) {
  auxsum=0
  for (j in 184:244) {
    KLdist_Four <- KLD(norm[,j],norm[,(i-1)*(days-1)+j])
    auxsum=auxsum+KLdist_Four$sum.KLD.px.py
  }
  KLlist_Four[,i-1] <- auxsum
}

KLlist_Five <- matrix(nrow=1,ncol=currency-1)
for (i in 2:currency) {
  auxsum=0
  for (j in 245:305) {
    KLdist_Five <- KLD(norm[,j],norm[,(i-1)*(days-1)+j])
    auxsum=auxsum+KLdist_Five$sum.KLD.px.py
  }
  KLlist_Five[,i-1] <- auxsum
}

KLlist_Six <- matrix(nrow=1,ncol=currency-1)
for (i in 2:currency) {
  auxsum=0
  for (j in 306:366) {
    KLdist_Six <- KLD(norm[,j],norm[,(i-1)*(days-1)+j])
    auxsum=auxsum+KLdist_Six$sum.KLD.px.py
  }
  KLlist_Six[,i-1] <- auxsum
}
tableKL <- matrix(c(boxplotnames[2:15], KLlist_One, KLlist_Two, KLlist_Three, KLlist_Four, KLlist_Five, KLlist_Six, KLlist),ncol=8)

tableKLdf <- as.data.frame(tableKL)
  

format(tableKLdf, scientific = FALSE)
colnames(tableKLdf)<- c("Names","First", "Second", "Third","Fourth", "Fifth","Sixth","Total")

library(writexl)
write_xlsx(tableKLdf,"E:/Universidad/PI 3/table KL.xlsx")

correlat <- matrix(nrow=7,ncol=currency-1)

#Correlacion
for (i in 2:currency-1) {
  correlat[1,i]= cor(logretdata2[1:61,1],logretdata2[1:61,i])
}

for (i in 2:currency-1) {
  correlat[2,i]= cor(logretdata2[62:122,1],logretdata2[62:122,i])
}
for (i in 2:currency-1) {
  correlat[3,i]= cor(logretdata2[123:183,1],logretdata2[123:183,i])
}
for (i in 2:currency-1) {
  correlat[4,i]= cor(logretdata2[184:244,1],logretdata2[184:244,i])
}
for (i in 2:currency-1) {
  correlat[5,i]= cor(logretdata2[245:305,1],logretdata2[245:305,i])
}
for (i in 2:currency-1) {
  correlat[6,i]= cor(logretdata2[306:366,1],logretdata2[306:366,i])
}
for (i in 2:currency-1) {
  correlat[7,i]= cor(logretdata2[,1],logretdata2[,i])
}

tablecor <- matrix(c(boxplotnames2,t(correlat)),ncol = 8)
tablecordf <- as.data.frame(tablecor)


format(tablecordf, scientific = FALSE)
colnames(tablecordf)<- c("Names","First", "Second", "Third","Fourth", "Fifth","Sixth","Total")

write_xlsx(tablecordf,"E:/Universidad/PI 3/table cor.xlsx")

plot(realdatadtf$Date[1:367],logret(data.matrix(realdatadtf$`EUR/USD`)),type="l",col="red",main = "Volatility EUR vs NOK",xlab="Date",ylab = "LogReturn")
lines(realdatadtf$Date[1:367],logret(data.matrix(realdatadtf$`NOK/USD`)),col="green")
