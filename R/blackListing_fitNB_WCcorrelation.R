library(MASS)
library(stats4)
library(ggplot2)

y = read.table("firstMateCoverage-L100000.bed") # count reads-- first mate non-dup passed QC
z = read.table("firstMateUnqCoverage-L100000.bed") # count reads-- unique first mate non-dup passed QC

head = y[1,]

# removing the first row
y = y[2:dim(y)[1],]
z = z[2:dim(z)[1],]

numCells = (dim(y)[2] - 3)/2

RC = matrix(nrow = dim(y)[1], ncol = numCells) # Column j contains the coverage of bins for cell j
RCAllCells = rep(0,dim(y)[1]) # count the total coverage of all cells in each bin

RCunq = matrix(nrow = dim(z)[1], ncol = numCells) # Column j contains the number of unique reads (passed QC quality) mapped to each bin from cell j
RCunqAllCells = rep(0,dim(z)[1]) # count the total number of reads (passed QC quality) from all cells mapped to each bin

WC = matrix(nrow = dim(y)[1], ncol = numCells) # Column j contains the Watson coverage of bins for cell j
CC = matrix(nrow = dim(y)[1], ncol = numCells) # Column j contains the Crick coverage of bins for cell j

for (i in 1:numCells)
{
  WC[,i] = as.numeric(as.character(y[,2*i+2]))
  CC[,i] = as.numeric(as.character(y[,2*i+3]))
  
  RC[,i] = WC[,i] + CC[,i]
  RCAllCells = RCAllCells + RC[,i]
  
  RCunq[,i] = as.numeric(as.character(z[,2*i+2])) + as.numeric(as.character(z[,2*i+3]))
  RCunqAllCells = RCunqAllCells + RCunq[,i]
}

#blacklisting
# excluding zero bins
RC = RC[which(RCAllCells != 0),]
WC = WC[which(RCAllCells != 0),]
CC = CC[which(RCAllCells != 0),]
RCunq = RCunq[which(RCAllCells != 0),]
y = y[which(RCAllCells != 0),]
RCunqAllCells = RCunqAllCells[RCAllCells != 0]
RCAllCells = RCAllCells[RCAllCells != 0]

# subset bins with low fraction of multiply mapped reads
multFrac = (RCAllCells - RCunqAllCells) / RCAllCells
alpha = 0.05 # the cutoff fraction of multiply mapped reads
hist(multFrac)
length(multFrac[multFrac < alpha])

RCAllCells = RCAllCells[multFrac < alpha]
RC = RC[which(multFrac < alpha),]
WC = WC[which(multFrac < alpha),]
CC = CC[which(multFrac < alpha),]
y = y[which(multFrac < alpha),]

# output blacklisted bins to a bed file

y = rbind(head, y)
write.table(y, "normalBinsCoverage-L100000.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# mean-var plot

avg = NULL # sample means
s = NULL # sample vars

for (i in 1:numCells)
{
  avg = c(avg, mean(RC[,i]))
  s = c(s, var(RC[,i])) 
}

plot(avg, s, xlab = "mean", ylab = "var")
c = sum(s*avg)/sum(avg*avg)
lines(avg, c*avg, col = "red")

# estime p based on the coefficient of this line fitting

estP = 1/c # estimated p
estR = NULL

estR = avg*estP/(1-estP)

pdf("readCount.pdf")

par(mfrow=c(2,2))
for (i in 1:numCells)
{
  hist(RC[,i],breaks = (-1:max(RC[,i]))+0.5, freq = FALSE, xlab = "", ylab = "", 
       main = paste("cell",i))
  x = 0:max(RC[,i])
  y = dnbinom(x, size = estR[i], prob = estP)
  lines(x,y,col="red")
  
  p = 1:99/100
  plot(qnbinom(p, size = estR[i], prob = estP), quantile(RC[,i], probs = p), 
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
  lines(c(0,max(RC[,i])), c(0,max(RC[,i])), col = "red")
}
dev.off()

# output the estimated parameters in a file

sink("NBpars.data")
cat(paste("p\t", estP, sep = ""))
for (i in 1:length(estR))
  cat(paste("\nr", i, "\t", estR[i], sep = ""))
sink()


# W and C read counts relation

pdf("WCrelation.pdf")

plots <- list()
i = 1
for (i in 1:numCells)
  local({
    i <- i
    p1 <<- ggplot(data.frame(WC[,i], CC[,i]), aes(WC[,i], CC[,i])) + 
      xlab(paste("W", i)) + ylab(paste("C", i)) + geom_point() + geom_density2d()
    plots[[i]] <<- p1
  })

for(i in 1:numPages)
{
  grid.arrange(grobs = plots[((i-1)*6+1): min(i*6,numCells)], nrow = 3, ncol = 2)
}

dev.off()
