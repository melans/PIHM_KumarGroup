options(warn=-1)
#install.packages('measurements', 'hydroGOF', repos='http://cran.us.r-project.org', dependencies = TRUE, lib="/glade/u/home/elansary/R/x86_64-pc-linux-gnu-library/3.6")
#install.packages('measurements', 'hydroGOF', repos='http://cran.us.r-project.org', dependencies = TRUE, lib="$HOME/R/lib")

library(measurements)
library(hydroGOF)

args <- commandArgs(TRUE);
dir <- args[1];

# outputs <- file ("outputs.csv",'w')
# writeLines("Run,Temporal,From,To,Spearman,R^2,NSE,mNSE",outputs)

setwd(dir)
# setwd('/home/mr/tamuk/phd/')
# getwd()

obs_data <- paste0(dir,"/obs.data")
sim_data <- paste0(dir,"/somerville.rivFlx1")
dA <- 1009   # drainage area, mi2
start <- '1951-01-01'
cfs <- conv_unit(1,'m','ft') * conv_unit(dA,'mi2', 'ft2') / conv_unit(1,'day','sec') 

#paste0(dir,'rslts.csv')

rslts <- file(paste0(dir,'rslts.csv'),'a')

# writeLines(c(dir),rslts)
# writeLines(paste('','T','spearman2','R2','NSE','mNSE',sep=','),rslts)

d <- read.csv(sim_data, header = FALSE, check.names=FALSE , sep="\t")
sim <- data.frame(date = as.Date(d$V1/24/60 , origin="1949-12-31") , sim = d$V83)
sim <- sim[sim$date >= start,]

d <- read.csv(obs_data, header = FALSE, check.names=FALSE, skip = 33 , sep="\t")
obs <- data.frame(date = as.Date(d$V3) , obs = d$V8)
obs <- obs[obs$date >= start & obs$date <= max(sim$date),]

daily <- data.frame(date = sim$date, obs = obs$obs, sim = sim$sim, sim2 = sim$sim*cfs )

writeLines(paste(
  'Daily',nrow(daily),
  round(cor(daily[,2], daily[,3])^2, 2),
  round(cor(daily[,2], daily[,3] , method = "spearman")^2, 2),
  formatC(NSE(daily[,2], daily[,3]), format="e", digits=2),
  formatC(mNSE(daily[,2], daily[,3]), format="e", digits=2),
  formatC(NSE(daily[,2], daily[,4]), format="e", digits=2),
  formatC(mNSE(daily[,2], daily[,4]), format="e", digits=2),
  dir,sep='\t'
),rslts)

m <- daily
m$month <- format(m$date,format="%m")
m$year <- format(m$date,format="%Y")

m <- aggregate(cbind(m[,2],m[,3],m[,4]) ~ year + month, data = m, mean, na.rm = TRUE)
m$date <- as.yearmon(paste0(m$year,"-",m$month))
monthly <- m[,c(6,3:5)]
monthly <- monthly[with(monthly, order(date)), ]

writeLines(paste(
  'Monthly',nrow(monthly),
  round(cor(monthly[,2], monthly[,3])^2, 2),
  round(cor(monthly[,2], monthly[,3] , method = "spearman")^2, 2),
  formatC(NSE(monthly[,2], monthly[,3]), format="e", digits=2),
  formatC(mNSE(monthly[,2], monthly[,3]), format="e", digits=2),
  formatC(NSE(monthly[,2], monthly[,4]), format="e", digits=2),
  formatC(mNSE(monthly[,2], monthly[,4]), format="e", digits=2),
  dir,sep='\t'
),rslts)

close(rslts)


