# BIGMEMORY
library(bigmemory)
library(biganalytics)


x <- read.big.matrix("clusters_micro_smaller.csv", header=TRUE,
                     backingfile="micro_s.bin",
                     descriptorfile="micro_s.desc")
library(sqldf)
tic()
f1 <- file('clusters_micro_smaller.csv')
toc()

tic()
mmicro <- sqldf("select * from f1", dbname = tempfile(), file.format = list(header = T, row.names = F))
toc()

a <- read.csv.sql("clusters_micro.csv", sql = "select * from file")

library(RSQLite)
con <- dbConnect("SQLite", dbname = "sample_db")

        
# read csv file into sql database
dbWriteTable(con, name="sample_data", value="clusters_macro.csv", 
             row.names=FALSE, header=TRUE, sep = ",")


library(ffbase)

tic()
ffx <- read.csv.ffdf(file="clusters_micro", header=TRUE)
toc()
          


        
a <- c(1,2,3)
b <- c(10,2,1)
c <- c(100,9,2)
all <- c(a,b,c)
mean(all)

all <- c(1,1,5,9,9)

all <- c(rep(29,20),c(32,27))



all2 <- all^2

sqrt((sum(all2) - ((sum(all))^2/length(all)))/(length(all)-1))

sd(all)

mean(c(mean(a),mean(b),mean(c)))
sd(c(sd(a),sd(b),sd(c)))
          

          
          
sqrt((18573 - (639/22))/21)
