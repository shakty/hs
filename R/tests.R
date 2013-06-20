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

tic()
fit <- lm(

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


        
