#' Dates for BIMM-143 classes  
#' 
#' Here we setup the dates for our BIMM-143 classes over our 10 week UCSD quarter
#' 

library(lubridate)


## First Tuesday and Thur of term date
tue.1 <- as_date("2018-10-02")
thu.1 <- as_date("2018-10-04")

tue <- tue.1 + weeks(0:9)
thu <- thu.1 + weeks(0:9)

x <- rbind( paste(month(tue), day(tue), sep = "/" ),
    paste(month(thu), day(thu), sep = "/" ))

c(x)

#' **N.B.**: Watch for vacations like Thanksgiving on a Thur in Nov and various Mon vacations
#' 
#' See UCSD accadamic calander for start, end and vacation dates.

#(today())+ days(7)

tue.1 + days(7)
