#load all needed packages
library(ggplot2)
library(pspline)
library(readxl)
library(matrixStats)
library(writexl)
#library(openxlsx)
library(rJava)
library(xlsx)


#INPUT WORKING DIRECTORY
#setwd("C:/Users/Helen/My Drive/Smith/Thesis/graphs/021922")
setwd("C:/Users/chatt/Google Drive/Smith/Thesis/graphs/concentrations")

#INPUT FILE NAME
data_file <- "concentrations input.xlsx"

raw_data <- read_excel(data_file) #load data

#INPUT DATE/TITLE (FOR GRAPHS)
graph_title = "3-04"

#get vector of strain names
strains_dupes <- raw_data$`Time (m)`
strains <- strains_dupes[!duplicated(strains_dupes)]
#get vector of time headers
t <- as.numeric(colnames(raw_data))
t <- t[-1:-1]
t <- signif(t,4) #set sig figs

#INPUT BLANK VALUE
blank = 0.277681667

#INPUT NUMBER OF REPLICATES
replicates <- 2

#define average.data
average.data <- function() {
  row_max <- nrow(raw_data) #get number of rows
  
  #generate iterable for loop
  numbers <- seq(1:(row_max-1))
  numbers_filter <- sapply (numbers, function(x) {x %% replicates == 0})
  iterable <- numbers[numbers_filter]
  iterable <- c(0, iterable) #-1 is first in iterable (for first strain in loop)
  matrix_raw <- data.matrix(raw_data) #convert to matrix
  
  #create empty matrix, with full dimensions
  averages <- matrix(NA, nrow = length(strains), ncol = length(t))
  
  #get average for all strains
  for (counter in iterable){
    if (counter == 0) {matrix_crop <- matrix_raw} #for first strain (iteration)
    else {matrix_crop <- matrix_raw[-1:-counter,]} #trim previous strain, except for first strain
    matrix_crop <- matrix_crop[-(replicates+1):-row_max,-1] #trim to just one strain
    averaged_data <- colMeans(matrix_crop) #get average of each column
    averages[(counter/replicates)+1, ] <- averaged_data # fill in data in row in matrix
  }
  
  #INPUT OUTPUT FILENAME
  output_file_name <- paste(graph_title, "average data.xlsx")
  
  averages <- averages - blank #subtract blank value from all cells
  averages_df <- as.data.frame(averages) #convert to dataframe
  rownames(averages_df) <- strains #change row names
  colnames(averages_df) <- t #change column names
  write.xlsx(averages_df, output_file_name,row.names=TRUE) # write to excel file
  return(averages_df) #return dataframe to global
  
}

#define od.sd
od.sd <- function() {
  row_max <- nrow(raw_data) #get number of rows
  #generate iterable for loop
  numbers <- seq(1:(row_max-1))
  numbers_filter <- sapply (numbers, function(x) {x %% replicates == 0})
  iterable <- numbers[numbers_filter]
  iterable <- c(0, iterable) #-1 is first in iterable (for first strain in loop)
  matrix_raw <- data.matrix(raw_data) #convert to matrix
  
  #INPUT OUTPUT FILENAME
  output_file_name <- paste(graph_title, "od sd.xlsx")

  sd_matrix <- matrix(NA, nrow = length(strains), ncol = length(t)) #create empty matrix, with full dimensions
  straincounter = 0 #create variable to go down through all strains (each row)
  dir.create("OD Error") #create new folder for od error graphs (in working directory)
    
  #get sd and graph od w/ error for all strains
  for (counter in iterable){
    straincounter <- straincounter + 1 # increase counter by one to go to next strain
    if (counter == 0) {matrix_crop <- matrix_raw} #for first strain/iteration
    else {matrix_crop <- matrix_raw[-1:-counter,]} # trim, except for first strain
    matrix_crop <- matrix_crop[-(replicates+1):-row_max,-1] # trim
    sd_data <- colSds(matrix_crop) # get sd
    strain_name <- strains[straincounter] # get strain name
    sd_matrix[straincounter, ] <- t(sd_data) # fill in matrix row with strain data
    
    od_data <- averaged_data[straincounter,] #get average ods for strain
    min_error <- od_data - sd_data #get min error values
    max_error <- od_data + sd_data #get max error values
    colnames(min_error)=t #fix match.names error
    colnames(max_error)=t
    colnames(od_data)=t
    od_df <- rbind.data.frame(time = t, od = od_data, min = min_error, max = max_error) #create df
    od_df <- data.frame(t(od_df)) #transpose (flip rows and columns)
    rownames(od_df) <- c() #delete row names
    
    #graph OD error over time
    ggplot(od_df, aes(x=time, y=od)) + 
      geom_ribbon(
        aes(ymin=min, ymax=max), fill="grey70") +
      geom_line(aes(y=od), lwd=0.1) +
      labs(title=paste(strain_name, "OD Error,", graph_title),
           x = "Time (min)",
           y = "OD")
    #save graph as png
    ggsave(paste0("OD error/",strain_name,".png"),
           width=20, height=10, unit="cm")
  }
  
  sd_df <- data.frame(sd_matrix) # convert to df
  rownames(sd_df) <- strains #row names are strains
  colnames(sd_df) <- t #column names are time
  write.xlsx(sd_df,output_file_name, row.names=TRUE, append=TRUE) #write to new sheet in excel
}

#define get.deriv
get.deriv <- function() {
  
  #INPUT OUTPUT FILE NAME
  output_file_name <- "derivs.xlsx"
  
  row <- colnames(averaged_data_flipped) #iterable for loop
  max_rates <- matrix(NA, nrow = length(strains), ncol = 2) #create empty matrix to track max rates
  row_num = 0 # counter for matrix row
  
  for (strain in row) {
    row_num = row_num + 1
    x <- averaged_data_flipped[ ,strain] #get od values for strain
    #penalize spline to fourth derivative (to later get 2nd)
    smoothed <- smooth.Pspline(t,x, norder=4)
    deriv1 <- predict(smoothed, t, nderiv=1) #first deriv
    deriv2 <- predict(smoothed, t, nderiv=2) #second deriv
    max_time_data <- deriv1[-1:-1] #first derivs, minus the first value
    max_x<-which.max(max_time_data) #position of max, -1
    max_growth <- signif(deriv1[max_x+1], digits=4) #get max growth
    max_time <- (max_x) * (16+(5/8)) #get time of max
    
    #get dataframe of derivatives and export
    derivs <- as.data.frame(cbind(t, deriv1=c(deriv1), deriv2=c(deriv2)))
    
    write.xlsx(derivs, output_file_name, sheetName=strain, col.names=TRUE, row.names=FALSE, append=TRUE)
    
    max_rate_row = c(max_growth, max_time) # create list of strain and its max rate
    max_rates[row_num,] <- max_rate_row #add strain's max rate and time to df
  }
  
  #export df of max rates as new sheet
  max_rates <- data.frame(max_rates) # convert to df
  colnames(max_rates) <- c("Max Growth Rate", "Max Rate Time") # rename columns
  rownames(max_rates) <- strains # rename rows
  write.xlsx(max_rates, output_file_name, sheetName="Max Rates", row.names=TRUE, append=TRUE)
  
  return (output_file_name) #return file name of deriv data
}
  
  
#define deriv.sd
deriv.sd <- function() {
  row_max <- nrow(raw_data) #get number of rows
  
  #subtract blank
  raw_data_flipped[,1:row_max] <- lapply(raw_data_flipped[,1:row_max], function(x) {as.numeric(as.character(x))})
  raw_data_flipped <- raw_data_flipped - blank
  
  #generate iterable for loop
  numbers <- seq(1:(row_max-1))
  numbers_filter <- sapply (numbers, function(x) {x %% replicates == 0})
  iterable <- numbers[numbers_filter]
  iterable <- c(0, iterable) #0 is first in iterable (for first strain in loop)

  deriv1s <- matrix(NA, nrow = length(strains_dupes), ncol = length(t)) # create empty matrix for 1st derivs
  deriv2s <- matrix(NA, nrow = length(strains_dupes), ncol = length(t)) # create empty matrix for 2nd derivs
  row_num = 0 # for matrix row
  
  #loop for all strains
  for (counter in iterable){
    #trim to just that one strain
    if (counter == 0) {data_cropped <- raw_data_flipped} # for first strain
    else {data_cropped <- raw_data_flipped[,-replicates:-counter]} # for all other strains - trim previous strains
    data_cropped <- data_cropped[,-(replicates+1):-row_max]
    #get deriv values for each well
    for (column in 1:replicates) {
      row_num = row_num + 1
      od <- data_cropped[,column] #get OD values for first well
      smoothed <- smooth.Pspline(t,od, norder=4) #smooth w/ pspline, penalize to 4th deriv
      deriv1 <- predict(smoothed, t, nderiv=1) #first deriv
      deriv2 <- predict(smoothed, t, nderiv=2) #second deriv
      deriv1s[row_num, ] <- c(deriv1) #fill matrix row with first derivs
      deriv2s[row_num, ] <- c(deriv2) # fill matrix row with second derivs
    }
  }
  
  #save derivative value for each run (not averages)
  deriv1s_df <- data.frame(deriv1s) #make into df
  deriv2s_df <- data.frame(deriv2s) 
  colnames(deriv1s_df) <- t #rename columns to time
  colnames(deriv2s_df) <- t
  rownames(deriv1s_df) <- make.names(strains_dupes, unique=TRUE) #rename rows to strains (dupes)
  rownames(deriv2s_df) <- make.names(strains_dupes, unique=TRUE)
  write.xlsx(deriv1s_df, "1st derivative all wells.xlsx") #write to new excel file
  write.xlsx(deriv2s_df, "2nd derivative all wells.xlsx")
  
  #INPUT OUTPUT FILE NAME
  output_file_name <- "derivative sd.xlsx"
  
  #get derivative sds for all strains
  for (counter in iterable){
    #1st deriv
    if (counter == 0) {deriv1_crop <- deriv1s  # for first strain
    } else {deriv1_crop <- deriv1s[-1:-counter,]} # trim previous strains, not for first strain
    deriv1_crop <- deriv1_crop[-(replicates+1):-row_max,] # trim to one strain
    deriv1_sd_data <- colSds(deriv1_crop) # get sd for each column
    #2nd deriv
    if (counter == 0) {deriv2_crop <- deriv2s} # for first strain
    else {deriv2_crop <- deriv2s[-1:-counter,]} # trim previous strains, not for first strain
    deriv2_crop <- deriv2_crop[-(replicates+1):-row_max,] #trim
    deriv2_sd_data <- colSds(deriv2_crop) #get sd
    #
    sd_df <- data.frame(t, deriv1_sd_data, deriv2_sd_data) #make dataframe
    colnames(sd_df) <- c("time", "1st deriv sd", "2nd deriv sd")
    strain_name <- strains_dupes[counter+1] #get strain name
    write.xlsx(sd_df,output_file_name, sheetName = as.character(strain_name), row.names=FALSE, append=TRUE) #write to new sheet in excel
  }
  
  return (output_file_name) #return file name of deriv sd xlsx
}

#define deriv.sd.graph
deriv.sd.graph <- function() {
  
  iterable <- (1:length(strains)) #create new iterable
  
  dir.create("1st deriv error graphs") #create new folder for 1st deriv error graphs in wd
  dir.create("2nd deriv error graphs") #create new folder for 2nd deriv error graphs in wd
  
  #loop for all strains (each sheet in excel file)
  for (counter in iterable) {
    deriv_data <- read_excel(deriv_data_file, sheet = counter) #load data for 1st and 2nd derivs
    deriv_sd_data <- read_excel(deriv_sd_file, sheet = counter) #load data for deriv sds
    deriv1 <- deriv_data$deriv1 #get 1st deriv values
    deriv2 <- deriv_data$deriv2 #get 2nd deriv values
    
    min1 <- deriv1 - deriv_sd_data$`1st deriv sd` #get 1st deriv min error
    max1 <- deriv1 + deriv_sd_data$`1st deriv sd` #get 1st deriv max error
    min2 <- deriv2 - deriv_sd_data$`2nd deriv sd` #get 2nd deriv min error
    max2 <- deriv2 + deriv_sd_data$`2nd deriv sd` #get 2nd deriv max error
    strain_name <- strains[counter] #get strain name
    
    max_rate_time <- which.max(deriv1[-1:-1]) #position of max growth rate, minus the first value
    full_time <- length(t) #get number of time points
    t_crop <- t#[-(max_rate_time+10:full_time)] #time values, 10 past the max growth rate, -1
    deriv1_crop <- deriv1#[-(max_rate_time+10:full_time)] #1st deriv values, 10 past max, -1
    min1_crop <- min1#[-(max_rate_time+10:full_time)] #1st deriv min error values, 10 past, -1
    max1_crop <- max1#[-(max_rate_time+10:full_time)] #1st deriv max error values, 10 past, -1
    deriv2_crop <- deriv2#[-(max_rate_time+10:full_time)] #2nd deriv values, 10 past max, -1
    min2_crop <- min2#[-(max_rate_time+10:full_time)] #2nd deriv min error values, 10 past, -1
    max2_crop <- max2#[-(max_rate_time+10:full_time)] #2nd deriv max error values, 10 past, -1
    
    deriv1_error_data <- data.frame(cbind(t_crop, deriv1_crop, min1_crop, max1_crop)) #new df for 1st deriv error
    deriv2_error_data <- data.frame(cbind(t_crop, deriv2_crop, min2_crop, max2_crop)) #new df for 2nd deriv error
    
    marker <- (max_rate_time) * t[2] #get time of max, 997=time interval
    max_growth <- signif(deriv1[max_rate_time+1], digits=4) #get max growth rate
    
    #graph 1st deriv error
    ggplot(deriv1_error_data, aes(x=t_crop, y=deriv1_crop)) + 
      geom_ribbon(
        aes(ymin=min1_crop, ymax=max1_crop), fill="grey70") +
      geom_line(aes(y=deriv1_crop), lwd=0.1) +
      labs(title=paste(strain_name, "1st derivative error,", graph_title),
           x = "Time (min)",
           y = "1st derivative") +
      geom_vline(xintercept=marker) +
      geom_text(mapping = aes(x=marker+40, y=0, label=max_growth))
    #save 1st deriv error graph as png
    ggsave(paste0("1st deriv error graphs/",strain_name,".png"),
           width=20, height=10, unit="cm")
    #graph 2nd deriv error
    ggplot(deriv2_error_data, aes(x=t_crop, y=deriv2_crop)) + 
      geom_ribbon(
        aes(ymin=min2_crop, ymax=max2_crop), fill="grey70") +
      geom_line(aes(y=deriv2_crop), lwd=0.1) +
      labs(title=paste(strain_name, "2nd derivative error,", graph_title),
           x = "Time (min)",
           y = "2nd derivative") +
      geom_vline(xintercept=marker)
    #save 2nd deriv error graph as png
    ggsave(paste0("2nd deriv error graphs/",strain_name,".png"),
           width=20, height=10, unit="cm")
  
    }
}

get.lagtime <- function() {
  strain_iterable <- (1:length(strains)) #get iterable to loop through xl file sheets
  lagtimes <- matrix(NA, nrow = length(strains), ncol = 4) #create empty matrix
  dir.create("Lag time graphs") #create folder for graphs in wd
  max_rate_df <- read_excel(deriv_data_file, sheet="Max Rates") #load max growth rate data
  max_rates <- max_rate_df$`Max Growth Rate` #isolate max growth rate column
  row_num = 0
  
  #loop for all strains
  for (strain_counter in strain_iterable){
    row_num = row_num + 1
    data <- read_excel(deriv_data_file, sheet = strain_counter) #load deriv data for one strain
    deriv1 <- data$deriv1 #get 1st deriv data
    deriv2 <- data$deriv2 #get 2nd deriv data
    max_rate <- as.numeric(max_rates[strain_counter]) #get max growth rate for strain
    od_data <- t(averaged_data[strain_counter,]) #get od data
    iterable <- 1:(length(deriv2)-2) #set iterable for datapoints, -2
    #loop to find where two consecutive measurements are positive, and growth rate > 15% of max
    for (measure in iterable){
      if (deriv2[measure] > 0 & deriv2[measure+1] > 0 &
          #od_data[measure] > .03 & od_data[measure+1] > 0.03 & od_data[measure+2] > 0.03) 
          deriv1[measure] > (max_rate)*0.15) {
        lagtime <- t[measure]
        break
      }
    }
    lagtimes[row_num, ] <- c(lagtime, od_data[measure], deriv1[measure], deriv2[measure]) # fill in matrix row with info
    
    lag_od_df <- data.frame(cbind(t,od_data)) #make into df
    lagtime_label <- round(lagtime, digits=1) #round lag time to 2 digits, for label
    strain_name <- strains[strain_counter] #get strain name
    
    #graph od data with lagtime marker
    ggplot(lag_od_df, aes(x=t, y=od_data)) + 
      geom_line(aes(y=od_data), lwd=.8) +
      labs(title=paste(strain_name, "Lag Time,", graph_title),
           x= "Time (min)",
           y= "OD") +
      geom_vline(xintercept=lagtime) +
      geom_text(mapping = aes(x=lagtime+50, y=.5, label = lagtime_label))
    #save graph as png
    ggsave(paste0("Lag time graphs/",strain_name,".png"),
           width=20, height=10, unit="cm")
    }
  
  lagtimes <- data.frame(lagtimes) # convert to df
  colnames(lagtimes) <- c("Lag time (min)", "OD", "1st deriv", "2nd deriv") #rename columns
  rownames(lagtimes) <- strains # rename rows
  
  #INPUT OUTPUT FILE NAME FOR LAG TIME
  output_file_name <- "lag time.xlsx"
  
  write.xlsx(lagtimes, output_file_name) #export to excel file
  
}

get.lagtime.all.wells <- function() {  
  # get iterable to loop through wells
  row_max <- nrow(raw_data)
  numbers <- seq(1:(row_max-1))
  numbers_filter <- sapply (numbers, function(x) {x %% replicates == 0})
  iterable_raw <- numbers[numbers_filter]
  iterable_raw <- c(0, iterable_raw) #-1 is first in iterable (for first strain in loop)
  iterable <- iterable_raw + 1
  
  # get iterable to loop through strains
  strain_iterable_single <- (1:(length(row_max)/3)) 
  strain_iterable <- sort(c(strain_iterable_single, strain_iterable_single, strain_iterable_single))
  
  lagtimes <- matrix(NA, nrow = length(strains), ncol = 6) # create empty matrix
  dir.create("Lag time graphs all wells") # create folder for graphs in wd
  max_rate_df <- read_excel(deriv_data_file, sheet="Max Rates") # load max growth rate data
  max_rates <- max_rate_df$`Max Growth Rate` # isolate max growth rate column
  
  deriv1_data <- read_excel("1st derivative all wells.xlsx") # open data file
  deriv2_data <- read_excel("2nd derivative all wells.xlsx")
  
  #loop for all wells
  for (counter in iterable) {
    strain_counter = ((counter-1) / replicates) + 1
    od_data <- t(averaged_data[strain_counter,]) #get od data
    for (loop in (1:replicates)) {
      counter2 = (counter + loop) - 1
      deriv1 <- t(deriv1_data[counter2, ]) #get 1st deriv data
      deriv1 <- deriv1[-1,]
      deriv2 <- t(deriv2_data[counter2, ]) #get 2nd deriv data
      deriv2 <- deriv2[-1,]
      max_rate <- as.numeric(max_rates[strain_counter]) #get max growth rate for strain
      measurements <- 1:(length(deriv2)-2) #set iterable for datapoints, -2
      #loop to find where two consecutive measurements are positive, and growth rate > 15% of max
      for (measure in measurements) {
        if (deriv2[measure] > 0 & deriv2[measure+1] > 0 &
            deriv1[measure] > (max_rate)*0.15) {
          assign(paste0("lagtime", loop), t[measure])
          break
        }
      }
    }
    
    lagtimes[strain_counter, ] <- c(lagtime1, lagtime2, lagtime3, od_data[measure], deriv1[measure], deriv2[measure]) # fill in matrix row with info
    
    lag_od_df <- data.frame(cbind(t,od_data)) #make into df
    lagtime_label1 <- round(lagtime1, digits=1) #round lag time to 2 digits, for label
    lagtime_label2 <- round(lagtime2, digits=1)
    lagtime_label3 <- round(lagtime3, digits=1)
    
    strain_name <- strains[strain_counter] #get strain name
    
    #graph od data with lagtime marker
    ggplot(lag_od_df, aes(x=t, y=od_data)) + 
      geom_line(aes(y=od_data), lwd=.8) +
      labs(title=paste(strain_name, "Lag Time,", graph_title),
           x= "Time (min)",
           y= "OD") +
      geom_vline(colour = "blue4", xintercept=lagtime1) +
      geom_vline(colour = "red4", xintercept=lagtime2) +
      geom_vline(colour = "green4", xintercept=lagtime3) +
      geom_text(colour = "blue4", mapping = aes(x=lagtime1+50, y=.7, label = lagtime_label1)) +
      geom_text(colour = "red4", mapping = aes(x=lagtime2+50, y=.5, label = lagtime_label2)) +
      geom_text(colour = "green4", mapping = aes(x=lagtime3+50, y=.3, label = lagtime_label3))
    #save graph as png
    ggsave(paste0("Lag time graphs all wells/", strain_name,".png"),
           width=20, height=10, unit="cm")
  }
  
  lagtimes <- data.frame(lagtimes) # convert to df
  colnames(lagtimes) <- c("Lag time 1 (min)", "Lag 2", "Lag 3", "OD", "1st deriv", "2nd deriv") #rename columns
  rownames(lagtimes) <- strains # rename rows
  
  #INPUT OUTPUT FILE NAME FOR LAG TIME
  output_file_name <- "lag time all wells.xlsx"
  
  write.xlsx(lagtimes, output_file_name) #export to excel file

}



#call functions:

averaged_data <- average.data() #call average.data function, return df
averaged_data_flipped <- data.frame(t(averaged_data)) #flip columns and rows

od.sd() #call od.sd function

deriv_data_file <- get.deriv() #call get.deriv function

raw_data_flipped <- data.frame(t(raw_data)) #transpose
colnames(raw_data_flipped) <- raw_data_flipped[1,] #rename columns
raw_data_flipped <- raw_data_flipped[-1,] #delete

deriv_sd_file <- deriv.sd() #call deriv.sd function

deriv.sd.graph() #call deriv.sd.graph function

get.lagtime() #call get.lagtime function

# get.lagtime.all.wells() #call get.lagtime.all.wells function
