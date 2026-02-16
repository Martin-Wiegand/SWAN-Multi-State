library("tidyverse")

main_path <- "D:/Work/OneDrive - University College London/SWAN/JMstateModel-master/JMstateModel-master/"

data_path <- paste0(main_path,"Example/")
code_path <- main_path
rds_path <- paste0(main_path,"RDS/") 
plot_path <- paste0(main_path,"Plots/")

# csv_file <- function(x){paste0(data_path,x,".csv")}
code_file <- function(x)paste0(code_path,x)
rds_file <- function(x)paste0(rds_path,x,".rds")
plot_file <- function(x)paste0(plot_path,x)
data_file <- function(x){paste0(data_path,x)}

# source(code_file("0-Utils"))
