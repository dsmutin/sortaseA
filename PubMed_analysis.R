setwd(" ")

##---

if (!require("purrr")) install.packages("purrr"); library(purrr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

##---

convert_files <- function (file_names, N_means = 70, 
                           draw = 1, normalized = T, process_plurals = T,
                           cex = 10, height = 3000, width = 2200, path = "./") {
  
  convert <- function (file_names, N_means, Num, draw = 1) {
    
    Title_n <- str_remove(file_names[Num], ".txt")
    PubMT <- readLines (file_names[Num])
    
    print ("file read, starting conversion")
    
    PubMT <- str_split(PubMT, "- ")
    PubMT <- PubMT[map(PubMT, 1) == "MH  "]
    PubMT <- unlist(map(PubMT, 2))
    PubMT <- str_to_lower(PubMT)
    PubMT <- str_remove_all(PubMT, "\\*")
    PubMT <- str_replace_all(PubMT, ", ", "/")
    PubMT <- str_replace_all(PubMT, "\\\n", "/") 
    PubMT <- str_replace_all(PubMT, "\\\t", "/")
    PubMT <- str_split(PubMT, "/")
    
    PubMT <- unlist(PubMT, use.names = F)
    
    PubMT <- data.frame (table (factor(PubMT)))
    
    if (process_plurals == T) {
      print ("depluralizing...")
      j <- 1
      rm_indexes <- c()
      while (j < length (PubMT[,1])) {
        plural_index <- PubMT[j,1] == str_c(PubMT[,1], "s") |
          PubMT[j,1] == str_c(PubMT[,1], "es") |
          PubMT[j,1] == str_c(str_remove_all(PubMT[,1],"ies"), "y")
        if (is.na(table(plural_index)[2]) == F) {
          PubMT[j,2] <- PubMT[j,2] + PubMT[which(plural_index == T),2]
          rm_indexes <- c (rm_indexes, which(plural_index == T))
        }
        j <- j + 1
      }
      if (length(rm_indexes) > 0) {
        PubMT <- PubMT[-rm_indexes,]
      }
      print ("succesfull")
    }
    
    PubMT <- arrange (PubMT, -Freq)
    colnames (PubMT) <- c ("title", Title_n)
    res <- PubMT[1:N_means,]
    
    
    if (draw == 1) {
      plot_res <- ggplot(data = res, mapping = aes(y = reorder(factor(res[,1]), -res[,2]), x = res[,2])) +
        geom_col(fill = "cadetblue4") +
        labs (title = Title_n) +
        xlab ("Number of articles") +
        ylab("") +
        theme_classic(base_size = cex)
      ggsave(str_c(Title_n, ".png"), plot_res, height = height, 
             width = width, units = "px", path = path)
    }
    
    return(PubMT)
  }
  
  result <- convert (file_names, N_means, 1, draw)
  i <- 2
  
  while (i <= length (files)) {
    result <- left_join (result, convert (file_names, N_means, i, draw), by = "title")
    result[is.na(result[,i+1]),i+1] <- 0
    
    if (draw == 2) {
      
      if (normalized == T) {
        res <- data.frame ("title" = result [, 1], "all" = (result [, 2]/max(result [, 2])),
                           mean = (result [, i+1])/max(result [, i+1]))
        res <- arrange(res, by = -res[,3])
        res <- res[1:N_means,]
        
        plot_res <- ggplot (data = res,  mapping = aes(y = reorder(factor(res[,1]), -res[,2]))) +
          geom_col(data = res, mapping = aes(x = 1),
                   fill = "cadetblue2") +
          geom_col(data = res, mapping = aes(x = res[,3]),
                   fill = "cadetblue4") +
          labs (title = colnames(result [i + 1])) +
          xlab ("normalized % of articles") +
          ylab ("key word") +
          theme_classic(base_size = cex)
        
      } else {
        res <- data.frame ("title" = result [1:N_means, 1], "all" = result [1:N_means, 2],
                           "mean" = (result [1:N_means, i+1]))
        
        plot_res <- ggplot (data = res,  mapping = aes(y = reorder(factor(res[,1]), -res[,3]))) +
          geom_col(data = res, mapping = aes(x = log10(res[,2]+1)),
                   fill = "cadetblue2") +
          geom_col(data = res, mapping = aes(x = log10(res[,3]+1)),
                   fill = "cadetblue4") +
          labs (title = colnames(result [i + 1])) +
          xlab ("log10 fraction of articles") +
          ylab ("key word") +
          theme_classic(base_size = cex)
      }
      ggsave(str_c(colnames(result [i + 1]), ".png"), plot_res, height = height, 
             width = width, units = "px", path = path)
    }
    i <- i + 1
  }
  
  return (result)
}
##---

no_data <- function (Name, Min_all, Min_smp) {
  index <- which(colnames(result) == Name)
  nm <- (result[,2] >= Min_all) & (result[, index] <= Min_smp)
  no_means <- as.character (result[nm,1])
  writeLines(sort(no_means), str_c("no_means_", Name, ".txt"))
  print (str_c("writing no_means_", Name, ".txt"))
}
  
##---

files <- list.files()[str_detect(list.files(), ".txt")]

#draw - to draw .png: 
#  0 - no images, 
#  1 - for raw data, 
#  2 - normalized pictures regarding the files[1]
#normalized - images regarding 100% or the number of articles

result <- convert_files (files, 120, draw = 2, normalized = T, path = "./plots/big", height = 4200)

#Name - name of the term column
#Min_all - minimum quantity of all articles
#Min_smp - minimum quantity of the term articles 
no_data ("output", 200, 10)
