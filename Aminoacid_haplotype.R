# Load libraries
library(dplyr)
library(tidyr)

# Set your work directory to the folder containing your data
setwd("C:/YOUR/WORK/DIRECTORY") 
df <- read.csv(file= "YOUR_DATA.csv", stringsAsFactors=T)

## PART 1 ##
# Sort sequences based on pattern of polymorphic site
# reference_exist = TRUE, if there is a reference sequence.
#   If reference sequence exist in the data, ensure that the reference 
#   sequence is in the first entry row
# group_exist = TRUE, to define population grouping
# sort_by_group = TRUE, sort the data by population group
# By default, the arguments are all FALSE

variation_sort <- function(x, reference_exist = FALSE, group_exist = FALSE, sort_by_group = FALSE){
  if (reference_exist){
    if (!group_exist) {
      print("Group does not exist")
      x <- x %>% arrange(x[2:nrow(x),2:length(x)])
    }
    if (group_exist) {
      print("Groups info available")
      if (!sort_by_group) {
        print("Do not sort by group")
        x <- x %>% arrange(x[2:nrow(x),3:length(x)])
      }
      if (sort_by_group) {
        print("Sorted by group")
        x <- x %>% arrange((x[2:nrow(x),2]),x[2:nrow(x),3:length(x)])
      }
    }  
  }
  if (!reference_exist){
    if (!group_exist) {
      print("Group does not exist")
      x <- x %>% arrange(x[,2:length(x)])
    }
    if (group_exist) {
      print("Groups info available")
      if (!sort_by_group) {
        print("Do not sort by group")
        x <- x %>% arrange(x[,3:length(x)])
      }
      if (sort_by_group) {
        print("Sorted by group")
        x <- x %>% arrange((x[,2]),x[,3:length(x)])
      }  
    }
  }
  return(x)
}

sorted_data <- variation_sort(df, reference_exist = F,
                              group_exist = F, sort_by_group = F)
sorted_data
write.csv(sorted_data, "outputs/sorted_data.csv", row.names=F)

## PART 2 ##
# List unique haplotypes and their frequency
# Important: this step is case sensitive
# ref_seq_as_H1 =  TRUE, if there is a reference sequence. The reference 
#   sequence and all sequences with same haplotype will be labelled as "H1" 
   

haplocount <- function(x, ref_seq_as_H1 = FALSE, group_exist = FALSE, by_group = FALSE){
  haplofunct1 <- function(selected_df){
    comb <- data.frame(lapply(selected_df[1:(length(selected_df)-1)], 
                              function(x) paste0(x,"$")))
    comb2 <- cbind(comb, selected_df[,length(selected_df)])
    comb3<- data.frame(table(apply(comb2, 1, paste,  collapse = "")))
    colnames(comb3) <- c("Combination", "Frequency")
    unique_comb <- comb3 %>%
      separate("Combination", sep = "\\$", into = colnames(selected_df))
    return(unique_comb)
  }
  if (!group_exist) {
    print("Group does not exist")
    selected_df <- x[,2:length(x)]
    analyzed <- haplofunct1(selected_df)
    output <- analyzed %>% arrange(desc(analyzed$Frequency))
    first_matching_row <- x[1,2:length(x)] # for ref seq later
  }
  if (group_exist) {
    print("Group info available")
    if (!by_group) {
      print("Do not analyze by group")
      selected_df <- x[,3:length(x)]
      analyzed <- haplofunct1(selected_df)
      output <- analyzed %>% arrange(desc(analyzed$Frequency))
      first_matching_row <- x[1,3:length(x)] # for ref seq later
    }
    if (by_group) {
      print("Analyzed by group")
      selected_df <- x[,2:length(x)]
      analyzed <- haplofunct1(selected_df)
      output <- analyzed %>% arrange(analyzed[,1], desc(analyzed$Frequency))
      first_matching_row <- x[1,2:length(x)] # for ref seq later
    }
  }
  if (ref_seq_as_H1){
    matching_row_index <- which(apply(output[, 1:(length(output)-1)], 1,
                                      function(row) all(row == first_matching_row)))
    output <- rbind(output[matching_row_index, ], 
                    output[-matching_row_index, ])
  }
  output$haplo_label <- paste("H", 1:nrow(output), sep="")
  return(output)
}

result <- haplocount(df, ref_seq_as_H1=F, group_exist = F, by_group = F)
head(result)
write.csv(result, "outputs/unique_haplotype_only_count.csv", row.names=F)

#Add H label into sorted_data
column_join_name <- list(colnames(sorted_data)[2:length(sorted_data)])
column_join_name

merged <- left_join(sorted_data, result)
merged <- merged[, -which(names(merged) == "Frequency")]
merged
write.csv(merged, "outputs/individual_haplotype_info.csv", row.names=F)
