data <- read.csv("file.csv")

# Sort seqeunces based on pattern on polymorphic site
sorted_data <- data[order(data[,3:10]), ]
write.csv(sorted_data, "sorted_file.csv")

# List unique haplotypes and their frequency
selected_data <- data[,2:10]
selected_data_strings <- apply(selected_data, 1, paste, collapse="")
frequency_table <- table(selected_data_strings)
sorted_frequency_table <- sort(frequency_table, decreasing = TRUE)
