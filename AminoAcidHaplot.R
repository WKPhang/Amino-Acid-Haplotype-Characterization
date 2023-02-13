data <- read.csv("file.csv")

# Sort seqeunces based on pattern on polymorphic site
sorted_data <- data[order(data[,2:30]), ]
write.csv(sorted_data, "sorted_file.csv")
