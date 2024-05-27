setwd("/Users/mortezaabyadeh/Desktop")
install.packages("readxl")
library(readxl)

excel_file <- "Sup.xlsx"
sheet_names <- excel_sheets(excel_file)

# Create an empty list to store datasets
sheet_datasets <- list()

# Loop through each sheet
for (sheet_name in sheet_names) {
  # Read data from the sheet
  sheet_data <- read_excel(excel_file, sheet = sheet_name)
  
  # Store the dataset in the list
  sheet_datasets[[sheet_name]] <- sheet_data
}


# Now, sheet_datasets contains separate datasets for each sheet
# You can access them using sheet_datasets$sheet_name
# For example: sheet_datasets$Sheet1, sheet_datasets$Sheet2, etc.

dim(sheet_datasets$`Set-1`)
head(sheet_datasets$`Set-1`, 10)



separate <- function(sheet_datasets, start_sheet = 2, end_sheet = 8) {
  result_list <- list()
  
  for (i in start_sheet:end_sheet) {
    sheet_name <- names(sheet_datasets)[i]
    
    # Extract data from the specified columns (1:3)
    data_x <- sheet_datasets[[sheet_name]][, 1:3]
    
    # Filter rows where the second column is less than 0.05
    # data_x <- data_x[data_x[, 2] < 0.05, ]
    
    # Assign names to the columns
    names(data_x) <- c("Name", "P-value", "Log2fc")
    
    # Store the result in a new dataset (e.g., Set1, Set2, ...)
    assign(paste0("Set", sheet_name), data_x)
    
    # Store the result in a list
    result_list[[sheet_name]] <- data_x
  }
  
  return(result_list)
}

# Call the function with your sheet_datasets
result <- separate(sheet_datasets, start_sheet = 2, end_sheet = 8)


dim(result$`Set-1`)
######
# Assuming you've already run the separate function and have 'result' with the seven lists
# If not, run the separate function first


###Included total number of proteins by inserting data in the return(c(nrow(data, ..)))
# Function to count rows based on Log2fc sign
count_rows_by_log2fc <- function(data, threshold = 0, p_value_threshold = 0.05) {
  # Subset the data based on both Log2fc and P-value conditions
  subset_data <- data[data$`Log2fc` < threshold & data$`P-value` < p_value_threshold, ]
  
  # Count the number of rows for both conditions
  num_rows_less_than_threshold <- nrow(subset_data)
  
  # Repeat the process for the greater than threshold condition
  subset_data_1 <- data[data$`Log2fc` > threshold & data$`P-value` < p_value_threshold, ]
  num_rows_greater_than_threshold <- nrow(subset_data_1)
  
  return(c(nrow(data), less_than_threshold = num_rows_less_than_threshold,
           greater_than_threshold = num_rows_greater_than_threshold))
}



# Get the number of rows for each Log2fc condition in each list
num_rows <- sapply(result, count_rows_by_log2fc)
dim(result$`Set-1`)

head(num_rows)
# Create a stacked bar plot
barplot(num_rows, beside = TRUE, col = c("grey", "skyblue", "lightcoral"),
        legend.text = c("Total number of proteins", "Lower abundant", "Higher abundant"),
        #args.legend = list(title = "DAPs"),
        main = "Number of lower and higher abundant proteins in each dataset",
        xlab = "Dataset", ylab = "Number of Proteins")





  

######################Number of DEG in EVs

EV <- data.frame(
  Total_number_of_proteins=c(2919, 2919, 2919, 2919, 2919, 2919), Lower_abundant=c(368, 401, 741, 392,1103,558), Higher_abundant=c(741, 1103, 368, 558, 401, 392)
)
dim(EV)
head(EV)
EV_1 <- t(EV)
head(EV_1)
colnames(EV_1) <- c("HS vs. SU", "HS vs. UC", "SU vs. HS", "SU vs. UC", "UC vs. HS", "UC vs. SU")

head(EV_1)

barplot(EV_1, beside = TRUE, col = c("grey", "skyblue", "lightcoral"),
        legend.text = c("Total number of proteins", "Lower abundant", "Higher abundant"),
        #args.legend = list(title = "DAPs"),
        main = "Number of lower and higher abundant proteins in different EV subpopulations",
        xlab = "Dataset", ylab = "Number of Proteins")






####################################################################Fisher test

###import your datasets; #####Please find the the end of codes###I have changed this function to work for dataframe


ad_up <- read.table("Up-AD.txt", header = T, stringsAsFactors = FALSE)
ad_down <- read.table("Down-AD.txt", header = T, stringsAsFactors = FALSE)


HS_SU_up <- read.table("HS-SU-Up.txt", header = T, stringsAsFactors = FALSE)
HS_SU_down <- read.table("HS-SU-Down.txt", header = T, stringsAsFactors = FALSE)


HS_UC_up <- read.table("HS-UC-Up.txt", header = T, stringsAsFactors = FALSE)
HS_UC_down <- read.table("HS-UC-Down.txt", header = T, stringsAsFactors = FALSE)


SU_HS_up <- read.table("SU-HS-Up.txt", header = T, stringsAsFactors = FALSE)
SU_HS_down <- read.table("SU-HS-Down.txt", header = T, stringsAsFactors = FALSE)

SU_UC_up <- read.table("SU-UC-Up.txt", header = T, stringsAsFactors = FALSE)
SU_UC_down <- read.table("SU-UC-Down.txt", header = T, stringsAsFactors = FALSE)


UC_HS_up <- read.table("UC-HS-Up.txt", header = T, stringsAsFactors = FALSE)
UC_HS_down <- read.table("UC-HS-Down.txt", header = T, stringsAsFactors = FALSE)

UC_SU_up <- read.table("UC-SU-Up.txt", header = T, stringsAsFactors = FALSE)
UC_SU_down <- read.table("UC-SU-Down.txt", header = T, stringsAsFactors = FALSE)


#########################################################################function for fisher's exact test
perform_fisher_test <- function(group_up, group_down, ad_up, ad_down) {
  # Get common genes within the group
  common_genes <- Reduce(intersect, list(group_up, group_down))
  
  # Calculate total_genes for all four groups
  total_genes <- length(unique(c(group_up, group_down, ad_up, ad_down)))
  
  # Contingency table
  contingency_table <- matrix(c(length(common_genes),
                                length(group_up) - length(common_genes),
                                length(ad_up) + length(ad_down) - length(common_genes),
                                total_genes - length(common_genes)),
                              nrow = 2, byrow = TRUE)
  
  # Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Print the result
  cat("Fisher's Exact Test Result:\n")
  print(fisher_test_result)
}


# Usage example for HS-SU group
perform_fisher_test(HS_SU_up, HS_SU_down, ad_up, ad_down)

###############################################################for this data fisher exact test did not work because of small sample size

#####fisher for single dataset

####import data
ad_up <- read.table("Up-AD.txt", header = T, stringsAsFactors = FALSE)$V1
ad_down <- read.table("Down-AD.txt", header = T, stringsAsFactors = FALSE)$V1


ev_up <- read.table("HS-SU-Up.txt", header = T, stringsAsFactors = FALSE)$V1
ev_down <- read.table("HS-SU-Down.txt", header = T, stringsAsFactors = FALSE)$V1

###get the common genes

common_genes <- Reduce(intersect, list(ad_up, ad_down, HS_SU_up, HS_SU_down))
total_genes <- 1715
contingency_table <- matrix(c(length(common_genes),
                              length(ad_up) + length(ad_down) - length(common_genes),
                              length(HS_SU_up) + length(HS_SU_down) - length(common_genes),
                              total_genes - length(common_genes)),
                            nrow = 2, byrow = TRUE)
###fisher test

fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)

### I did only for up_you don't need to...!

common_genes_up <- Reduce(intersect, list(ad_down, ev_up))
total_genes_up <- 1000
contingency_table <- matrix(c(length(common_genes_up),
                              length(ad_down) - length(common_genes_up),
                              length(ev_up) - length(common_genes_up),
                              total_genes_up - length(common_genes_up)),
                            nrow = 2, byrow = TRUE)
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)




###############################################################for this data fisher exact test did not work because of small sample size


###############I used hypergeometric test
ad_up <- read.delim("Up-AD.txt", header = T)
ad_down <- read.delim("Down-AD.txt", header = T)


HS_SU_up <- read.delim("HS-SU-Up.txt", header = T)
HS_SU_down <- read.delim("HS-SU-Down.txt", header = T)


HS_UC_up <- read.delim("HS-UC-Up.txt", header = T)
HS_UC_down <- read.delim("HS-UC-Down.txt", header = T)


SU_HS_up <- read.delim("SU-HS-Up.txt", header = T)
SU_HS_down <- read.delim("SU-HS-Down.txt", header = T)

SU_UC_up <- read.delim("SU-UC-Up.txt", header = T)
SU_UC_down <- read.delim("SU-UC-Down.txt", header = T)


UC_HS_up <- read.delim("UC-HS-Up.txt", header = T)
UC_HS_down <- read.delim("UC-HS-Down.txt", header = T)

UC_SU_up <- read.delim("UC-SU-Up.txt", header = T)
UC_SU_down <- read.delim("UC-SU-Down.txt", header = T)



#####Do not use the following code for data.frame, instead use the function I have provided below
# Calculate total number of unique proteins in both groups
total_proteins_all <- length(unique(c(group_up, group_down, ad_up, ad_down)))

# Calculate the number of common proteins
common_proteins <- length(intersect(c(group_up, group_down), c(ad_up, ad_down)))

# Calculate the number of proteins in each group
num_proteins_group <- length(unique(c(group_up, group_down)))
num_proteins_ad <- length(unique(c(ad_up, ad_down)))

# Perform hypergeometric test
hypergeo_test_result <- phyper(common_proteins - 1, num_proteins_group, total_proteins_all - num_proteins_group, num_proteins_ad, lower.tail = FALSE)

# Print the result
cat("Hypergeometric Test Result:\n")
print(hypergeo_test_result)

################################################ use the below function######I made

perform_hypergeometric_test <- function(group_up, group_down, ad_up, ad_down) {
  # Calculate total number of unique proteins in both groups
  total <- nrow(rbind(group_up, group_down, ad_up, ad_down))
  total_proteins_all <- nrow(unique(rbind(group_up, group_down, ad_up, ad_down)))
  common_proteins <- total - total_proteins_all
  num_proteins_group <- nrow(unique(rbind(group_up, group_down)))
  num_proteins_ad <- nrow(unique(rbind(ad_up, ad_down)))
  
  # Perform hypergeometric test
  hypergeo_test_result <- phyper(common_proteins - 1, common_proteins, total_proteins_all - common_proteins, num_proteins_ad, lower.tail = FALSE)
  
  # Return the p-value
  return(hypergeo_test_result)
}

# Usage example
p_value <- perform_hypergeometric_test(UC_SU_up, UC_SU_down, ad_up, ad_down)
cat("Hypergeometric Test Result:\n", "p-value =", p_value, "\n")

###You can make a loop through your data for the above function

p_value <- perform_hypergeometric_test(HS_SU_up, HS_SU_down, ad_up, ad_down)
cat("Hypergeometric Test Result:\n", "p-value =", p_value, "\n")
p_value

dataset_names <- c("HS_SU", "HS_UC", "SU_HS", "SU_UC", "UC_HS", "UC_SU")
file_names <- c("Up-AD.txt", "Down-AD.txt", "HS-SU-Up.txt", "HS-SU-Down.txt",
                "HS-UC-Up.txt", "HS-UC-Down.txt", "SU-HS-Up.txt", "SU-HS-Down.txt",
                "SU-UC-Up.txt", "SU-UC-Down.txt", "UC-HS-Up.txt", "UC-HS-Down.txt",
                "UC-SU-Up.txt", "UC-SU-Down.txt")


######Below loop is not showing the last two one_reason:i number that loop iterate through it
for (i in seq(1, length(dataset_names)*2, by = 2)) {
  # Read dataset files
  dataset_up <- read.delim(file_names[i], header = TRUE)
  dataset_down <- read.delim(file_names[i + 1], header = TRUE)
  
  # Perform hypergeometric test
  p_value <- perform_hypergeometric_test(dataset_up, dataset_down, ad_up, ad_down)
  
  # Print the result
  cat("Hypergeometric Test Result for", dataset_names[i/2], ":\n", "p-value =", p_value, "\n\n")
}

##########use the below function

dataset_names <- c("HS_SU", "HS_UC", "SU_HS", "SU_UC", "UC_HS", "UC_SU")
file_names <- c("Up-AD.txt", "Down-AD.txt", "HS-SU-Up.txt", "HS-SU-Down.txt",
                "HS-UC-Up.txt", "HS-UC-Down.txt", "SU-HS-Up.txt", "SU-HS-Down.txt",
                "SU-UC-Up.txt", "SU-UC-Down.txt", "UC-HS-Up.txt", "UC-HS-Down.txt",
                "UC-SU-Up.txt", "UC-SU-Down.txt")



for (i in seq(3, (length(dataset_names)+1)*2, by = 2)) {
  # Read dataset files
  dataset_up <- read.delim(file_names[i], header = TRUE)
  dataset_down <- read.delim(file_names[i + 1], header = TRUE)
  
  # Perform hypergeometric test
  p_value <- perform_hypergeometric_test(dataset_up, dataset_down, ad_up, ad_down)
  
  # Print the result
  cat("Hypergeometric Test Result for", dataset_names[i/2], ":\n", "p-value =", p_value, "\n\n")
}


################################################################################corrceted fisher's exact test for these datasets

perform_fisher_test <- function(group_up, group_down, ad_up, ad_down) {
  # Get common genes within the group
  total <- nrow(rbind(group_up, group_down, ad_up, ad_down))
  total_proteins_all <- nrow(unique(rbind(group_up, group_down, ad_up, ad_down)))
  common_proteins <- total - total_proteins_all

  # Contingency table
  contingency_table <- matrix(c(length(common_proteins),
                                length(group_up) - length(common_proteins),
                                length(ad_up) + length(ad_down) - length(common_proteins),
                                total_proteins_all - length(common_proteins)),
                              nrow = 2, byrow = TRUE)
  
  # Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Print the result
  cat("Fisher's Exact Test Result:\n")
  print(fisher_test_result)
}


# Usage example for HS-SU group
perform_fisher_test(HS_SU_up, HS_SU_down, ad_up, ad_down)

####somehow the perform_fisher_test function is not working in the loop##########use the next loop
for (i in seq(3, (length(dataset_names)+1)*2, by = 2)) {
  dataset_up <- read.delim(file_names[i], header = TRUE)
  dataset_down <- read.delim(file_names[i + 1], header = TRUE)
  
  p_value <- perform_fisher_test(dataset_up, dataset_down, ad_up, ad_down)
  
  cat("Fisher Test Result for", dataset_names[i/2], ":\n", "p-value =", p_value, "\n\n")
}

#################following code used the function and used the fisher results from fisher.test(contingency_table), means used function in the middle and directly got the P-value from function
for (i in seq(3, (length(dataset_names)+1)*2, by = 2)) {
  # Read dataset files
  dataset_up <- read.delim(file_names[i], header = TRUE)
  dataset_down <- read.delim(file_names[i + 1], header = TRUE)
  
  # Perform Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Extract the p-value from the result
  p_value <- fisher_test_result$p.value
  
  # Print the result
  cat("Fisher Test Result for", dataset_names[i/2], ":\n", "p-value =", p_value, "\n\n")
}

############################## fisher test #############use the next function

perform_fisher_test <- function(group_up, group_down, ad_up, ad_down) {
  # Get common genes within the group
  total <- nrow(rbind(group_up, group_down, ad_up, ad_down))
  total_proteins_all <- nrow(unique(rbind(group_up, group_down, ad_up, ad_down)))
  common_proteins <- total - total_proteins_all
  
  # Contingency table
  contingency_table <- matrix(c(length(common_proteins),
                                length(group_up) + length(group_down) - length(common_proteins),
                                length(ad_up) + length(ad_down) - length(common_proteins),
                                total_proteins_all - length(common_proteins)),
                              nrow = 2, byrow = TRUE)
  
  # Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  return(fisher_test_result)
}

perform_fisher_test(HS_SU_up, HS_SU_down, ad_up, ad_down)
perform_fisher_test(HS_UC_up, HS_UC_down, ad_up, ad_down)
perform_fisher_test(SU_UC_up, SU_UC_down, ad_up, ad_down)


##################USED above individual fisher's exact test#####################


########Its returning same P-value!
################### Use the followign loop as included both up and down regulated proteins, others only had up regulated.
dataset_names <- c("HS_SU", "HS_UC", "SU_HS", "SU_UC", "UC_HS", "UC_SU")
file_names <- c("Up-AD.txt", "Down-AD.txt", "HS-SU-Up.txt", "HS-SU-Down.txt",
                "HS-UC-Up.txt", "HS-UC-Down.txt", "SU-HS-Up.txt", "SU-HS-Down.txt",
                "SU-UC-Up.txt", "SU-UC-Down.txt", "UC-HS-Up.txt", "UC-HS-Down.txt",
                "UC-SU-Up.txt", "UC-SU-Down.txt")


for (i in seq(3, (length(dataset_names)+1)*2, by = 2)) {
  dataset_up <- read.delim(file_names[i], header = TRUE)
  dataset_down <- read.delim(file_names[i + 1], header = TRUE)
  
  total <- nrow(rbind(file_names[i], file_names[i+1], ad_up, ad_down))
  total_proteins_all <- nrow(unique(rbind(file_names[i], file_names[i+1], ad_up, ad_down)))
  common_proteins <- total - total_proteins_all
  
  # Contingency table
  contingency_table <- matrix(c(length(common_proteins),
                                length(file_names[i]) + length(file_names[i+1]) - length(common_proteins),
                                length(ad_up) + length(ad_down) - length(common_proteins),
                                total_proteins_all - length(common_proteins)),
                              nrow = 2, byrow = TRUE)
  
  # Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  p_value <- fisher_test_result$p.value
  
  # Print the result
  cat("Fisher Test Result for", dataset_names[i/2], ":\n", "p-value =", p_value, "\n\n")
}



##########Need to work on it. ############It's showing negative value!

dataset_names <- c("HS_SU", "HS_UC", "SU_HS", "SU_UC", "UC_HS", "UC_SU")
file_names <- c("Up-AD.txt", "Down-AD.txt", "HS-SU-Up.txt", "HS-SU-Down.txt",
                "HS-UC-Up.txt", "HS-UC-Down.txt", "SU-HS-Up.txt", "SU-HS-Down.txt",
                "SU-UC-Up.txt", "SU-UC-Down.txt", "UC-HS-Up.txt", "UC-HS-Down.txt",
                "UC-SU-Up.txt", "UC-SU-Down.txt")

for (i in seq(3, (length(dataset_names)+1)*2, by = 2)) {
  dataset_up <- read.delim(file_names[i], header = TRUE)
  dataset_down <- read.delim(file_names[i + 1], header = TRUE)
  
  total <- nrow(rbind(dataset_up, dataset_down, ad_up, ad_down))
  total_proteins_all <- nrow(unique(rbind(dataset_up, dataset_down, ad_up, ad_down)))
  common_proteins <- total - total_proteins_all
  
  # Contingency table
  contingency_table <- matrix(c(common_proteins,
                                nrow(dataset_up) + nrow(dataset_down) - common_proteins,
                                length(ad_up) + length(ad_down) - common_proteins,
                                total_proteins_all - common_proteins),
                              nrow = 2, byrow = TRUE)
  
  # Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  p_value <- fisher_test_result$p.value
  
  # Print the result
  cat("Fisher Test Result for", dataset_names[i/2], ":\n", "p-value =", p_value, "\n\n")
}



###################Modification through most of above codes#####################

total <- nrow(rbind(group_up, group_down, ad_up, ad_down))
total_proteins_all <- nrow(unique(rbind(group_up, group_down, ad_up, ad_down)))
common_proteins <- total - total_proteins_all
num_proteins_group <- nrow(unique(rbind(group_up, group_down)))
num_proteins_ad <- nrow(unique(rbind(ad_up, ad_down)))



####################################################################Data mining

install.packages("easyPubMed")
library(easyPubMed)
myQuery <- "POU Class 3 Homeobox 2 AND Alzheimer"
myIdList <- get_pubmed_ids(myQuery)
all_steps <- seq(1, myIdList$Count, by = 50)
results <- lapply(all_steps, function(i) {
  y <- fetch_pubmed_data(pubmed_id_list = myIdList, retmax = 50, retstart = i)  
  yy <- table_articles_byAuth(y, included_authors = "first", getKeywords = TRUE)
  yy[, c("title", "pmid", "doi", "jabbrv", "lastname", "keywords")]
})
results <- do.call(rbind,results)
nrow(results)
head(results)
library(openxlsx)
write.xlsx(results, "/Users/mortezaabyadeh/Documents/AD proteome/EV-AD/Data mining/POU3F2.xlsx")



######################################### Fisher for individuals not a function ############################################

perform_fisher_test <- function(group_up, group_down, ad_up, ad_down) {
  # Get common genes within the group
  total <- nrow(rbind(group_up, group_down, ad_up, ad_down))
  total_proteins_all <- nrow(unique(rbind(group_up, group_down, ad_up, ad_down)))
  common_proteins <- total - total_proteins_all
  
  # Contingency table
  contingency_table <- matrix(c(length(common_proteins),
                                length(group_up) + length(group_down) - length(common_proteins),
                                length(ad_up) + length(ad_down) - length(common_proteins),
                                total_proteins_all - length(common_proteins)),
                              nrow = 2, byrow = TRUE)
  
  # Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  return(fisher_test_result)
}

perform_fisher_test(HS_SU_up, HS_SU_down, ad_up, ad_down)
perform_fisher_test(HS_UC_up, HS_UC_down, ad_up, ad_down)
perform_fisher_test(SU_UC_up, SU_UC_down, ad_up, ad_down)

####################### Vikas Data ######################

perform_fisher_test <- function(Data_1, Data_2, Data_4) {
  # Get common genes within the group
  total_proteins_all <- Data_1 + Data_2
  common_proteins <- Data_4
  
  # Contingency table
  contingency_table <- matrix(c(length(common_proteins),
                                length(Data_1) - length(common_proteins),
                                length(Data_2) - length(common_proteins),
                                total_proteins_all - length(common_proteins)),
                              nrow = 2, byrow = TRUE)
  
  # Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  return(fisher_test_result)
}

perform_fisher_test(234, 456, 23)
