setwd("\CPTAC analysis 06-04-2025")
#################################################### Normal Weight###########################################################################################


library(TCGAbiolinks)

#barcodes
barcodes <- c(
  "CPT0209640011", "CPT0065570006", "CPT0139000007", "CPT0247090013", "CPT0238040011",
  "CPT0078090010", "CPT0124410010", "CPT0124360009", "CPT0162660010", "CPT0208690010",
  "CPT0123990006", "CPT0093860006", "CPT0064520007", "CPT0094130007", "CPT0094030006",
  "CPT0236700011", "CPT0208840010", "CPT0161990011", "CPT0161770010", "CPT0081470010",
  "CPT0123750009", "CPT0126160007", "CPT0088180007", "CPT0109110007", "CPT0019300007",
  "CPT0124980010", "CPT0093930007", "CPT0124500010", "CPT0183750011", "CPT0011520009",
  "CPT0248590014", "CPT0247040006", "CPT0246210007", "CPT0226570010", "CPT0238820011",
  "CPT0237910010", "CPT0237840010", "CPT0197420010", "CPT0207780007", "CPT0264150006",
  "CPT0246840013", "CPT0254680009", "CPT0238710010", "CPT0209310014", "CPT0276010009",
  "CPT0166650012", "CPT0226680011", "CPT0198700011", "CPT0236560010", "CPT0183450011",
  "CPT0236780011"
)

# Query 
normal_query <- GDCquery(
  project = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = barcodes
)

# Download the data
GDCdownload(normal_query)
normal_data <- GDCprepare(normal_query)


library(SummarizedExperiment)


# Create results directory
dir.create("results", showWarnings = FALSE)
dir.create("results/expression_matrix", showWarnings = FALSE)

# Convert the expression matrix to CSV format
write.csv(assay(normal_data, "unstranded"), file = "results/expression_matrix/pdac_normal_weight.csv")

# Print the first few rows of the data to confirm
print(head(assay(normal_data, "unstranded")))

# View the first few rows of the data
print(head(normal_data))

#################################################### Overweight###########################################################################################

# Define the barcodes
overweight_barcodes <- c(
  "CPT0063970010", "CPT0008750009", "CPT0109300007", "CPT0094940007", "CPT0078000009",
  "CPT0078250009", "CPT0125090009", "CPT0123440009", "CPT0239100015", "CPT0246720009",
  "CPT0227930010", "CPT0092630006", "CPT0109200007", "CPT0088260007", "CPT0001730008",
  "CPT0083270009", "CPT0087890007", "CPT0077860009", "CPT0025750006", "CPT0226760010",
  "CPT0238610015", "CPT0226490010", "CPT0127270009", "CPT0065640006", "CPT0008800009",
  "CPT0108400006", "CPT0218200010", "CPT0221470011", "CPT0208550011", "CPT0124610010",
  "CPT0123650006", "CPT0078160010", "CPT0086190007", "CPT0064090009", "CPT0108930006",
  "CPT0081260010", "CPT0000250008", "CPT0125130009", "CPT0198420010", "CPT0094740007",
  "CPT0081360010", "CPT0017040010", "CPT0089720006", "CPT0246810013", "CPT0238440011",
  "CPT0185920010", "CPT0238920010", "CPT0237270008", "CPT0238530011", "CPT0187550010",
  "CPT0226600010", "CPT0187930011", "CPT0241320011", "CPT0239910010", "CPT0207980011",
  "CPT0224460010", "CPT0170240010", "CPT0217700010"
)

# Query 
overweight_query <- GDCquery(
  project = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = overweight_barcodes
)

# Download the data
GDCdownload(overweight_query)
overweight_data <- GDCprepare(overweight_query)


library(SummarizedExperiment)

# Convert the expression matrix to CSV format
write.csv(assay(overweight_data, "unstranded"), file = "results/expression_matrix/pdac_over_weight.csv")

# Print the first few rows of the data to confirm
print(head(assay(overweight_data, "unstranded")))

# View the first few rows of the data
print(head(overweight_data))

#################################################### Obese ###########################################################################################
#barcodes
Obese_barcodes <- c(
  "CPT0077930010", "CPT0109010007", "CPT0094230007", "CPT0186040010", "CPT0166760010",
  "CPT0091660006", "CPT0063570006", "CPT0123860006", "CPT0014040009", "CPT0094180006",
  "CPT0174700010", "CPT0081060007", "CPT0094630007", "CPT0162760009", "CPT0246940009",
  "CPT0198910011", "CPT0248410010", "CPT0241500010"
)

# Query
Obese_query <- GDCquery(
  project = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = Obese_barcodes
)

# Download the data
GDCdownload(Obese_query)
Obese_data <- GDCprepare(Obese_query)

# Convert the expression matrix to CSV format
write.csv(assay(Obese_data, "unstranded"), file = "results/expression_matrix/pdac_Obese.csv")

# Print the first few rows of the data to confirm
print(head(assay(Obese_data, "unstranded")))