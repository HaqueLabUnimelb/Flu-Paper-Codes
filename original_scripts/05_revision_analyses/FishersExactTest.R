# Fisher's Exact Test to determine if the degree of overlap in tow genelists is significnat or occurs by chance

# Load background gene list of the spatial data - all genes detected

total_gene_universe <- read.table(file = '/data/projects/punim2394/lung_curio_mouse/mousespatialde/backgroundgenelist.txt', stringsAsFactors = FALSE)[,1]
gene_list_1 <- read.table(file = '/data/projects/punim2394/lung_curio_mouse/mousespatialde/genelist_mouse4_spatialde.txt', stringsAsFactors = FALSE)[,1]
gene_list_2 <- read.table(file = '/data/projects/punim2394/lung_curio_mouse/mousespatialde/genelist_mouse5_spatialde.txt', stringsAsFactors = FALSE)[,1]


# Inputs
A <- gene_list_1
B <- gene_list_2
N <- length(total_gene_universe)

k <- length(intersect(A, B))
a <- length(A)
b <- length(B)


# Build contingency table
mat <- matrix(c(
  k,
  a - k,
  b - k,
  N - (a + b - k)
), nrow = 2)

# Fisher's test

fisher.test(mat, alternative = "greater")

