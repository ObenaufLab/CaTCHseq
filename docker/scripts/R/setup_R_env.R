con <- file("/tools/scripts/R/requirements.txt")
requirements.file <- readLines(con = con, warn = FALSE)
close(con)

# Remove commented lines (the ones starting with a '#')
packages <- requirements.file[!grepl(pattern = "^#", x = requirements.file)]
# Remove empty lines
packages <- packages[nchar(packages) > 0]
# Remove duplicates
packages <- unique(packages)

# Install standalone packages
setRepositories(ind = 1:3, addURLs = c("https://satijalab.r-universe.dev", "https://bnprks.r-universe.dev/"))
install.packages(packages[!grepl(pattern = "^BioConductor::", x = packages) & !grepl(pattern = "^Seurat", x = packages) & !grepl(pattern = "^remotes::", x = packages)],
    dependencies = TRUE,
    repos = "http://cran.rstudio.com/"
)

# Install Bioconductor packages
bc.packages <- packages[grepl(pattern = "^BioConductor::", x = packages)]
bc.packages <- sub(pattern = "^BioConductor::", replacement = "", x = bc.packages)
BiocManager::install(bc.packages)

# Install Seurat helpers
re.packages <- packages[grepl(pattern = "^remotes::", x = packages)]
re.packages <- sub(pattern = "^remotes::install_github", replacement = "", x = bc.packages)
remotes::install_github(c(re.packages))
