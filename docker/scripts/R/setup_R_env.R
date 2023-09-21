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
install.packages(packages[!grepl(pattern = "^BioConductor::", x = packages)], 
                 dependencies = TRUE, 
                 repos = 'http://cran.rstudio.com/')

# Install Bioconductor packages
bc.packages <- packages[grepl(pattern = "^BioConductor::", x = packages)]
bc.packages <- sub(pattern = "^BioConductor::", replacement = "", x = bc.packages)
BiocManager::install(bc.packages)