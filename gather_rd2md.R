#!/usr/bin/env Rscript
# Load required libraries
library(Rd2md)

# Check if any .Rd files were provided as command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript rd_to_markdown.R *.Rd\n")
  quit(status = 1)
}

# Process each .Rd file
for (rd_file in args) {
  if (!file.exists(rd_file)) {
    warning(paste("File not found:", rd_file))
    next
  }
  
  # Read the .Rd file
  rd_content <- read_rdfile(rd_file)
  
  # Convert to markdown
  md_content <- as_markdown(rd_content)
  newtitle=paste0('# ', sub('\\.\\w+$','', basename(rd_file)), "\n\n")
  md_content <- sub('^#', newtitle, md_content)
  
  md_content <- gsub('\\[(.*?)\\]\\(\\1\\)', '\\1', md_content, perl=TRUE)
  md_content <- gsub('\\[(.*?)\\]\\(\\1\\-class\\)', '\\1', md_content, perl=TRUE)
  md_content <- gsub('## Seealso', '### See also', md_content)
  md_content <- gsub('\\n## Author\\s+\\n+[^\\n]+\\n+', '', md_content, perl=TRUE)
  # Print the markdown content
  cat(md_content)
  
  # Add a separator between files
  cat("\n---\n\n")
}
