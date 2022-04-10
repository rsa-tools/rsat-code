#############################
## Load required libraries ##
#############################

## List of packages to install from CRAN
required.packages = c("dplyr",
                      "gt", 
                      "lubridate",
                      "optparse")

for (lib in required.packages) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}


option_list = list(
  
  make_option(c("-m", "--motif_table"), type = "character", default = NULL, 
              help = "RSAT motif DB table: $RSAT/public_html/motif_databases/db_matrix_files.tab (Mandatory). ", metavar = "character"),
  
  make_option(c("-o", "--output_folder"), type = "character", default = NULL, 
              help = "Folder to save the results (Mandatory)", metavar = "character")
);

message("; Reading arguments from command-line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


## Output file prefix
rsat.motif.db.file <- opt$motif_table
output.folder      <- opt$output_folder

dir.create(output.folder, recursive = TRUE, showWarnings = FALSE)

## Example: Rscript $RSAT/R-scripts/motif_databases/RSAT_motif_DB_table.R -m /home/rsat/packages/rsat-2021/rsat/public_html/motif_databases/db_matrix_files.tab -o $RSAT/public_html/motif_databases
# rsat.motif.db.file <- "db_matrix_files.tab"
# output.folder      <- "."


###############
## Functions ##
###############

count.motifs.in.db <- function(motif.collection.lines = NULL) {
  
  ## All motifs are in transfac, se we count the number of separators
  return(sum(grepl(motif.collection.lines, pattern = "^//\\s*$")))
}


##########
## Main ##
##########
message("; Reading RSAT DB table: ", rsat.motif.db.file)
rsat.motif.db.tab <- read.table(rsat.motif.db.file, header = TRUE, comment.char = ";", sep = "\t") %>%
  rename(Collection  = X.COLLECTION,
         Format      = FORMAT,
         File        = FILE,
         Description = DESCR,
         Version     = VERSION,
         Category    = CATEGORY,
         DataBase    = DATABASE)


## Read the tf files
rsat.motif.db.content <- purrr::map(.x = as.vector(rsat.motif.db.tab$File),
                                    .f = ~readLines(con = .x))

## Count the number of motifs per collection
## We count the number of separators '//' in each file
message("; Counting motifs in each collection")
rsat.motif.db.nb.motifs <- purrr::map_dbl(.x = rsat.motif.db.content,
                                          .f = ~count.motifs.in.db(motif.collection.lines = .x))

# save.image("TEST.Rdata")
rsat.motif.db.tab$Nb_motifs <- rsat.motif.db.nb.motifs

rsat.motif.db.tab <- rsat.motif.db.tab %>% 
                      arrange(Category, desc(Nb_motifs))



## Export html document
message("; Generating html file")
rsat.motif.db.tab %>%
  # mutate(Nb_motifs_log10 = log10(Nb_motifs)) %>% 
  select(Category, DataBase, Collection, Version, Description, Nb_motifs, URL) %>% 
  group_by(Category) %>% 
  gt(rowname_col = "DataBase") %>%
  tab_header(
    title    = md(paste0("Motif collections in RSAT by ", today())),
    subtitle = md(paste0("List of the ", nrow(rsat.motif.db.tab), " motif collections integrated in RSAT. This table lists the collection names, date of integration or version, a brief description, and the URL indicating the origin of each database."))) %>% 
  data_color(
    columns = Nb_motifs,
    colors = scales::col_numeric(palette = c("white", "yellow", "darkred"),
                                 domain = range(Nb_motifs))) %>% 
  # tab_spanner_delim(delim = ".") %>% 
  cols_align(align = "center") %>%
  tab_options(
    table.width = pct(100),
    row_group.as_column = TRUE
  ) %>% 
  grand_summary_rows(
    columns = c(Nb_motifs),
    fns = list(Total = ~sum(.)),
    formatter = fmt_number,
    decimal = 0,
    use_seps  = TRUE
  ) %>% 
  gtsave(file.path(output.folder, "RSAT_motif_DB_table.html"))


