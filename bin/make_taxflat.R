#!/usr/bin/env Rscript

# make_taxflat.R
#
# Reads nodes.dmp and names.dmp from NCBI's taxdump.tar.gz and reformats to a
# flat file with ncbi_taxon_id, rank, and recognized taxa from domain to 
# strain.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(stringr))

SCRIPT_VERSION = "0.9"

options(warn = 1)

# Get arguments
# For testing opt = list(options = list(verbose = TRUE), args = c('nodes.dmp', 'names.dmp', 'new_taxflat.tsv.gz'))
option_list = list(
  make_option(
    c("-t", "--threads"), default=1, 
    help="Set number of cpu threads"
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  ),
  make_option(
    c("-V", "--version"), action="store_true", default=FALSE, 
    help="Print program version and exit"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] nodes.dmp names.dmp output.tsv.gz",
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

if ( opt$options$version ) {
  write(SCRIPT_VERSION, stdout())
  quit('no')
}

DEBUG   = 0
INFO    = 1
WARNING = 2
LOG_LEVELS = list(
  DEBUG   = list(n = 0, msg = 'DEBUG'),
  INFO    = list(n = 1, msg = 'INFO'),
  WARNING = list(n = 2, msg = 'WARNING'),
  ERROR   = list(n = 3, msg = 'ERROR')
)
logmsg    = function(msg, llevel='INFO') {
  if ( opt$options$verbose | LOG_LEVELS[[llevel]][["n"]] >= LOG_LEVELS[["INFO"]][["n"]] ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}

setDTthreads(opt$options$threads)

logmsg(sprintf("Reading nodes from %s", opt$args[1]))
nodes <- fread(
  opt$args[1],
  col.names = c('tax_id', 'parent', 'rank', 'embl', 'division', 'inherited_div_flag', 'genetic_code_id', 'inherited_gc_flag', 'mitochondrial', 'inherited_mgc_flag', 'GenBank_hidden', 'hidden_subtree', 'comments', 'x'),
  sep = '|'
) %>% lazy_dt() %>%
  transmute(tax_id = tax_id, parent = parent, rank = str_remove_all(rank, '\t')) %>%
  as.data.table()

logmsg(sprintf("Reading names from %s", opt$args[2]))
names <- fread(
  opt$args[2],
  col.names = c('tax_id', 'name_txt', 'unique_name', 'name_class', 'x'),
  sep = '|'
) %>% lazy_dt() %>%
  transmute(tax_id, name = str_remove_all(name_txt, '\t'), name_class = str_remove_all(name_class, '\t')) %>%
  as.data.table()

# 0. Create a table joining nodes and names, and with a name column containing the rank
nodes_names <- lazy_dt(nodes) %>%
  filter(tax_id != 1) %>%       # This is just the parent of the actual top level
  inner_join(lazy_dt(names) %>% filter(name_class == 'scientific name') %>% select(-name_class), by = 'tax_id') %>%
  mutate(
    rname = case_when(
      rank == 'superkingdom' ~ sprintf("d__%s", name),
      rank == 'kingdom'      ~ sprintf("k__%s", name),
      rank == 'phylum'       ~ sprintf("p__%s", name),
      rank == 'class'        ~ sprintf("c__%s", name),
      rank == 'order'        ~ sprintf("o__%s", name),
      rank == 'family'       ~ sprintf("f__%s", name),
      rank == 'genus'        ~ sprintf("g__%s", name),
      rank == 'species'      ~ sprintf("s__%s", name),
      TRUE                   ~ ''
    )
  ) %>%
  as.data.table()

# 1. Create a table containing the lowest level, i.e with parent == 1; create a thier column
logmsg("Joining children with parents")
flat <- lazy_dt(nodes_names) %>%
  filter(parent == 1) %>%
  mutate(thier = ifelse(grepl('^[a-z]__', rname), rname, '')) %>%
  select(-rname) %>%
  as.data.table()

# Remove rows in t that are now in flat
nodes_names <- lazy_dt(nodes_names) %>% anti_join(lazy_dt(flat), by = 'tax_id') %>% as.data.table()

i = 0
# 2. Loop until we've emptied the nodes_names, and build the hierarchical name
while(nrow(nodes_names) > 0) {
  logmsg(sprintf("\t--> %3d: %7d rows remaining in nodes_names, %7d in flat <--", i, nrow(nodes_names), nrow(flat)))
  flat <- funion(
    flat,
    lazy_dt(nodes_names) %>%
      inner_join(
        lazy_dt(flat) %>% select(tax_id, thier),
        by = c('parent' = 'tax_id')
      ) %>%
      mutate(thier = sprintf("%s<-->%s", thier, rname)) %>%
      select(-rname) %>%
      as.data.table()
    )

  nodes_names <- lazy_dt(nodes_names) %>% anti_join(lazy_dt(flat), by = 'tax_id') %>% as.data.table()

  # If there are any nodes without parents, or other complications, break out if we did 100 iterations (30-40 is normal)
  i <- i + 1
  if ( i > 100 ) break
}

# Create the rank fields, and write to file
logmsg(sprintf("Writing the final table to %s", opt$args[3]))
lazy_dt(flat) %>%
  transmute(
    tax_id, rank, name,
    domain   = ifelse(grepl('d__', thier), str_replace(thier, '.*d__(.*)', '\\1') %>% str_remove('<-->.*'), NA),
    kingdom  = ifelse(grepl('k__', thier), str_replace(thier, '.*k__(.*)', '\\1') %>% str_remove('<-->.*'), NA),
    phylum   = ifelse(grepl('p__', thier), str_replace(thier, '.*p__(.*)', '\\1') %>% str_remove('<-->.*'), NA),
    class    = ifelse(grepl('c__', thier), str_replace(thier, '.*c__(.*)', '\\1') %>% str_remove('<-->.*'), NA),
    order    = ifelse(grepl('o__', thier), str_replace(thier, '.*o__(.*)', '\\1') %>% str_remove('<-->.*'), NA),
    family   = ifelse(grepl('f__', thier), str_replace(thier, '.*f__(.*)', '\\1') %>% str_remove('<-->.*'), NA),
    genus    = ifelse(grepl('g__', thier), str_replace(thier, '.*g__(.*)', '\\1') %>% str_remove('<-->.*'), NA),
    species  = ifelse(grepl('s__', thier), str_replace(thier, '.*s__(.*)', '\\1') %>% str_remove('<-->.*'), NA),
  ) %>%
  as.data.table() %>%
  fwrite(opt$args[3], sep = '\t', row.names = FALSE)

logmsg("Done")
