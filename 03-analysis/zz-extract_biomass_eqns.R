# extracting the REF_SPECIES table from the fiapnw database, which we'll 
# use for biomass estimation for trees and snags. The whole database is 3.5gb
# and I have a local copy of it on the UC laptop, which I don't want to 
# copy into this project folder just for one table. This script 
# extracts the relevant table and copies it into this project folder, and
# only needs to be run once
# this file was downloaded from
# https://data.fs.usda.gov/research/pnw/rma/databases/FIAPNW.zip
# on 07/29/2021
fiadb_location = 
  'C://Users/danfoster/Documents/PILA/02-data/00-source/FIAPNW/fiapnw.db'

library(here)
library(RSQLite)
library(tidyverse)


# read in the relevant tables
fiadb = dbConnect(RSQLite::SQLite(),
                  fiadb_location)

REF_SPECIES = 
  dbReadTable(fiadb, 'REF_SPECIES')

dbDisconnect(fiadb)

# write out a csv
write.csv(REF_SPECIES,
          here::here('02-data',
                     '00-source',
                     'fia',
                     'REF_SPECIES.csv'),
          row.names = FALSE)
