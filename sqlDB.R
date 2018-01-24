library(DBI)
library(RMariaDB)

#run mysqld --console

datadir <- 
con <- dbConnect(RMariaDB::MariaDB(), user = "root", dbname= "phospho")
dbListTables(con)

res <- dbSendQuery(con, "SELECT * FROM peptide_quantification")

conditions <- dbReadTable(con, "condition")
dbFetch(res)
