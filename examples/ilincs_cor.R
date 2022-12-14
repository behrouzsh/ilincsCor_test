library(httr)
setwd("/Users/behrouzshamsaei/Documents/ilincs_cor/")
urlupload <- 'https://shiny.ilincs.org/devcor/api/SignatureMeta/upload'
r <- POST(urlupload, body = list(file = upload_file("tests/LINCSCP_100.tsv")))
filename <- content(r)[[1]]$fileName[[1]]
filename
url <- "https://shiny.ilincs.org/devcor/ilincs_cor"
lib <- "LIB_5"
limit <- 0.3
workers <- 20
arg <- paste0("file=", filename,"&lib=",lib,"&limit=",limit,"&workers=", workers )
upload_response <- POST(url, body = arg, encode = "form")
#upload_response <- POST(url, body = "file=LINCSCP_100.tsv&lib=LIB_5&limit=0.3&workers=16", encode = "form")
content(upload_response)