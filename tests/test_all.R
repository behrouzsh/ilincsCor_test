source("tests/create_test_lib.R")

create_LIB(lib_name="LIB_test", total_sigs=1000, seed_sig_csv="tests/LINCSCP_1000.tsv")

library(iLincsCor)
lib <- new(iLincsCor,"data/LIB_test/")
vec <- lib$read_input("tests/LINCSCP_1000.tsv")
x<-bench::mark(cor_vec <-lib$cor(vec,4))
print(x)
cat("Result: ")
cat(lib$signature_names[which(cor_vec>0.9)])
cat("\n")

