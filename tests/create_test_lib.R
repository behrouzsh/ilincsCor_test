create_LIB <- function(lib_name="LIB_test", total_sigs=1000, seed_sig_csv="tests/LINCSCP_1000.tsv") {

sig<-read.csv(seed_sig_csv,,sep='\t')

chunk_size <- 2^27
output_dir <- paste0("data/",lib_name,"/")

# FIXME: if file exists info.dat read content and compare 

# skip if already exists
if (file.exists(paste0(output_dir,"info.dat"))) {
   cat("Library already exists. Skipping creation....\n")
   return(NULL)
}

cat(paste0("Creating ",lib_name,"\n"))

dir.create(output_dir)

src_data="generated_test.RData"
src_weight="generated_test.RData"

lincscpDiffExp<-matrix(rnorm(total_sigs*nrow(sig)),nrow=nrow(sig))
lincscpPValues<-matrix(rexp(total_sigs*nrow(sig)),nrow=nrow(sig))

lincscpDiffExp <- cbind(lincscpDiffExp,as.vector(unlist(sig['Value_LogDiffExp'])))
lincscpPValues <- cbind(lincscpPValues,as.vector(rep(1,nrow(lincscpPValues))))
# rownames(lincscpDiffExp) <- as.vector(unlist(sig['ID_geneid']))
# colnames(lincscpDiffExp) <- paste0("testsig_",seq(1,total_sigs))

# write.table(rownames(lincscpDiffExp), file=file.path(output_dir, "gene_ids.dat"),row.names=FALSE,col.names=FALSE)
# write.table(colnames(lincscpDiffExp), file=file.path(output_dir, "signature_names.dat"),row.names=FALSE,col.names=FALSE)
write.table(as.vector(unlist(sig['ID_geneid'])), file=file.path(output_dir, "gene_ids.dat"),row.names=FALSE,col.names=FALSE)
write.table(paste0("testsig_",seq(1,ncol(lincscpPValues))), file=file.path(output_dir, "signature_names.dat"),row.names=FALSE,col.names=FALSE)
mat_dim <- dim(lincscpDiffExp)
write.table(c(paste0("signatures=",mat_dim[2]),paste0("genes=",mat_dim[1]),paste0("src_data=",src_data),paste0("src_weight=",src_weight)), file=file.path(output_dir, "info.dat"), row.names=FALSE,col.names=FALSE,quote=FALSE)


# create data file 8
size <- length(as.vector(lincscpDiffExp))
chunks <- seq(1,size,chunk_size)
out_data <- file(file.path(output_dir, "data_8.dat"),"wb")
for(start in chunks) {
  writeBin(con=out_data, object=as.vector(lincscpDiffExp)[start:min(start+chunk_size-1,size)])
}
close(out_data)

# create data file 16
size <- length(as.vector(lincscpDiffExp))
chunks <- seq(1,size,chunk_size)
out_data <- file(file.path(output_dir, "data.dat"),"wb")
for(start in chunks) {
  writeBin(con=out_data, object=as.vector(lincscpDiffExp)[start:min(start+chunk_size-1,size)],size=4)
}
close(out_data)

# create weight file 8
size <- length(as.vector(lincscpPValues))
chunks <- seq(1,size,chunk_size)
out_data <- file(file.path(output_dir, "weight_8.dat"),"wb")
for(start in chunks) {
  writeBin(con=out_data, object=as.vector(lincscpPValues)[start:min(start+chunk_size-1,size)])
}
close(out_data)

# create weight file 16
size <- length(as.vector(lincscpPValues))
chunks <- seq(1,size,chunk_size)
out_data <- file(file.path(output_dir, "weight.dat"),"wb")
for(start in chunks) {
  writeBin(con=out_data, object=as.vector(lincscpPValues)[start:min(start+chunk_size-1,size)],size=4)
}
close(out_data)

lincscpDiffExp <- NULL
lincscpPValues <- NULL
}
