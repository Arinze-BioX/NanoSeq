## Read from stdin
library(data.table)
input = data.table(read.table(file("stdin"), header = T))

genotypetoid = function(gt, ref){
  sapply(gt, function(x){
    abs(1 - sum(unlist(strsplit(x, "")) == ref)/nchar(x))
  })
}


nm = names(input)[29]
input[,geno:=gsub("/", "", get(nm))]
input[,num_gt:=genotypetoid(geno, REF), by = 1:length(geno)]


## This file is identical for all hg versions, since it's computed per snp, which are invariant across reference genomes (at least for this set)
references = suppressWarnings(fread("/gpfs/data/petljaklab/resources/hg19/pipeline_resources/somatic_celline/qc/1k_cells_genotyping/reference_numeric_genotypes.txt"))
references$V1 = NULL

res = unlist(lapply(references, function(x)sum(abs(x - input$num_gt), na.rm = T)))

res = data.frame(res)
res$line = row.names(res)
res = data.table(res)
res = res[order(res)]
names(res)[1] = "distance"

write.table(res, file = "", sep = "\t", quote = F, row.names = F)
