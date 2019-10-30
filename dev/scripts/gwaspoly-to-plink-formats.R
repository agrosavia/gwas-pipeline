#!/usr/bin/Rscript

# Convert gwaspoly genotype format to plink transposed format (tped, tfam)

args = commandArgs(trailingOnly = TRUE)

args = c("genotype-diplo.tbl")
gwaspolyGenotypeFile = args [1]

gwaspolyGenotype = read.csv (file=gwaspolyGenotypeFile, header=T)
gw = gwaspolyGenotype

iid     = gwaspolyGenotype [,1]
fid     = iid
chrm    = gwaspolyGenotype [,2]
dist    = 0
pos     = gwaspolyGenotype [3]
alleles = gwaspolyGenotype [,-c(1,2,3)]

alleles [alleles==0] = "A A"
alleles [alleles==1] = "A B"
alleles [alleles==2] = "B B"

plinkGenotypeTPED = cbind (iid,fid,chrm,dist,pos,alleles)
pl = plinkGenotype

plinkTPEDFilename = paste0 (strsplit (gwaspolyGenotypeFile, split="[.]")[[1]][1], "-plink.tped")
write.table (file=plinkTPEDFilename, plinkGenotype, col.names=F, row.names=F, quote=F)

snpsNames = colnames (gwaspolyGenotype)[-c(1,2,3)]
iid =fid = snpsNames
plinkFamilyTFAM = cbind (iid,fid,0,0,0,-9)
pf = plinkFamilyTFAM

plinkTFAMFilename = paste0 (strsplit (gwaspolyGenotypeFile, split="[.]")[[1]][1], "-plink.tfam")
write.table (file=plinkTFAMFilename, plinkFamilyTFAM, col.names=F, row.names=F, quote=F)

