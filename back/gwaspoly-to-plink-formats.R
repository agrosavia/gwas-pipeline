#!/usr/bin/Rscript

# Convert gwaspoly genotype format to plink transposed format (tped, tfam)

options(stringsAsFactors = FALSE)
#------------------------------------------------------------------------------
## Write and format tassel structure
#------------------------------------------------------------------------------
writeTasselStruct <- function (structureFile) {
	structData = read.csv (file=structureFile, header=T)
	idNames = as.character (structData [,1])
	Trait = sapply (strsplit (idNames, split="And_",fixed=T), function (x) x[[2]][1])
	struct = structData [,-1]
	newData = cbind ("<Trait>"=Trait, struct)

	message (">>> Writing tassel structure")
	outFilename = paste0 (strsplit (structureFile, split="[.]")[[1]][1], "-tassel.tbl")
	sink (outFilename)
	cat ("<Covariate>\n")
	write.table (file="", newData, col.names=T, row.names=F, quote=F, sep="\t")
	sink()
}

#------------------------------------------------------------------------------
## Format and write tassel phenotype
#------------------------------------------------------------------------------
createPlinkPhenotype <- function (gwaspolyPhenotypeFile) {
	phenotype = read.csv (file=gwaspolyPhenotypeFile, header=T)
	idNames = as.character (phenotype [,1])
	Samples = sapply (strsplit (idNames, split="And_",fixed=T), function (x) x[[2]][1])
	BLIGHT = phenotype [,2]
	plinkPheno = cbind (FID=Samples,IID=Samples, BLIGHT= BLIGHT)

	message (">>> Writing tassel phenotype...")
	outFilename = paste0 (strsplit (gwaspolyPhenotypeFile, split="[.]")[[1]][1], "-plink.tbl")
	sink (outFilename)
	write.table (file=, plinkPheno, col.names=T, row.names=F, quote=F, sep="\t")
}

#----------------------------------------------------------
#----------------------------------------------------------
createFileMAP <- function (genotype, plinkFilename) {
	snpsIds = gsub ("solcap_snp_", "", gwaspolyGenotype [,1])
	chrm    = genotype [,2]+1
	pos     = genotype [3]

	genoMAP    = cbind (chr=chrm, iid=snpsIds, dist=0, pos=pos)
	message (">>> Writing MAP files...")
	filenameMAP = paste0 (plinkFilename, "-plink.map")
	write.table (file=filenameMAP, genoMAP, col.names=F, row.names=F, quote=F, sep="\t")
}

#----------------------------------------------------------
#----------------------------------------------------------
createFilePED <- function (genotype, plinkFilename) {
	#namesGeno  = c("fid", "iid", "pid", "mid", "sex", "phe", rep ("X", ncol (talleles)))
	alleles    = genotype [,-c(1,2,3)]
	talleles   = tabAlleles (alleles)
	#talleles   = t (alleles)
	samplesIds = gsub ("And_", "", colnames (alleles))
	genoPED    = cbind (samplesIds, samplesIds, 0,0,0,0, talleles)

	message (">>> Writing PED files...")
	filenamePED = paste0 (plinkFilename, "-plink.ped")
	write.table (file=filenamePED, genoPED, col.names=F, row.names=F, quote=F, sep="\t")
}

#----------------------------------------------------------
# Add tabs to alleels changign AG --> A	G
#----------------------------------------------------------
tabAlleles <- function (alleles) {
	alleles  [is.na (alleles)] = "00"
	ncols = ncol (alleles)
	nrows = nrow (alleles)
	alleles = t (alleles)
	tabs <- function (x) {return (sprintf("%s\t%s", substr(x,1,1),substr(x,2,2)))}
	allelesTabbed =matrix (sapply (alleles,tabs), ncol=nrows, nrow=ncols, byrow=F)
	return (allelesTabbed)
}

#----------------------------------------------------------
# Main
#----------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)

args = c("agrosavia-genotype-diplo-AGs.tbl","agrosavia-phenotype.tbl", "structure-checked.tbl")
gwaspolyGenotypeFile  = args [1]
gwaspolyPhenotypeFile = args [2]
structureFile         = args [3]

plinkFilename = strsplit (gwaspolyGenotypeFile, split="[.]")[[1]][1]

gwaspolyGenotype = read.table (file=gwaspolyGenotypeFile, header=T)
gw = gwaspolyGenotype

# Write PED
createFilePED (gwaspolyGenotype, plinkFilename)
# Write MAP
createFileMAP (gwaspolyGenotype, plinkFilename)

## Write and format tassel phenotype
createPlinkPhenotype (gwaspolyPhenotypeFile)

## Write and format tassel structure
#writeTasselStruct (structureFile) {


