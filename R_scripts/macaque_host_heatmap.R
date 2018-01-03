# macaque_host_heatmap.R
# Created 11/28/2016, last updated 10/30/2017
# Run with --help flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input"), type="character", default="./",
              help="Input directory", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="host_heatmap", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("USAGE: $ run_DESeq_stats.R -I working_directory/ -O save.filename")

# check for necessary specs
if (is.null(opt$input)) {
  print ("WARNING: No working input directory specified with '-I' flag.")
  stop()
} else {  cat ("Working directory is ", opt$input, "\n")
  wd_location <- opt$input  
  setwd(wd_location)  }

if (is.null(opt$out)) {
  print ("WARNING: No save name for DESeq results specified; defaulting to 'DESeq_results.tab'.") 
  save_filename <- opt$out
} else { save_filename <- opt$out }

# import other necessary packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(genefilter)
  library(pheatmap)
})

# GET FILE NAMES
control_files <- list.files(
  pattern = "control_*", full.names = T, recursive = FALSE)
control_names = ""
for (name in control_files) {
  control_names <- c(control_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
control_names <- control_names[-1]
control_names_trimmed = ""
for (name in control_names) {
  control_names_trimmed <- c(control_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
control_names_trimmed <- control_names_trimmed[-1]

exp_files <- list.files(
  pattern = "experimental_*", full.names = T, recursive = FALSE)
exp_names = ""
for (name in exp_files) {
  exp_names <- c(exp_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
exp_names <- exp_names[-1]
exp_names_trimmed = ""
for (name in exp_names) {
  exp_names_trimmed <- c(exp_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
exp_names_trimmed <- exp_names_trimmed[-1]

# READ IN FILES
# loading the control table
y <- 0
for (x in control_files) {
  y <- y + 1
  if (y == 1) {
    control_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(control_table) = c("DELETE", x, "ID", "V3")
    control_table <- control_table[,c(3,4,2)]
    control_table <- control_table[complete.cases(control_table), ]
    control_table$merged <- paste(control_table$ID, control_table$V3)
    control_table <- control_table[,c(4,3)]
  }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "ID", "V3")
    temp_table <- temp_table[complete.cases(temp_table), ]
    temp_table$merged <- paste(temp_table$ID, temp_table$V3)
    temp_table <- temp_table[,c(5,2)]
    print(x)
    control_table <- merge(control_table, temp_table, by = "merged", all = T)  }
}
control_table[is.na(control_table)] <- 0
rownames(control_table) = control_table$merged
control_table_trimmed <- control_table[,-1]

# loading the experimental table
y <- 0
for (x in exp_files) {
  y <- y + 1
  if (y == 1) {
    exp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(exp_table) = c("DELETE", x, "ID", "V3")
    exp_table <- exp_table[,c(3,4,2)]
    exp_table <- exp_table[complete.cases(exp_table), ]
    exp_table$merged <- paste(exp_table$ID, exp_table$V3)
    exp_table <- exp_table[,c(4,3)]
  }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "ID", "V3")
    temp_table <- temp_table[complete.cases(temp_table), ]
    temp_table$merged <- paste(temp_table$ID, temp_table$V3)
    temp_table <- temp_table[,c(5,2)]
    print(x)
    exp_table <- merge(exp_table, temp_table, by = "merged", all = T)  }
}
exp_table[is.na(exp_table)] <- 0
rownames(exp_table) = exp_table$merged
exp_table_trimmed <- exp_table[,-1]

# getting the column names simplified
colnames(control_table_trimmed) = control_names_trimmed
colnames(exp_table_trimmed) = exp_names_trimmed

complete_table <- merge(control_table_trimmed, exp_table_trimmed, by=0, all = TRUE)
complete_table[is.na(complete_table)] <- 1
rownames(complete_table) <- complete_table$Row.names
complete_table <- complete_table[,-1]
completeCondition <- data.frame(condition=factor(c(
  rep(paste("control", 1:length(control_files), sep=".")), 
  rep(paste("ICD", 1:length(exp_files), sep=".")))))
completeCondition1 <- t(completeCondition)
colnames(complete_table) <- completeCondition1
completeCondition2 <- data.frame(condition=factor(c(
  rep("control", length(control_files)), 
  rep("ICD", length(exp_files)))))

# this stage is to trim off the gene IDs and get just names
complete_table2 <- complete_table
complete_table2_names <- data.frame(do.call('rbind', strsplit(as.character(row.names(complete_table)),'|',fixed=TRUE)))
complete_table2$names <- complete_table2_names$X5
complete_table2 <- data.table(complete_table2)
complete_table2 <- complete_table2[, lapply(.SD, sum), by=names]
saved_names = complete_table2$names
complete_table2 = data.frame(complete_table2)
row.names(complete_table2) = saved_names
complete_table2$names = NULL

dds2 <- DESeqDataSetFromMatrix(complete_table2, completeCondition2, ~condition)
dds2 <- DESeq(dds2)

rld=rlog(dds2, blind=F)

res2 = results(dds2)
org_results2 = data.frame(res2)

top_org_results = org_results2[org_results2$log2FoldChange > 5, ]
t = assay(rld)
tf = subset(t, rownames(t) %in% row.names(top_org_results))

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 10)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- tf - rowMeans(tf)
anno <- as.data.frame(colData(rld)[, c("condition","sizeFactor")])
anno2 <- as.data.frame(colData(rld)[,"condition"])
rownames(anno2) <- rownames(anno)

pheatmap(mat, annotation_col = anno2)

# removing "Predicted" from labels
labels = data.frame(do.call('rbind', strsplit(as.character(rownames(mat)),'PREDICTED:  ',fixed=TRUE)))
rownames(mat) = labels$do.call..rbind...strsplit.as.character.rownames.mat.....PREDICTED......

# adjusting the color range
paletteLength <- 100
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(-5), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat)/paletteLength, max(mat), length.out=floor(paletteLength/2)))

library(gplots) # for color setup

cat ("Saving heatmap as ", save_filename, ".pdf now.\n")
pdf(file = paste(save_filename, ".pdf", sep = ""), width=10, height=7)
pheatmap(mat, annotation_col = anno2, breaks=myBreaks,
         fontsize_row=16, fontsize_col = 16,
         color=bluered(100))
dev.off()