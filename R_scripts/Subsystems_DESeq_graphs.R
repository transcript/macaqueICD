# Subsystems_DESeq_graphs.R
# Created 11/06/2017, by Sam Westreich
# Last updated 12/08/2017

args <- commandArgs(TRUE)

library(optparse)
option_list = list(
  make_option(c("-o", "--out"), type="character", default="combined_graph.pdf", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-L", "--level"), type="integer", default=1,
              help="level of Subsystems hierarchy for DESeq stats [default=%default]", metavar="character"),
  make_option(c("-d", "--directory"), type="character", default=NULL,
              help="working directory location", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print("USAGE: $ run_DESeq_graphs.R -I working_directory/ -O save.filename -L level (1,2,3,4)")

# check for necessary specs
if (is.null(opt$directory)) {
  print ("WARNING: No working directory specified with '-d' flag.")
  stop()
} else {
  cat ("Working directory is ", opt$directory, "\n")
  wd_location <- opt$directory
}

if (is.null(opt$out)) {
  print ("WARNING: No save name for plot specified; defaulting to 'combined_graph.pdf'.") }

library("DESeq2")
library("data.table")

setwd(wd_location) # OR
#setwd("~/Desktop/Projects/Lab Stuff/Macaque project/DIAMOND_results/Subsystems_results/reduced_files/renamed_files/")

# get list of files
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
  pattern = "experiment_*", full.names = T, recursive = FALSE)
exp_names = ""
for (name in exp_files) {
  exp_names <- c(exp_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
exp_names <- exp_names[-1]
exp_names_trimmed = ""
for (name in exp_names) {
  exp_names_trimmed <- c(exp_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
exp_names_trimmed <- exp_names_trimmed[-1]

# loading the control files in
y <- 0
for (x in control_files) {
  y <- y + 1
  if (y == 1) {
    control_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(control_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    control_table <- control_table[,-1]
    rownames(control_table) <- control_table$Level4 }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    rownames(temp_table) <- temp_table$Level4
    temp_table <- temp_table[,c(2,3)]
    control_table <- merge(temp_table, control_table, by = "Level4")  
  }
}
control_table <- control_table[,-ncol(control_table)]

# loading the experimental files in
y <- 0
for (x in exp_files) {
  y <- y + 1
  if (y == 1) {
    exp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(exp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    exp_table <- exp_table[,-1]
    rownames(exp_table) <- exp_table$Level4 }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    rownames(temp_table) <- temp_table$Level4
    temp_table <- temp_table[,c(2,3)]
    exp_table <- merge(temp_table, exp_table, by = "Level4")  
  }
}
exp_table <- exp_table[,-ncol(exp_table)]

# At this point, the whole table is read in.  Next step (for statistical comparison) is to 
# get just the level we want to compare.

# for Level 1 comparisons:
if (opt$level == 1) {
  l1_control_table <- data.table(control_table[, !names(control_table) %in% c("Level2", "Level3", "Level4")])
  l1_exp_table <- data.table(exp_table[, !names(exp_table) %in% c("Level2", "Level3", "Level4")])
}
if (opt$level == 2) {
  l1_control_table <- data.table(control_table[, !names(control_table) %in% c("Level1", "Level3", "Level4")])
  l1_exp_table <- data.table(exp_table[, !names(exp_table) %in% c("Level1", "Level3", "Level4")])
  names(l1_control_table)[names(l1_control_table) == 'Level2'] <- 'Level1'
  names(l1_exp_table)[names(l1_exp_table) == 'Level2'] <- 'Level1'
}
if (opt$level == 3) {
  l1_control_table <- data.table(control_table[, !names(control_table) %in% c("Level1", "Level2", "Level4")])
  l1_exp_table <- data.table(exp_table[, !names(exp_table) %in% c("Level1", "Level2", "Level4")])
  names(l1_control_table)[names(l1_control_table) == 'Level3'] <- 'Level1'
  names(l1_exp_table)[names(l1_exp_table) == 'Level3'] <- 'Level1'
}
if (opt$level == 4) {
  l1_control_table <- data.table(control_table[, !names(control_table) %in% c("Level1", "Level3", "Level2")])
  l1_exp_table <- data.table(exp_table[, !names(exp_table) %in% c("Level1", "Level3", "Level2")])
  names(l1_control_table)[names(l1_control_table) == 'Level4'] <- 'Level1'
  names(l1_exp_table)[names(l1_exp_table) == 'Level4'] <- 'Level1'
  l1_control_table <- data.table(as.data.frame(l1_control_table)[,c(2:(ncol(l1_control_table)),1)])
  l1_exp_table <- data.table(as.data.frame(l1_exp_table)[,c(2:(ncol(l1_exp_table)),1)])
}

# rename blank spots to (no hierarchy)
l1_control_table$Level1 <- ifelse(l1_control_table$Level1 == "", "NO HIERARCHY", 
                               as.character(l1_control_table$Level1))
l1_exp_table$Level1 <- ifelse(l1_exp_table$Level1 == "", "NO HIERARCHY", 
                           as.character(l1_exp_table$Level1))

# reducing stuff down to avoid duplicates
colnames(l1_control_table) <- c(control_names_trimmed, "Level1")
colnames(l1_exp_table) <- c(exp_names_trimmed, "Level1")
l1_control_table <- l1_control_table[, lapply(.SD, sum), by=Level1]
l1_exp_table <- l1_exp_table[, lapply(.SD, sum), by=Level1]
l1_control_table$Level1[is.na(l1_control_table$Level1)] <- "NO HIERARCHY"

l1_table <- merge(l1_control_table, l1_exp_table, by="Level1", all.x = T)
rownames(l1_table) <- l1_table$Level1
l1_names <- l1_table$Level1
l1_table$Level1 <- NULL
l1_table[is.na(l1_table)] <- 0
l1_data_table = data.frame(l1_table)
colnames(l1_data_table)=colnames(l1_table)
rownames(l1_data_table)=rownames(l1_table)

# now the DESeq stuff
completeCondition <- data.frame(condition=factor(c(rep("control", length(control_files)), 
  rep("ICD", length(exp_files)))))
dds <- DESeqDataSetFromMatrix(l1_data_table, completeCondition, ~ condition)
dds <- DESeq(dds)
baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( 
  counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )

# for normalized counts
normalized_counts <- as.data.frame( t( t(counts(dds)) / sizeFactors(dds)))
normalized_counts$Rownames <- row.names(l1_table)
write.table(normalized_counts, "macaque_subsys_l4_normalized_table.tsv", quote=F, row.names=T, col.names=T, sep="\t")

# continuing with normal results
res <- results(dds, contrast = c("condition", "ICD", "control"))
l1_results <- data.frame(res)
rownames(l1_results) <- rownames(l1_table)
rownames(baseMeanPerLvl) <- rownames(l1_table)
l1_results <- merge(l1_results, baseMeanPerLvl, by="row.names")
l1_results <- l1_results[,c(1,2,8,9,3,4,5,6,7)]
colnames(l1_results)[c(3,4)] <- c("controlMean", "experimentalMean")
l1_results <- l1_results[order(-l1_results$baseMean),]
l1_results[,c("lfcSE", "stat", "pvalue")] = NULL

write.table(l1_results, file = "Macaque_Subsys_l2_DESeq_8Nov.tab", 
  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# The following is to make a graph of these categories...
library(ggplot2)
graph_table <- l1_results[,c(1,5,6)]
graph_table$Row.names <- as.factor(graph_table$Row.names)
graph_table$Row.names <- factor(graph_table$Row.names, levels = graph_table$Row.names[order(graph_table$log2FoldChange)])
colnames(graph_table) <- c("Category", "log2FoldChange", "padj")
graph_table$color <- ifelse(graph_table$log2FoldChange>0, "green", "red")

# adjusting bar widths
widths = log(l1_results$baseMean)
graph_table$widths = widths
graph_table = graph_table[order(graph_table$log2FoldChange),]
graph_table$order = c(1:nrow(graph_table))
graph_table$cumwidths = cumsum(graph_table$widths/max(widths))

# the graph!
p <- ggplot(data=graph_table, aes(order, log2FoldChange, fill=color)) + 
  geom_bar(width=graph_table$widths/max(widths), stat="identity", position="identity") + 
  theme(legend.position = "none", axis.text.x = element_text(angle=90, vjust=1, hjust=1)) + 
  scale_x_continuous(breaks=graph_table$order, labels=graph_table$Category) +
  scale_fill_manual("legend", values=c("green"="green", "red"="red"))

# adding stars
graph_table <- graph_table[order(graph_table$log2FoldChange),]
label.df <- data.frame(Category=graph_table$Category, 
                       log2FoldChange = graph_table$log2FoldChange,
                       padj = graph_table$padj,
                       color=graph_table$color,
                       widths=graph_table$widths,
                       order=graph_table$order,
                       cumwidths=graph_table$cumwidths)
label.df <- subset(label.df, padj < 0.01)
label.df$padj = NULL

cat ("Saving Subsystems barplot as ", opt$out, " now.\n")
pdf(file = opt$out, width=9, height=12)
p + geom_text(data=label.df, label="*")
dev.off()
