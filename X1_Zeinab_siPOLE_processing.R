#Written by Cuyler Luck
#Contact: cuyler.luck@ucsf.edu / cuylerluck@gmail.com or ross.okimoto@ucsf.edu

#This analysis is on unstranded RNA-seq data from X1 cells transfected with either a control siRNA or an siRNA targeting POLE.

#Sample key:
#A1_NC_X1 - Control X1
#A2_siPE_X1 - siPOLE X1
#B1_NC_X1 - Control X1
#B2_siPE_X1 - siPOLE X1

#First load packages and set working directory:
library(data.table) #1.14.8
library(dplyr) #1.1.1
library(edgeR) #3.40.2
library(ggplot2) #3.4.1
library(pheatmap) #1.0.12
library(biomaRt) #2.54.1
library(ggrepel) #0.9.3
library(tidyr) #1.3.0
library(patchwork) #1.1.2

#set working directory
setwd("/Volumes/cuyler/ucsf_okimoto_lab/zeinab_siPOLE") #this will change by user & location of data


#Now read in the data files
#Skip the first four lines, they just give info on odd-mapping statistics
A1 = fread("alignments/star/A1_NC_X1ReadsPerGene.out.tab", skip = 4)
A2 = fread("alignments/star/A2_siPE_X1ReadsPerGene.out.tab", skip = 4)
B1 = fread("alignments/star/B1_NC_X1ReadsPerGene.out.tab", skip = 4)
B2 = fread("alignments/star/B2_siPE_X1ReadsPerGene.out.tab", skip = 4)


#This is unstranded data generated using the NEBNext Ultra II RNA Library Prep kit for Illumina
#so, we can take column 2 from STAR as the output
#see https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

A1_unstrand = dplyr::select(A1, c(1,2))
A2_unstrand = dplyr::select(A2, c(1,2))
B1_unstrand = dplyr::select(B1, c(1,2))
B2_unstrand = dplyr::select(B2, c(1,2))


#Let's put names on the gene count columns that correspond to the sample so that we can merge these into one big data frame
colnames(A1_unstrand) = c("gene","X1_control_1")
colnames(A2_unstrand) = c("gene","X1_siPOLE_1")
colnames(B1_unstrand) = c("gene","X1_control_2")
colnames(B2_unstrand) = c("gene","X1_siPOLE_2")


#joining them in an order that makes more sense to group samples
master = inner_join(A1_unstrand, B1_unstrand, by = "gene")
master = inner_join(master, A2_unstrand, by = "gene")
master = inner_join(master, B2_unstrand, by = "gene")



#move ENSG ID to rownames for edgeR syntax
master = tibble::column_to_rownames(master, "gene")

#We can now pass this master data frame into an edgeR pipeline
#This is a simple single-factor experiment, so we will use the qCML method with the Exact test.

groups = c("X1_ctrl","X1_ctrl","X1_siPOLE", "X1_siPOLE")
dg = DGEList(counts = master, group = groups)
keep = filterByExpr(dg)
dg = dg[keep, , keep.lib.sizes = FALSE]
dg = calcNormFactors(dg)
dg = estimateDisp(dg)
#plotBCV(dg) #to see what the BCV looks like if desired. 

#now we can compare siPOLE vs control.
exact_X1_siPOLE_control = exactTest(dg)

#now we can pull out tables of results to use elsewhere
X1_siPOLE_control = exact_X1_siPOLE_control$table

#Let's add FDR-adjusted p-values to this
X1_siPOLE_control = dplyr::mutate(X1_siPOLE_control, fdr = p.adjust(PValue, method = "fdr"))

#let's also pull out the TMM-normalized log(cpm) values from the edgeR pipeline
#this will be useful later for looking at expression of specific genes
logcpm = as.data.frame(cpm(dg, log = T))


#and we can ask edgeR to plot an MDS plot to see how well samples cluster
pdf(file = "plots/plotMDS_siPOLE_control_samples.pdf", width = 10, height = 6)
plotMDS(dg)
dev.off()


#before we do anything more, I want to add on the gene symbols (not just ENSG IDs) for simplicity

#Translate gene names from ENSG IDs to gene symbols using biomaRt
#Using Ensembl genes for GRCh38.p14, which is one version newer than what I used to align samples, but should be OK (same overall genome)
#This will almost definitely lose some genes in the key for which multiple symbols match or for which no symbols exist.

#ensembl = useEnsembl(biomart = "ensembl")
#listDatasets(ensembl) #this shows that at the time of analysis, the 'hsapiens_gene_ensembl' dataset is for GRCh38.p14.
ensembl = useEnsembl(biomart="ensembl", dataset='hsapiens_gene_ensembl')

#listFilters(ensembl) #I will filter by ensembl_gene_id
#listAttributes(ensembl) #I will retrieve hgnc_symbol and ensembl_gene_id
#actually get the key from biomaRt
#I'm getting symbols for all genes that came out of STAR originally, using A1_unstrand
symbol_key = getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), filters="ensembl_gene_id",values=A1_unstrand$gene, mart=ensembl)
#trim out entries that have no gene symbol
symbol_key = symbol_key[symbol_key$hgnc_symbol!="",]

#make sure there are unique matches
length(symbol_key$hgnc_symbol)
length(unique(symbol_key$hgnc_symbol))
length(unique(symbol_key$ensembl_gene_id))
#there are duplicate symbols and keys, let's get rid of them
symbol_key = symbol_key[!duplicated(symbol_key$hgnc_symbol),]
symbol_key = symbol_key[!duplicated(symbol_key$ensembl_gene_id),]
#check again
length(symbol_key$hgnc_symbol)
length(unique(symbol_key$hgnc_symbol))
length(unique(symbol_key$ensembl_gene_id))
#now all lengths match, good.

#now we can add gene symbol to the results DFs
#for these I will left join -- this will save rows even if they have no symbol
X1_siPOLE_control = tibble::rownames_to_column(X1_siPOLE_control, "ensembl_gene_id")
X1_siPOLE_control_names = left_join(X1_siPOLE_control, symbol_key, by = "ensembl_gene_id")


#and lets add the names of genes in the logcpm dataframe. for this I only want genes that have symbols because I'm probably going to need symbols later
#so I will use inner_join

logcpm = tibble::rownames_to_column(logcpm, "ensembl_gene_id")
logcpm = inner_join(logcpm, symbol_key, by = "ensembl_gene_id")
logcpm = tibble::column_to_rownames(logcpm, "hgnc_symbol")
#then get rid of ENSMBL IDs (currently in column 1)
logcpm = logcpm[,-1]


#Plot some quick volcano plots from the data

X1_volcano = ggplot(data = X1_siPOLE_control_names, aes(x = logFC, y = -log10(fdr))) +
  geom_point(color = "gray")+
  geom_point(data = X1_siPOLE_control_names[-log10(X1_siPOLE_control_names$fdr) > 5 & X1_siPOLE_control_names$logFC<0,], color = "dark blue")+
  geom_point(data = X1_siPOLE_control_names[-log10(X1_siPOLE_control_names$fdr) > 5 & X1_siPOLE_control_names$logFC>0,], color = "maroon")+
  theme_classic()+
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = X1_siPOLE_control_names[X1_siPOLE_control_names$hgnc_symbol %in% c("POLE", "CDKN1A"),],
                  aes(label = hgnc_symbol), box.padding = 0.5, nudge_y = 5, nudge_x = -0.25) +
  ggtitle("X1 siPOLE vs control")


pdf(file = "plots/X1_volcano.pdf",width = 4, height = 4)
X1_volcano
dev.off()


#Write the full table to .csv files:

write.table(X1_siPOLE_control_names, file = "output_data/X1_siPOLE_control.csv", sep = ",", row.names = F)


#Zeinab also wants to look into specific p53 / senescence pathway genes

p53_genes = c("FAS","TNFRSF10B", "CCNB1", "CCNE1", "CDKN1A", "DDB2", "PIDD1", "SESN1", "SESN3")
sen_genes = c("ZFP36L1", "CCNA2", "CCNB1", "CCNE1", "CDKN1A", "PPP1CA", "PPP1CC", "SLC25A6")

#Let's make sure these are all in the data still
p53_genes_retained = p53_genes[(p53_genes %in% X1_siPOLE_control_names$hgnc_symbol)]
sen_genes_retained = sen_genes[(sen_genes %in% X1_siPOLE_control_names$hgnc_symbol)]

#Some overlap, so just give me unique ones
unique_p53_sen = unique(c(p53_genes_retained, sen_genes_retained))

#Yes, they are all there in the DE dataframes for both cell lines.

#Now let's make volcano plots like earlier, but highlighting these genes.

X1_volcano_pole_sen = ggplot(data = X1_siPOLE_control_names, aes(x = logFC, y = -log10(fdr))) +
  geom_point(color = "gray")+
  geom_point(data = X1_siPOLE_control_names[X1_siPOLE_control_names$hgnc_symbol %in% unique_p53_sen,], color = "dark blue")+
  theme_classic()+
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = X1_siPOLE_control_names[X1_siPOLE_control_names$hgnc_symbol %in% unique_p53_sen,],
                  aes(label = hgnc_symbol), box.padding = 0.5) +
  ggtitle("X1 siPOLE vs control, p53 or sen = blue")

pdf(file = "plots/X1_p53_sen_volcano.pdf",width = 6, height = 6)
X1_volcano_pole_sen
dev.off()


#Zeinab wants some log2(cpm) plots too

#make a new logcpm DF for plotting
logcpm_plotting = tibble::rownames_to_column(logcpm, "Gene")
logcpm_pivot = pivot_longer(logcpm_plotting, names_to = "sample", cols = c(2:5), values_to = "logcpm")

#add a variable about what kind of sample each is
logcpm_pivot = dplyr::mutate(logcpm_pivot, Condition = ifelse(sample %in% c("X1_control_1","X1_control_2"), "X1_Control", "X1_siPOLE"))

#set the order of the condition variable
logcpm_pivot$Condition = factor(logcpm_pivot$Condition, levels = c("X1_Control","X1_siPOLE"))


#This function will make a nice logcpm plot for any gene you pass into it.
make_logcpm_plot = function(name){
  gene = name
  ggplot(logcpm_pivot[logcpm_pivot$Gene %in% paste(gene),], aes(x = Condition, y = logcpm, fill = Condition)) + 
    stat_summary(fun = mean, geom = "bar") +
    geom_point() +
    theme_classic() +
    ggtitle(paste(gene)) +
    ylab("log2(cpm)") +
    coord_cartesian(ylim = c(0,12)) +
    theme(axis.title = element_text(size = 14), 
          axis.text = element_text(size = 12, hjust = 1, angle = 30),
          legend.position = "none",
          plot.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = c("red","pink"))
}

#Define genes to plot
genes_plot = c("CDKN1A","SESN3","POLE")
genes_plot = sort(genes_plot)

#then make a list containing all the plot objects
logcpm_genes_plot = lapply(genes_plot, make_logcpm_plot)


#and arrange the plots to be plotted nicely using patchwork
patchwork_logcpm_genes_plot=wrap_plots(logcpm_genes_plot[1:length(genes_plot)]) + plot_layout(ncol=5)

pdf(file="plots/patchwork_logcpm_selected_genes.pdf", width = 10, height = 3)
patchwork_logcpm_genes_plot
dev.off()

#and instead of recalculating p-values for differences here, we should us the FDR-adjusted p values from the edgeR pipeline.
X1_siPOLE_noNA = X1_siPOLE_control_names[!is.na(X1_siPOLE_control_names$hgnc_symbol),]
X1_siPOLE_noNA[X1_siPOLE_noNA$hgnc_symbol == "POLE",]$fdr
X1_siPOLE_noNA[X1_siPOLE_noNA$hgnc_symbol == "CDKN1A",]$fdr
X1_siPOLE_noNA[X1_siPOLE_noNA$hgnc_symbol == "SESN3",]$fdr


#writing out the aggregated gene counts, and the log2(cpm) data.

master_write = tibble::rownames_to_column(master, "ENSG_ID")
write.table(master_write, file = "output_data/gene_counts.txt", sep = "\t", quote = F, row.names = F)

logcpm_write = tibble::rownames_to_column(logcpm, "Gene")
write.table(logcpm_write, file = "output_data/log2cpm.txt", sep = "\t", quote = F, row.names = F)


