rm(list = ls())

#### PRELIMINARIES ############################################################################################# 

#*Uploads the needed libraries --------------------------------------------------------------------------------

require(DESeq2)

require(ggplot2)

require(data.table)

require(plotly)

require(DT)

require(htmlwidgets)

require(R2HTML)

require(biomaRt)

soft_version <- packageVersion("DESeq2")

#*Set the number of significant digits for the output --------------------------
sig_dig = 4

#*Set the working directory ---------------------------------------------------------------------------------

#get the default wd
default_wd <- getwd()

#Set the directory where all the saved outputs will be stored
setwd("~/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Minoru_project/Second_batch/DE_Analysis_feature_counts/Tables_and_figures/")# <--- insert here the path to the working (output) directory
#setwd("/media/davide/Dropbox/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Nadia_Schoenmaker_project/DE_analysis_feature_counts/Tables_and_figures_HE_vs_KO/")# <--- insert here the path to the working (output) directory

new_wd <- getwd()
#### DATA UPLOAD ###############################################################################################################

annotation_table <- read.csv("~/Dropbox (Cambridge University)/Bioinformatics_core/Useful_bash_scripts_and_references/Mus_musculus.GRCm38.95_gene_annotation_table.txt", sep = "\t")
#annotation_table <- read.csv("/media/davide/Dropbox/Dropbox (Cambridge University)/Bioinformatics_core/Useful_bash_scripts_and_references/Homo_sapiens.GRCh38.95_gene_annotation_table.txt", sep = "\t")

# inputs the list of the count files
#input_files <- list.files(path = "/media/davide/Dropbox/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Pani_Gopaluni_project/DE_analysis_feature_counts/Counts/", pattern = "*CD04_act", full.names = TRUE) #<--- insert here the path to the working directory
input_files <- list.files(path = "~/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Minoru_project/Second_batch/DE_Analysis_feature_counts/Counts/", pattern = "*tab$", full.names = TRUE) #<--- insert here the path to the working directory
#input_files <- input_files[grep(input_files, pattern="b_WT|b_KO")]

# inputs the list of the summary files
summary_files <-  list.files(path = "~/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Minoru_project/Second_batch/DE_Analysis_feature_counts/Counts/", pattern = "summary$", full.names = TRUE)
#summary_files <- summary_files[grep(summary_files, pattern="b_WT|b_KO")]

#Reads the count files (input_files)
#create a list; each element of a list (named ''sample'' in the following) is a count table
counts_tables <- lapply(input_files, fread, header = FALSE, stringsAsFactors = FALSE, skip=2, select=c(1,7), col.names=c("GeneID","Counts"))

#Reads the summary files (summary_files)
#create a list; each element of a list is a summary table
summary_tables <- lapply(summary_files, read.delim, header = FALSE, stringsAsFactors = FALSE, skip=1, quote = "")

#takes the names of the genes from the counts_tables
genes_names <- counts_tables[[1]]$GeneID
genes_number = length(genes_names)

#take the names of each element (sample) of the list from the input files 
samples_names <- substr(noquote(unlist(lapply(basename(input_files), tools::file_path_sans_ext))) ,1,9)

#assign the names to the elements of the counts_tables list, composed by the counts tables; NOTE: each replicate has its ID
names(counts_tables) <- samples_names

#Creates a single data frame with all the samples as columns, for reporting and clustering purposes -- see heatmaps below
counts_tables_dataframe <- sapply(counts_tables, '[[', 2)
#assignes the genes names to the rows of the counts_tables_dataframe
rownames(counts_tables_dataframe) <- genes_names

#assign the names to the elements of the summary_tables list, composed by the summary tables; NOTE: each replicate has its ID
names(summary_tables) <- samples_names

#take the roots of the samples_names; replicates cannot be distinguished here -- set the name of the control 
name_control <- "Group_1"                                                                        #<- insert here the NAME OF THE CONTROL
names_mutants <- setdiff(unique((substr(samples_names,1,7))), name_control)

#Creates the datasets_matrix, containing all the names of the datasets and replicates for reporting purposes
datasets_matrix <- matrix(samples_names, nrow=length(grep(name_control, samples_names)), ncol=length(names_mutants)+1)
colnames(datasets_matrix) <- c(name_control, names_mutants)

####PAIRWISE COMPARISONS ####################################################################################################
# *Set the tresholds ----------------------------------------------------------------------------------------------------

lfc = 1.0 #treshold for the ADJUSTED fold changes
pval = 0.01 #treshold for the ADJUSTED pvalues

num_comparisons=1 #loop's counter initialisation

###COMPARISONS, PLOTS AND TABLES LOOP ======================================================================================

# *Start of the loop ------------------------------------------------------------------------------------------------------
while(num_comparisons <= length(names_mutants)){    

  # **Grabbing the data ----------------------------------------------------------------------------------------------------  
  name_mutant <- names_mutants[num_comparisons]  
  
  #extract the needed samples from the whole list
  assign( paste0("counts_tables_", name_control), counts_tables[grep(name_control, samples_names, value = TRUE)]) #extract the controls counts tables 

  assign( paste0("counts_tables_", name_mutant), counts_tables[grep(name_mutant, samples_names, value = TRUE)]) #extract the mutants counts tables 

  assign( paste0("counts_list_", name_control, "_", name_mutant) , c( get( paste0("counts_tables_", name_control)),  get(paste0("counts_tables_", name_mutant)) )) #merges the two counts tables above in one list

  # **Creating the counts matrix needed for DESeq2 AND ... ---------------------------------------------------------------------  

  #creates the counts matrix: each row is a gene, the first n columns are the counts coming from the control's replicates, the last p columns are the counts from the mutant replicates
  assign( paste0( "columns_list_", name_control, "_", name_mutant), sapply(get( paste0("counts_list_", name_control, "_", name_mutant) ), `[[` , 2) )  #takes one column each two (i.e. only the columns containing the counts) from the counts list

  assign( paste0( "counts_matrix_", name_control, "_", name_mutant),  matrix(unlist( get( paste0( "columns_list_", name_control, "_", name_mutant) ) , use.names = FALSE), ncol = length( get( paste0("counts_list_", name_control, "_", name_mutant) ) ) )  )#convert the list into a matrix, for convenience

  assign( paste0( "counts_matrix_", name_control, "_", name_mutant, "_with_names"), get(  paste0( "counts_matrix_", name_control, "_", name_mutant) ) )  #the counts matrix with names is also created; at this stage is the same of the counts matrix

  # #filters out low counts
  
  # matrix_to_filter <- get(paste0( "counts_matrix_", name_control, "_", name_mutant, "_with_names"))
  # summing_rows <- apply(matrix_to_filter, 1, sum)
  # rows_to_keep <- which(summing_rows > 10)
  # matrix_filtered <- matrix_to_filter[rows_to_keep,]
  # genes_names <- names(matrix_filtered)
  # assign( paste0( "counts_matrix_", name_control, "_", name_mutant, "_with_names"), matrix_filtered)

  #** ... AND Creating the summary matrix ---------------------------------------------------------------------
  #The summary matrix contains, for each comparison, the general information about the counts statistics
  
  #extracts the needed samples from the whole list
  assign( paste0("summary_tables_", name_control), summary_tables[grep(name_control, samples_names, value = TRUE)]) #extract the controls summary tables 
  
  assign( paste0("summary_tables_", name_mutant), summary_tables[grep(name_mutant, samples_names, value = TRUE)]) #extract the mutants summary tables 
  
  assign( paste0("summary_list_", name_control, "_", name_mutant) , c( get( paste0("summary_tables_", name_control)),  get(paste0("summary_tables_", name_mutant)) )) #merges the two summary tables above in one list
  
  #creates the summary matrix: each row is a gene, the first n columns are the counts coming from the control's replicates, the last p columns are the counts from the mutant replicates
  assign( paste0( "columns_list_", name_control, "_", name_mutant), sapply(get( paste0("summary_list_", name_control, "_", name_mutant) ), `[[` , 2) )  #takes one column each two (i.e. only the columns containing the counts) from the counts list
  
  assign( paste0( "summary_matrix_", name_control, "_", name_mutant),  matrix(unlist( get( paste0( "columns_list_", name_control, "_", name_mutant) ) , use.names = FALSE), ncol = length( get( paste0("summary_list_", name_control, "_", name_mutant) ) ) )  )#convert the list into a matrix, for convenience
  
  assign( paste0( "summary_matrix_", name_control, "_", name_mutant, "_with_names"), get(  paste0( "summary_matrix_", name_control, "_", name_mutant) ) )  #the counts matrix with names is also created; at this stage is the same of the counts matrix
  
  
  col_names <-c(grep(name_control, samples_names, value=TRUE),  grep(name_mutant, samples_names, value = TRUE))  # takes the right names for the columns of the counts matrix, from the sample names 

  dummy<- get( paste0( "counts_matrix_", name_control, "_", name_mutant, "_with_names") )         #workaround for assigning dynamically colnames and rownames to the counts_matrix - START
  colnames(dummy) <- col_names
  assign(  paste0( "counts_matrix_", name_control, "_", name_mutant, "_with_names")  ,  dummy)

  dummy<- get( paste0( "counts_matrix_", name_control, "_", name_mutant, "_with_names") )
  rownames(dummy) <- genes_names
  assign(  paste0( "counts_matrix_", name_control, "_", name_mutant, "_with_names")  ,  dummy)    #workaround for assigning dynamically colnames and rownames to the counts_matrix - END

  dummy<- get( paste0( "summary_matrix_", name_control, "_", name_mutant) )         #workaround for assigning dynamically colnames and rownames to the summary_matrix - START
  colnames(dummy) <- col_names
  rownames(dummy) <- summary_tables[[1]][,1]
  assign(  paste0( "summary_matrix_", name_control, "_", name_mutant)  ,  dummy) #workaround for assigning dynamically colnames and rownames to the summary_matrix - END

  # #Evaluating the number of the reads counted and uncounted by HTSEQ2
  total_uncounted <- apply(get(paste0( "summary_matrix_", name_control, "_", name_mutant))[-1,], 2, sum)
  total_counted <- get(paste0( "summary_matrix_", name_control, "_", name_mutant))[1,]
  total_number <- total_uncounted + total_counted
  fraction_counted <- signif(total_counted/total_number, digits=3)
  void_row<-rep(" ", length(total_counted))
  statistics_matrix <- rbind(total_number, total_counted, total_uncounted,fraction_counted,void_row)
  rownames(statistics_matrix)<-c("TOTAL READS MAPPED", "TOTAL READS COUNTED", "TOTAL READS NON COUNTED", "FRACTION COUNTED" , "DETAILS UNCOUNTED:")

    #merging summary matrix and statistics matrix
  dummy<- rbind(statistics_matrix, get( paste0( "summary_matrix_", name_control, "_", name_mutant) )) 
  assign(  paste0( "summary_matrix_", name_control, "_", name_mutant)  ,  dummy)
  
  #**Defining the "experimental design" -----------------------------------
  comparisons <- sapply(col_names, function(x) substr(x,1,7))  #takes only the first parts of the column names; these are the identifiers of the control and mutant's data, irrespectively of the replicates 

  experimental_design <- data.frame(row.names = col_names, comparisons=comparisons)
  experimental_design$correlated <- c(seq(1: length(grep(name_control, col_names)) ),  seq(1: length(grep(name_mutant, col_names)) ) )
  
  experimental_design$comparisons<- relevel(experimental_design$comparisons, ref=name_control)
  
  #**Calling DESeq2 ------------------------

  #Creating the dds data structure, needed from DESEQ2
  dds <- DESeqDataSetFromMatrix(countData= get( paste0( "counts_matrix_", name_control, "_", name_mutant, "_with_names")) , colData=experimental_design, design=~comparisons)

  #Calls DEseq2 on dds and store the results in de
  de<-DESeq(dds)

  #Uses the built-in DESEQ2 results operator for creating a matrix-like structure storing the DE analysis results; the structure is stored in res_raw and subsequently converted in a data frame (res) 
  
  res_raw <-results(de)
  res <- signif(as.data.frame(res_raw), digits = sig_dig) #converts into data frame and sets the number of digits
  res <- cbind(rownames(res), res) #adds one column with the genes ID to the res dataframe (it will be useful later on)
  colnames(res) <- c("ID","Mean of norm counts", "log2 FC (MLE)", "lFC Std Err" , "Wald Stat", "Wald test pval", "BH pval") # sets the column names
  
  
  ##*Plots and Tables -----------------------------------------------------------------------------------------------------------------
  
  #**Assigns dynamic names to the table to export, taking the res data frame defined above ---------------------------------- 
  assign(paste0("results_table_",name_control,"_", name_mutant),  res ) #takes the whole res data frame
  dummy<-merge(annotation_table, get(paste0("results_table_",name_control,"_", name_mutant)), by.x="gene_id", by.y="ID")
  assign(paste0("results_table_",name_control,"_", name_mutant), dummy)
  
  
  assign(paste0("results_table_topscore_",name_control,"_", name_mutant), res[which(abs(res$`log2 FC (MLE)`)>lfc & res$`BH pval` <pval),]) #takes only the rows of the res data frame that have the best log2FC and the best BH pval
  dummy<-merge(annotation_table, get(paste0("results_table_topscore_",name_control,"_", name_mutant)), by.x="gene_id", by.y="ID")
  assign(paste0("results_table_topscore_",name_control,"_", name_mutant), dummy)
  
  
  #**Prints the tables above in .csv---------------------------------------------------------------------------------
  
  #write.table(get(paste0("results_table_topscore_",name_control,"_", name_mutant)), file=paste0("results_table_topscore_",name_control,"_", name_mutant, ".tsv"),quote=F,sep="\t", row.names = FALSE)
  
  #write.table(get(paste0("results_table_",name_control,"_", name_mutant)), file=paste0("results_table_",name_control,"_", name_mutant, ".tsv"), quote=F,sep="\t", row.names = FALSE)
  
  #write.table(get(paste0("counts_matrix_",name_control,"_", name_mutant, "_with_names")), file=paste0("counts_matrix_",name_control,"_", name_mutant, ".tsv"), quote=F,sep="\t", row.names = FALSE)
  
  #**Building the datatables ---------------------------------------------------------------------------------------------------------
  #The datatables are widgets created through the saveWidget function; this can "automatically" be sorted clicking on it --> Useful for having a general overview
  
  #***Summingup datatable ------------------------------------------
  #Builds the datatable containing all the resuls of the DE analysis
  
  #Creates the summingup_matrix, which includes only some columns of the correspondent results_table (more handy for visualisation) 
  assign(paste0("summingup_matrix_",name_control, "_", name_mutant), get(paste0("results_table_",name_control,"_", name_mutant))[, c(1:7,10,11)])
  dummy <- get(paste0("summingup_matrix_",name_control, "_", name_mutant))
  rownames(dummy) <- c()
  assign( paste0("summingup_matrix_",name_control, "_", name_mutant), dummy )
  #Creates a dynamic table (widget) of all the results
  assign(paste0("summingup_datatable_",name_control, "_", name_mutant) , datatable(get(paste0("summingup_matrix_",name_control, "_", name_mutant))) )
  #saves the datatable widget to in the working directory
  saveWidget( get(paste0("summingup_datatable_",name_control, "_", name_mutant)), file= paste0("summingup_datatable_",name_control, "_", name_mutant, ".html") )

  #***Topscores datatable ----------------------------------------- 
  #Builds the matrix containing only the top log2FC (adjusted) top scores (independently on the pvalues)
  assign(paste0("summingup_matrix_topscores_",name_control, "_", name_mutant),  subset(get(paste0("summingup_matrix_",name_control, "_", name_mutant)),  (abs(get(paste0("summingup_matrix_",name_control, "_", name_mutant))[,"log2 FC (MLE)"]) >lfc &  get(paste0("summingup_matrix_",name_control, "_", name_mutant))[,"BH pval"] <pval )  ) ) 
  dummy <- get(paste0("summingup_matrix_topscores_",name_control, "_", name_mutant))
  rownames(dummy) <- c()
  assign( paste0("summingup_matrix_topscores_",name_control, "_", name_mutant), dummy )
  #Creates a dynamic table (widget) of the topscores
  assign(paste0("topscores_datatable_",name_control, "_", name_mutant) , datatable(get(paste0("summingup_matrix_topscores_",name_control, "_", name_mutant))) )
  #saves the datatable widget to in the working directory
  saveWidget( get(paste0("topscores_datatable_",name_control, "_", name_mutant)), file= paste0("topscores_datatable_",name_control, "_", name_mutant, ".html") )

  #**Dispersion plot---------------------------------------------------------------------------------------------
  #This will not be stored or printed; for sanity check purposes only
  plotDispEsts(de, main=paste(name_mutant, " vs ", name_control))

  #**MA plot -----------------------------------------------------------------------------------------------------
  #Plot of the mean of normalised (according to DESEq) counts of the control vs. the log2fold change "corrected as well"

  #***Plot static MA -------------------------------------------------------------------------------------------------- 

  #Creates the dataframe for ggplot
  MA_dataframe <- data.frame(Ln_mean=log(res$`Mean of norm counts`), foldchange = res$`log2 FC (MLE)`, pvaladj=res$`BH pval`)
  rownames(MA_dataframe) <- rownames(res)
  MA_dataframe <- MA_dataframe[which(!is.na(MA_dataframe$foldchange)), ]
  MA_dataframe$Diff_Exp <- rep(0, nrow(MA_dataframe))
  MA_dataframe$Diff_Exp[which(abs(MA_dataframe$foldchange)>=lfc  )] <- "Relevant log2 FC" 
  MA_dataframe$Diff_Exp[which(abs(MA_dataframe$foldchange)>=lfc & MA_dataframe$pvaladj<=pval)] <- "Relevant log2 FC and Pval" 
  MA_dataframe$Diff_Exp[which(abs(MA_dataframe$foldchange)<lfc)] = "Non significant" 

  #Creates the ggplot
  r <-ggplot(MA_dataframe, aes(x=Ln_mean, y=foldchange, text=rownames(MA_dataframe)))+
    geom_point(aes(colour= Diff_Exp), size=.5)+
    geom_hline(yintercept=0, linetype=1, color="green") + geom_hline(yintercept=lfc,linetype=3, color="green") + geom_hline(yintercept=-lfc, linetype=3, color="green")+
    scale_colour_manual(values = c("Relevant log2 FC and Pval" ="red", "Non significant" = "black", "Relevant log2 FC" ="blue"))+
    xlab("ln mean of norm. counts") + ylab("log2 fold change")+ scale_x_continuous(labels = function(x)as.integer(exp(x)))+
    #+ylim(-4, 4)
    labs(aesthetic='custom text')+
    ggtitle(paste("MA plot", name_control, "vs.", name_mutant))

  #Saves the plot in a variable
  assign(paste("static", "MAplot", name_control, name_mutant, sep = "_"), r)

  #Saves the pdf of the plot in the OUTPUT directory
  pdf(file= paste("static", "MAplot", name_control, name_mutant, sep = "_"))
  print(r)
  dev.off()

  #***Plot dynamic MA -------------------------------------------------------------------------------------------------- 

  #Creates the plot
  s<- ggplotly(r, tooltip=c("text", "x", "y") )

  #SAves the plot in a variable
  assign(paste("dynamic", "MAplot", name_control, name_mutant, sep = "_"), s)

  #Saves the WIDGET of the plot in the OUTPUT directory
  print(s)
  saveWidget(s, file= paste0("dynamic_MAplot",name_control, "_", name_mutant, ".html") )

  #**Volcano Plots --------------------------------------------------------------

  #***Plot static volcano--------------------------------------------------------

  #creates a data frame with some of the columns taken from res which, in turn, summarises the DESEq2 results
  tab = data.frame(logFC = res$`log2 FC (MLE)`, negLogPval = -log10(res$`BH pval`))
  rownames(tab) <- rownames(res)

  #Identifies the DE genes
  candidate_results <- subset(res, (abs(res$`log2 FC (MLE)`) > lfc & res$`BH pval` < pval))
  assign( paste0("candidate_results_tab_",name_control,"_",name_mutant), subset(tab, (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))) )

  #Identifies the not DE genes
  non_candidate_results <- subset(res, (abs(res$`log2 FC (MLE)`) <= lfc | res$`BH pval` >= pval))
  assign(paste0("non_candidate_results_tab_",name_control,"_",name_mutant), subset(tab, (abs(tab$logFC) <= lfc | tab$negLogPval <= -log10(pval))) )

  #Mark DE and non DE genes in the dataframe 
  dummy <- get(paste0("non_candidate_results_tab_",name_control,"_",name_mutant))
  dummy$Diff_Exp <- "DE -"
  assign(paste0("non_candidate_results_tab_",name_control,"_",name_mutant), dummy)

  dummy <- get(paste0("candidate_results_tab_",name_control,"_",name_mutant))
  if(nrow(dummy)>0){dummy$Diff_Exp <- "DE +"}
  assign(paste0("candidate_results_tab_",name_control,"_",name_mutant), dummy)

  #merge the "candidate" and "non candidate" dataframes
  assign(paste0("results_tab_",name_control,"_",name_mutant), rbind(get(paste0("candidate_results_tab_",name_control,"_",name_mutant)), get(paste0("non_candidate_results_tab_",name_control,"_",name_mutant)) ))

  #builds the ggplot
  h<-ggplot( get(paste0("results_tab_",name_control,"_",name_mutant)), aes(x=logFC, y=negLogPval, text=rownames(get(paste0("results_tab_",name_control,"_",name_mutant))) ) )+
   geom_point(aes(colour= Diff_Exp), size=.5)+
   geom_hline(yintercept=-log10(pval), linetype=3, color="green") + geom_vline(xintercept=-c(-lfc, lfc), linetype=3, color="blue")+
   scale_colour_manual(values = c("DE +" ="red", "DE -" = "black"))+
   xlab("log2 fold change") + ylab("-log10 pvalue")+ 
   ggtitle(paste("Volcano plot", name_control, "vs.", name_mutant))

  #Saves the plot in a variable
  assign(paste("static", "Vplot", name_control, name_mutant, sep = "_"), h)

  #Saves the pdf of the plot in the OUTPUT directory
  pdf(file= paste("static", "Vplot", name_control, name_mutant, sep = "_"))
  print(h)
  dev.off()

  #***Plot dynamic volcano----------------------------------------------------------------------------------------------------------

  #Creates the plot
  k<- ggplotly(h, tooltip=c("text", "x", "y") )

  #Saves the plot in a variable
  assign(paste("dynamic", "Vplot", name_control, name_mutant, sep = "_"), k)

  print(k)

  #Saves the WIDGET of the plot in the OUTPUT directory
  saveWidget(k, file= paste0("dynamic_Vplot",name_control, "_", name_mutant, ".html") )

  #Updates the loop counter
  print(num_comparisons)
  num_comparisons= num_comparisons +1
} 
#*End of the loop --------------------------------------------
 
# #PRODUCES ADDITIOANL TABLES =================================================

#**Produces the Raw Counts table ---------------------------------------------------------------
raw_counts <- cbind(rownames(counts(de)), counts(de))
raw_counts_tot <- merge(annotation_table, raw_counts, by.x="gene_id", by.y="V1")
write.csv(raw_counts_tot, file="raw_counts.csv", row.names=FALSE)

#**Produces the Raw Counts Normalised table ---------------------------------------------------------------
raw_counts_normalised <- cbind(rownames(counts(de,normalized=TRUE)), counts(de, normalized=TRUE))
raw_counts_normalised_tot <- merge(annotation_table, raw_counts_normalised, by.x="gene_id", by.y="V1")
write.csv(raw_counts_normalised_tot, file="raw_counts_normalised.csv", row.names=FALSE)

#** CPM, TPM and RPKM  ---------------------------------------------------------------

#takes the genes lengths from the first feature counts input file and creates the genes_lenghts table
feature_counts_table <- fread(input_files[[1]])
genes_lenghts <- as.matrix(feature_counts_table$Length)
rownames(genes_lenghts) <- feature_counts_table$Geneid
colnames(genes_lenghts) <- "Lenght"

#merges the raw_counts matrix (all the raw counts for each sample) with the genes_lenghts_table
raw_counts_with_lenght <- (merge(genes_lenghts, counts(de), by="row.names"))
raw_counts_with_lenght <- as.matrix(raw_counts_with_lenght)
rownames(raw_counts_with_lenght) <- raw_counts_with_lenght[,1]

#eliminates the genesID column from the raw_counts_with_lenght matrix
raw_counts_with_lenght<- raw_counts_with_lenght[, -c(1)]
class(raw_counts_with_lenght) <-  "numeric"

#takes from raw_counts_with_lenghth the genes lenghts and storeas them in yy; eliminates this column from the matrix
yy<- (raw_counts_with_lenght[,1])
raw_counts_with_lenght<- raw_counts_with_lenght[, -c(1)]

# computes the CPM
CPM <- apply(raw_counts_with_lenght,2, function(x){x*10^6/sum(x)})
CPM <- merge(annotation_table, CPM, by.x="gene_id", by.y="row.names")
write.csv(CPM, file="raw_counts_CPM.csv", row.names=FALSE)

# computes the TPM
RPK <-  (raw_counts_with_lenght)/(yy/1000)  #<- Reads Per Kilobase (Normalises the counts for the genes length, expressed in kilobases)
Scaling_factor <- apply(RPK,2,sum)/10^6
TPM <- sweep(RPK, 2, Scaling_factor, FUN = '/')
TPM <- merge(annotation_table, TPM, by.x="gene_id", by.y="row.names")
write.csv(CPM, file="raw_counts_TPM.csv", row.names=FALSE)

# computes the RPKM
#RPKM<-apply(CPM, 2, h<- function(x){x*10^3/yy})

# #computes the average RPKM for controls and mutants
# assign(paste0("average_",name_control), apply(RPKM[,grep(colnames(RPKM), pattern=name_control)],1,mean ) )
# dummy<-as.matrix(get(paste0("average_",name_control)))
# colnames(dummy) <- paste0("Average RPKM ",name_control )
# assign(paste0("average_",name_control), dummy)

#computes the average CPM for controls and mutants
assign(paste0("average_",name_control), apply(CPM[,grep(colnames(CPM), pattern=name_control)],1,mean ) )
dummy<-as.data.frame(get(paste0("average_",name_control)))
colnames(dummy) <- paste0("Average CPM ",name_control)
rownames(dummy) <- CPM$gene_id
assign(paste0("average_",name_control), dummy)
#
assign(paste0("average_",names_mutants), apply(CPM[,grep(colnames(CPM), pattern=names_mutants)],1,mean ) )
dummy<-as.data.frame(get(paste0("average_",names_mutants)))
colnames(dummy) <- paste0("Average CPM ",names_mutants)
rownames(dummy) <- CPM$gene_id
assign(paste0("average_",names_mutants), dummy)
#
#merges the CPM table with the results table and the genes_length table
dummy<- merge(get(paste0("average_",name_control)),  get(paste0("average_",names_mutants)), by="row.names") 
rownames(dummy) <- dummy$Row.names
dummy<- merge(dummy[,-1], genes_lenghts, by="row.names")
#
assign(paste0("results_table_", name_control,"_", names_mutants, "_" , "CPM"), merge( get(paste0("results_table_",name_control,"_", names_mutants)), dummy, by.y="Row.names", by.x="gene_id") ) 
#
assign(paste0("results_table_topscore_", name_control,"_", names_mutants, "_", "CPM"), merge( get(paste0("results_table_topscore_",name_control,"_", names_mutants)), dummy, by.y="Row.names", by.x="gene_id") ) 
#
write.csv(get(paste0("results_table_", name_control,"_", names_mutants, "_" , "CPM")), file=paste0("results_table_", name_control,"_", names_mutants, "_" , "CPM.csv"), row.names = FALSE)
#
write.csv(get(paste0("results_table_topscore_", name_control,"_", names_mutants, "_" , "CPM")), file=paste0("results_table_topscore_", name_control,"_", names_mutants, "_" , "CPM.csv"), row.names = FALSE)

#**Filters the counts matrix  according to the CPM ------------------------

# # Filters the counts matrix eliminating the rows where at least one element violates the condition (low countys threshold)

count_zeroes<- function(x){length(which(x<=0.8))}
numberofzeroes <- apply(CPM[,c(6:ncol(CPM))], 1, count_zeroes)
CPM_filtered <- CPM[which(numberofzeroes<=2),]
# write.csv(counts_matrix_cpm_filtered, file="cpm_counts_filtered.csv")
# 
#Filters the results_mtrx with the same criteria above 
results_mtrx <- get(paste0("results_table_",name_control,"_", name_mutant, "_CPM"))
results_mtrx_filtered <- results_mtrx[which(results_mtrx$gene_id %in% CPM_filtered$gene_id), ]
write.csv(results_mtrx, file="results_table.csv", row.names = FALSE)
write.csv(results_mtrx_filtered, file="results_table_filtered.csv", row.names = FALSE)

results_mtrx_topscore <- get(paste0("results_table_topscore_",name_control,"_", name_mutant, "_CPM"))
results_mtrx_topscore_filtered <- results_mtrx_topscore[which(results_mtrx_topscore$gene_id %in% CPM_filtered$gene_id), ]
write.csv(results_mtrx_topscore_filtered, file="results_table_topscore_filtered.csv", row.names = FALSE)

# #CREATES A HEATMAP =================================================
# library(gplots)
# class(CPM_filtered) <- "numeric"
# CPM_filtered_resorted <-CPM_filtered[order(CPM_filtered[,2], decreasing = TRUE),]
# class(CPM_filtered_resorted) <- "numeric"
# CPM_filtered_resorted_log <- log(CPM_filtered_resorted)
# 
# heatmap.2(CPM_filtered_resorted_log[c(1:3000),], distfun=function(x) dist(x, method="minkowski"), hclustfun=function(x) hclust(x, method="ward.D2"), scale="row", trace = 'none', labRow = FALSE, dendrogram="column", margins=c(8,1), cexCol = .75)

#SET THE WORKING DIRECTORY BACK TO DEFAULT ================================

setwd(default_wd)

# #CREATES THE TABLE FOR GSEA =====================================================
# #The data are taken from the res_raw table.
# 
# GSEA_table <- data.frame(Sign= sign(res_raw$log2FoldChange), pval=(res_raw$padj),stringsAsFactors = FALSE)
# GSEA_table$metric <- -log10(GSEA_table$pval)/GSEA_table$Sign
# rownames(GSEA_table) <- rownames(res_raw)
# 
# GSEA_input_ranked_table <- cbind(rownames(GSEA_table), GSEA_table$metric)
# colnames(GSEA_input_ranked_table) <- cbind("Gene name", "Metric")
# GSEA_input_ranked_table<-GSEA_input_ranked_table[-which(is.na(GSEA_input_ranked_table[,"Metric"])),]
# 
# setwd("~/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Stefania_Carobbio_project/Second_dataset/Downstream_analysis/GSEA_tables/")
# 
# write.table(GSEA_input_ranked_table, file=paste0("GSEA_input_ranked_table_",name_control,"_", name_mutant, ".rnk"), quote=F,sep="\t", row.names = FALSE)


# # #GAGE ANALYSIS =====================================================
# 
# require(gage)
# require(org.Mm.eg.db)
# data("kegg.gs")
# kg.mouse<- kegg.gsets(species = "mouse", check.new = TRUE)
# kegg.gs<- kg.mouse$kg.sets
# 
# 
# #setwd("/media/davide/Dropbox/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Dimitris_project/DE_analysis_feature_counts_filtering/GAGE_analysis_1_vs_2_filtered/")
# setwd("~/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Dimitris_project/DE_analysis_feature_counts_filtering/GAGE_analysis_1_vs_2_filtered/")
# 
# anntotation_table_all <- merge(annotation_table, annotation_table_entrez, by.x="gene_id", by.y="Locus", all.x=TRUE)
# 
# annotation_table_extract <- cbind(as.vector(anntotation_table_all$gene_id),as.vector(anntotation_table_all$GeneID)) 
# colnames(annotation_table_extract) <- c("gene_id", "GeneID")
# annotation_table_extract <- unique(annotation_table_extract[which(annotation_table_extract[,1] %in% results_mtrx_filtered$gene_id),])
# 
# results_mtrx_filtered_ENTREZIDs <- merge(results_mtrx_filtered, annotation_table_extract, by.x= "gene_id", by.y="gene_id", all.x=TRUE)
# results_mtrx_filtered_ENTREZIDs <- results_mtrx_filtered_ENTREZIDs[which(!is.na(results_mtrx_filtered_ENTREZIDs$GeneID)),]
# 
# gage_input<- results_mtrx_filtered_ENTREZIDs$`log2 FC (MLE)`
# names(gage_input)<-results_mtrx_filtered_ENTREZIDs$GeneID
# 
# pathways <- gage(gage_input, gsets = kegg.gs, ref = NULL, samp = NULL)
# 
# selected_upreg <- pathways$greater[, "q.val"] < 0.3 & !is.na(pathways$greater[, "q.val"])
# 
# path.IDs.upreg <- rownames(pathways$greater)[selected_upreg]
# path.IDs.upreg_short <- substr( path.IDs.upreg, 1, 8)
# write.csv(path.IDs.upreg, file="Upreg")
# 
# 
# selected.downreg <-pathways$less[, "q.val"] < 0.3 & !is.na(pathways$less[,"q.val"])
# 
# path.IDs.downreg <- rownames(pathways$less)[selected.downreg]
# path.IDs.downreg_short <- substr( path.IDs.downreg, 1, 8)
# write.csv(path.IDs.downreg, file="Downreg")
# 
# require(pathview)
# 
# #view first 10 pathways as demo>
# 
# pdf()
# 
# pv.out.list <- sapply(path.IDs.downreg_short, function(pid) pathview(gene.data =  gage_input, pathway.id = pid, species = "mmu", out.suffix="DeSeq2", kegg.native = T ))
# pv.out.list <- sapply(path.IDs.upreg_short, function(pid) pathview(gene.data =  gage_input, pathway.id = pid, species = "mmu", out.suffix="DeSeq2", kegg.native = T ))
# 
# 
# pv.out.list <- sapply(path.IDs.downreg_short, function(pid) pathview(gene.data =  gage_input, pathway.id = pid, species = "mmu", out.suffix="DeSeq2", kegg.native = F ))
# pv.out.list <- sapply(path.IDs.upreg_short, function(pid) pathview(gene.data =  gage_input, pathway.id = pid, species = "mmu", out.suffix="DeSeq2", kegg.native = F ))
# 
# dev.off()


