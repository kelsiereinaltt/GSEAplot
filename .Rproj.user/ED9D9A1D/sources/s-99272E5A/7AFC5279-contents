#' Execute GSEA Analysis
#'
#' Analyzes genetic expression data and determines whether defined gene sets show statistically significant differences with respect to two phenotypes.
#' @param input.ds.name Name of the gct expression file
#' @param input.cls.name Name of the cls phenotype file
#' @param gene.set.input Name of the geneset file
#' @param doc.string Name of the output folder for analysis and for naming output files
#' @param nperm number of permutations
#' @param fdr.q.val.threshold significance threshold for fdr q-values
#' @param bar_percent proportional height of tick mark to window
#' @param gap_percent proportional height between minimum enrichment score and top of tick mark to the window
#' @param under_percent proportional height of white space under tick marks to the window size
#' @param upper_percent proportional height of white space over enrichment graph to the window size
#' @param color_line color of enrichment score line in plot pdf
#' @param color_tick color to tick marks on plot
#' @param abs.val Default is false. Determines whether genes are ranked according to signal to noise or absolute value of signal to noise (when abs.val=T)
#' @param gs.size.threshold.max Default is 1000. Maximum matches between geneset and gene labels
#'
#' @details GSEA analysis is computed using the Broad Institute's R source code. Genes are ranked according to signal to noise ratio (difference in means/sum of standard deviations for the two phenotypes)
#' @return pp a list pp which includes : plots, gene.set.reference.matrix, gene.set.leading, report1, report2, ES
#' @return plots a list of ggplot objects, this can act as an input to the function plot.ES() to output a pdf
#' @return gene.set.reference.matrix a list of each gene set within a gene set database and the gene symbols corresponding to each set
#' @return gene.set.leading a similar structure to gene.set.reference.matrix but only contains the gene symbols within each gene set that are part of the leading edge set
#' @return report1 summary of GSEA analysis data for the first phenotype
#' @return report2 summary of GSEA analysis data for the second phenotype
#' @return ES this object contains the enrichment scores and enrichment tags used to create the plots described earlier. The user can use this information to customize plots as they wish
#' @export
#' @references Subramanian, Tamayo, et al. (2005), PNAS 102, 15545-15550, http://www.broad.mit.edu/gsea/

GSEAplots= function(input.ds.name="",input.cls.name="", gene.set.input="",
                    doc.string="", nperm=1000,fdr.q.val.threshold = 0.25,bar_percent=0.1,
                    gap_percent=0.1, under_percent=0.1,upper_percent=0.1,color_line="black", color_tick="black",abs.val=F,gs.size.threshold.max=1000){
  # GSEA 1.0 -- Gene Set Enrichment Analysis / Broad Institute

 # Inputs:
    #   program_location: the outermost directory. The one that contains the source code and code to run the file
    #   input.ds.name: Input gene expression Affymetrix dataset file in GCT format "name.gct"
    #   input.cls.name:  Input class vector (phenotype) file in CLS format "name.cls"
    #   gene.set.input: Gene set name for hallmark it is "h.all.v6.1.symbols.gmt"
    #   doc.string: Dataset_geneset
    #   nperm: number of permutations, default is 1000
    #   nom.p.val.threshold Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
    #   fdr.q.val.threshold   = 0.25,   Significance threshold for FDR q-vals for gene sets (default: 0.25)



  wd_new=getwd()
  if (file.exists(paste0(wd_new,"/", doc.string))==FALSE){
  dir.create(doc.string)
  }
  results_new=GSEA(                                                                    # Input/Output Files :-------------------------------------------
                                                                          input.ds=input.ds.name, # Input gene expression Affy dataset file in RES or GCT format
                                                                          input.cls=input.cls.name,
                                                                          #input.ds=paste(wd_new,datasets.folder,input.ds.name, sep="",collapse=NULL), # Input gene expression Affy dataset file in RES or GCT format
                                                                           #input.cls=paste(wd_new,datasets.folder,input.cls.name,sep="",collapse=NULL),
                                                                           # gs.db =   paste(genesets.folder,gene.set.input,sep="",collapse=NULL),         # Gene set database in GMT format
                                                                          gs.db =   gene.set.input,
                                                                          output.directory      = paste0(wd_new,"/", doc.string,"/"),
                                                                           output.directory2      =paste0(wd_new,"/"),
                                                                                # Directory where to store output and results (default: "")
                                                                           #  Program parameters :-------------------------------------------------------------------------------------------------------------------------
                                                                           doc.string            = doc.string,   # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
                                                                           non.interactive.run   = T,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
                                                                          reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
                                                                          nperm                 = nperm,            # Number of random permutations (default: 1000)
                                                                           weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
                                                                           nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
                                                                           fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
                                                                           fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
                                                                           topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
                                                                           adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
                                                                           gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
                                                                           gs.size.threshold.max = gs.size.threshold.max,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
                                                                           reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
                                                                           preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
                                                                           random.seed           = 3338,            # Random number generator seed. (default: 123456)
                                                                           perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
                                                                           fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
                                                                           replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
                                                                           save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
                                                                           OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
                                                                           use.fast.enrichment.routine = T,          # Use faster routine to compute enrichment for random permutations (default: T)
                                                                           abs.val=abs.val                               #rank by absolute value of signal to noise ratio
  )

  #-----------------------------------------------------------------------------------------------------------------------------------------------

  if (length(which(sapply(results_new$out5,is.null))) == 0){
    ES.tags.files <- results_new$out5
    ES.data.files <- results_new$out6
    ES.report.files <- results_new$out7
    gene.set.numbers <- results_new$out4
  } else{
    #this step removes the gene sets which did not generate ES.tags, ES.data, or ES.report files
  ES.tags.files <- results_new$out5[-which(sapply(results_new$out5,is.null))]
  ES.data.files <- results_new$out6[-which(sapply(results_new$out6,is.null))]
  ES.report.files <- results_new$out7[-which(sapply(results_new$out7,is.null))]
  gene.set.numbers <- results_new$out4[-which(sapply(results_new$out5,is.null))]
  }
  gene.set.reference.matrix <- results_new$gene.set.reference.matrix
  gene.set.leading <- rep(list("null"),length(gene.set.numbers))
  ES <- rep(list("null"),length(gene.set.numbers))
  enrichind <- rep(list("null"),length(gene.set.numbers))
  report1 <- results_new$report1
  report2 <- results_new$report2


  if (regexpr(pattern="HALLMARK_", gene.set.numbers[1]) == -1) {
  #nothing
  } else {
    for (i in 1: length(gene.set.numbers)){
      g <- strsplit(gene.set.numbers[[i]],split="HALLMARK_")
      h <- g[[1]][2]
      gene.set.numbers[[i]] <- h
    }
  }


  #plotting and finding leading edge set
  plots <- vector(mode="list",length=length(ES.data.files))
  for (i in 1:length(ES.tags.files)){
     dat1_name=paste(wd_new,"/",doc.string,"/",doc.string,".",ES.data.files[[i]],sep="",collapse=NULL)
     dat2_name=paste(wd_new,"/",doc.string,"/",doc.string,".",ES.tags.files[[i]],sep="",collapse=NULL)
     report_name=paste(wd_new,"/",doc.string,"/",doc.string,".",ES.report.files[[i]],sep="",collapse=NULL)
    dat1 = read.table(dat1_name,header=T,sep="\t")
    dat2 = read.table(dat2_name,header=T,sep="\t")
    report=read.table(report_name,sep="\t")
    datcomb <- cbind(dat1,dat2)
    datcomb= datcomb[,-5]
    ES[[i]]=datcomb
    enrich_ind=which(dat2$EStag==1)
    height=max(dat1$RES)-min(dat1$RES)
    bar_height=bar_percent*height
    gap_height=gap_percent*height
    under_height=under_percent*height
    upper_height=upper_percent*height
    window_height=height+bar_height+gap_height+under_height+upper_height
    y_lower=min(dat1$RES)-gap_height-bar_height
    window_low=y_lower-under_height
    window_high=max(dat1$RES)+upper_height
    d=data.frame(x=enrich_ind, y=matrix(y_lower,length(enrich_ind),1), vx=matrix(0,length(enrich_ind),1), vy=matrix(bar_height,length(enrich_ind),1))
    p <- ggplot(datcomb, aes(index,RES))+geom_line(col=color_line)
    p <- p+geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),col=color_tick)
    p <- p+theme_classic()
    p <- p+ylim(window_low,window_high)
    p <- p+ggtitle(gene.set.numbers[[i]])
    p <- p+theme(plot.title = element_text(size = 12))
    plots[[i]] <- p

    #location of largest RES absolute value
    leading_index <- which.max(abs(datcomb$RES))

    #depending on which phenotype the gene set is more related to the leading edge set is on different sides
    if (datcomb$RES[leading_index] > 0){
      #takes the report of all enriched genes and selects the ones before the max
      leading.edge.names <- report$V2[which(report$V5 < leading_index)]
      leading.edge.RES <- report$V7[which(report$V5 < leading_index)]
    }
    else {
      #takes the report of all enriched genes and selects the ones after the max
      leading.edge.names <- report$V2[which(report$V5 > leading_index)]
      leading.edge.RES <- report$V7[which(report$V5 > leading_index)]
    }

    leading.edge.set <- data.frame(leading.edge.names,leading.edge.RES)
    gene.set.leading[[i]] <- leading.edge.names

  }


  names(gene.set.leading) <- gene.set.numbers
  names(ES) <- gene.set.numbers
  gene.set.leading[] <- lapply(gene.set.leading, as.character)


  return(list(plots=plots,gene.set.reference.matrix=gene.set.reference.matrix,gene.set.leading=gene.set.leading,report1=report1,report2=report2,ES=ES))
}







