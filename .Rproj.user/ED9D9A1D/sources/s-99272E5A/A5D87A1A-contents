# library(GSEA.plot)
# data(hallmark.gs)
#
# setwd("C:/Users/Student/Documents/CivelekLab/August_Edits/")
# fname = "f_KLF14trans.txt"
# transf = read.table(fname,header=F,stringsAsFactors=F)
# use_data(transf)
# fname = "m_KLF14trans.txt"
# transm = read.table(fname,header=F,stringsAsFactors=F)
# use_data(transm)
# tx = c("HALLMARK_KLF14","http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_KLF14",transf[,1], transm[,1])
# tz=c("BS","this source",transf[,1],transm[,1])
# ty=rbind(tx,tz)
#
# test=add_to_database(database=hallmark.gs,addition=tx)
#
# database=hallmark.gs
# addition=ty
#
# test2=add_to_database(database=hallmark.gs, addition=ty)
# gene.set.input=test2
#
# #####
# # data
# fname = "expr_gtexSubq_gsea.txt"
# expr = read.table(fname,header=T,stringsAsFactors=F)
# fname = "lab_gtexSubq_gsea.txt"
# ann = read.table(fname,header=F,stringsAsFactors=F)[1,] %>% t
# ann = as.character(ann)
# pheno.input=create_phenoinput(ann)
#
# # format expr
# expr2 = expr[!(duplicated(expr[,1])==T),]
# expr.input = expr2[,3:ncol(expr2)]
# rownames(expr.input) = expr2[,1]
#
# pp = GSEAplots(input.ds.name = expr.input, input.cls.name = pheno.input,
#                gene.set.input = gene.set.input, doc.string = "gtex6p",
#                nperm = 1000, gs.size.threshold.max=1000, fdr.q.val.threshold = 0.25,
#                bar_percent = 0.1, gap_percent = 0.1, under_percent = 0.1,
#                upper_percent = 0.1, color_line = "black", color_tick = "black",
#                abs.val = F, gs.size.threshold.max =1000)
#
# plot.ES(list.of.plots=pp$plots,plotname="GSEA_plots")
#


#fixing get_gene_symbols
 # sets=get_genesets(hallmark.gs)
 # symbols=get_genesymbols(hallmark.gs)
 # entry1=c("description",symbols$HALLMARK_TNFA_SIGNALING_VIA_NFKB)
# entry2=c("description",symbols$HALLMARK_HYPOXIA)
# entry3=c("description",symbols$HALLMARK_CHOLESTEROL_HOMEOSTASIS)
# entry=list(entry1,entry2,entry3)
# names(entry)=c(sets[1],sets[2],sets[3])
# entry2=create_geneset_db(entry)
# gene.set.input=create_geneset_db(entry)
#
#
# library(GSEA.plot)
# data(transf)
# data(transm)
# tx = c("http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_KLF14",transf[,1], transm[,1])
# ty=list(tx)
# names(ty)="HALLMARK_KLF14"
# gene.set.input=create_geneset_db(ty)
#
# data(gtex_ann)
# pheno.input=create_phenoinput(gtex_ann)
#
# data(gtex_expr)
# expr.input=gtex_expr
#
# #test with one input  that is from hallmark
# sets=get_genesets(hallmark.gs)
# symbols=get_genesymbols(hallmark.gs)
# entry1=c("description",symbols$HALLMARK_TNFA_SIGNALING_VIA_NFKB)
# entry1=list(entry1)
# names(entry1)=sets[1]
# gene.set.input=create_geneset_db(entry1)
#
# gtex_example= GSEAplots(input.ds.name=expr.input,
#                         input.cls.name=pheno.input, gene.set.input=gene.set.input,
#                         doc.string="Gtex_hall", nperm=1000,gs.size.threshold.max = 1000,
#                         fdr.q.val.threshold = 0.25,abs.val=F,bar_percent=0.1, gap_percent=0.1,
#                         under_percent=0.1,upper_percent=0.1,color_line="black",
#                         color_tick="black")
#
# plot.ES(list.of.plots=gtex_example$plots,plotname="GSEA_plots")
