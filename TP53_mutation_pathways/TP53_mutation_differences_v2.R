setwd("/Users/davidshorthouse/OneDrive - MRC Cancer Unit at University of Cambridge/Metabolism_Fellowship/MCF7_removal/TP53_mutation_pathways/")

library("MetaboDiff")
library("phytools")
library("colorspace")
library("RColorBrewer")

### Load bulk data
assaydata <- read.csv("./assaydata_tp53mutation.csv")
row.names(assaydata) <- assaydata$ionIdx
assaydata[1] <- NULL
coldata <- read.csv("./coldata_tp53mutation.csv")
row.names(coldata) <- coldata$ionIdx
coldata[1] <- NULL
rowdata <- read.csv("./rowdata_tp53mutation.csv")
rowdata<- rowdata[!duplicated(rowdata$ionIdx), ]
met <- create_mae(assaydata, rowData = rowdata, colData = coldata)

### Normalise data
met <- create_mae(assaydata, rowData = rowdata, colData = coldata)
met <- get_SMPDBanno(met, column_hmdb = 5, column_kegg_id =  NA, column_chebi_id = NA)
met = knn_impute(met, cutoff = 0.4)
met <- normalize_met(met)

dff_met = diff_test(met,
                    group_factors = c("TP53_status"))
met_pathways <- dff_met %>%
  diss_matrix %>%
  identify_modules(min_module_size=5) %>%
  name_modules(pathway_annotation="SMPDB.Pathway.Name") %>%
  calculate_MS(group_factors=c("TP53_status"))

##### dff_met now contains the data to generate trees - we want pvalues for individual mutations

# Conf_mutants
### Do conformation vs WT
conf_assaydata <- read.csv("./Conformation_assaydata_tp53mutation.csv")
row.names(conf_assaydata) <- conf_assaydata$ionIdx
conf_assaydata[1] <- NULL
conf_coldata <- read.csv("./Conformation_coldata_tp53_mutation.csv")
row.names(conf_coldata) <- conf_coldata$ionIdx
conf_coldata[1] <- NULL
conf_met <- create_mae(conf_assaydata, rowData = rowdata, colData = conf_coldata)
conf_met <- get_SMPDBanno(conf_met, column_hmdb = 5, column_kegg_id =  NA, column_chebi_id = NA)
conf_met = knn_impute(conf_met, cutoff = 0.4)
conf_met <- normalize_met(conf_met)

conf_met_diff = diff_test(conf_met,
                          group_factors = c("TP53_status"))

# Contact mutants
### Do contact mutants vs WT
cont_assaydata <- read.csv("./Contact_assaydata_tp53mutation.csv")
row.names(cont_assaydata) <- cont_assaydata$ionIdx
cont_assaydata[1] <- NULL
cont_coldata <- read.csv("./Contact_coldata_tp53_mutation.csv")
row.names(cont_coldata) <- cont_coldata$ionIdx
cont_coldata[1] <- NULL
cont_met <- create_mae(cont_assaydata, rowData = rowdata, colData = cont_coldata)
cont_met <- get_SMPDBanno(cont_met, column_hmdb = 5, column_kegg_id =  NA, column_chebi_id = NA)
cont_met = knn_impute(cont_met, cutoff = 0.4)
cont_met <- normalize_met(cont_met)

cont_met_diff = diff_test(cont_met,
                          group_factors = c("TP53_status"))

####### Now have pvalues for individual mutants and tree from bulk data #######

calculate_MS2 = function(met, group_factors){
  for (i in 1:length(group_factors)) {
    id = grep(group_factors[i],names(metadata(met)))[1]
    df = metadata(met)[[id]]
    df$modules = metadata(met_pathways)$modules
    res_df = plyr::ddply(df, "modules", plyr::summarise, av_pval=mean(pval),
                         av_adj_pval = mean(adj_pval))
    res_df$av_dm = rep(0,nrow(res_df))
    for(x in 1:nrow(df)){
      if(sum(df$modules==(x-1)&df$pval<0.05)>0){
        res_df$av_dm[x] = median(df[df$modules==(x-1)&df$pval<0.05,"dm"],na.rm=TRUE)
      }
    }
    metadata(met)[[paste0("MS_",group_factors[i])]] = res_df
  }
  met
}

### Calculate avg pvalues but with bulk data pathway assignments
conf_met_pathways <- conf_met_diff %>%
  diss_matrix %>%
  identify_modules(min_module_size=5) %>%
  name_modules(pathway_annotation="SMPDB.Pathway.Name") %>%
  calculate_MS2(group_factors=c("TP53_status"))


cont_met_pathways <- cont_met_diff %>%
  diss_matrix %>%
  identify_modules(min_module_size=5) %>%
  name_modules(pathway_annotation="SMPDB.Pathway.Name") %>%
  calculate_MS2(group_factors=c("TP53_status"))

#### Try to plot ####
## Function for desaturating colors by specified proportion
desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

rainbowdesat <- desat(rainbow(100,start=0, end =0.65, rev= TRUE), 0.5)
par(mar = rep(0, 4))
pie(rep(1, length(rainbowdesat)), col = rainbowdesat, border = NULL)



p_value_cutoff = 0.05

## Prepare the data
image(1:nrow(rainbowdesat), 1, as.matrix(1:nrow(rainbowdesat)), 
      col=rainbowdesat,
      xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")

id1 = grep("TP53_status",names(metadata(conf_met_pathways)))[2]
id2 = grep("TP53_status",names(metadata(cont_met_pathways)))[2]
tree = ape::as.phylo(metadata(met_pathways)$METree)

ungrouped <- tree$tip.label[c(26, 27, 32, 42, 47, 49, 50, 55, 56, 62, 64, 68, 70, 72, 76, 78, 86, 88, 90, 91, 92, 93, 95, 96)]
pruned.tree <-drop.tip(tree,tree$tip.label[match(ungrouped, tree$tip.label)])

pruned_metadata <- metadata(conf_met_pathways)[[id1]][-c(26, 27, 32, 42, 47, 49, 50, 55, 56, 62, 64, 68, 70, 72, 76, 78, 86, 88, 90, 91, 92, 93, 95, 96),]
x1 = -log10(metadata(conf_met_pathways)[[id1]]$av_pval)
x2 = -log10(metadata(cont_met_pathways)[[id2]]$av_pval)
names(x1) = tree$tip.label
names(x2) = tree$tip.label

lims<-range(c(min(x1,x2),max(x1,x2))) ## Get plot limits

lims <- range(c(0,2))


map1<-contMap(tree,x1,plot=FALSE, lims = lims)
map1<-setMap(map1, colors = rev(brewer.pal(n = 11, name = "RdYlBu")))## invert color map

map2<-contMap(tree,x2,plot=FALSE, lims = lims)
map2<-setMap(map2, colors = rev(brewer.pal(n = 11, name = "RdYlBu")))## invert color map

pruned.map1 <-drop.tip.contMap(map1,ungrouped)
pruned.map2 <-drop.tip.contMap(map2,ungrouped)



## dimensions are in inches
pdf("Multilevel_phyloplot_pruned_v2.pdf", width = 8, height = 6)
layout(matrix(1:3,1,3),widths=c(0.3,0.4,0.3))
par(cex=1.1) ## make sure the correct font size is used in subplots


plot(pruned.map1,fsize=c(0,0.4),ftype=c("off","reg"),sig=1,legend=5,
     xlim=c(-0,1)*max(nodeHeights(pruned.tree)),mar=c(1.1,0.1,4.1,0.1))



ylim<-c(1-0.12*(length(pruned.tree$tip.label)-1),length(pruned.tree$tip.label))
plot.new()
plot.window(xlim=c(-0.1,0.1),ylim=ylim)
text(rep(0,length(pruned.tree$tip.label)), 1:length(pruned.tree$tip.label),sapply(strsplit(pruned.tree$tip.label, '\\|'), function(x) paste(x[-1], collapse = '|')),
     font=2, cex = 0.3)
### Need to remove the first element from pruned.tree$tip.label split on |


ylim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim
plot(pruned.map2,fsize=c(0,0.4),ftype=c("off","reg"),
     direction="leftwards", sig=1,legend=F,
     xlim=c(-0,1)*max(nodeHeights(tree)),mar=c(1.1,0.1,4.1,0.1), ylim = ylim)


dev.off()




pdf("Multilevel_phyloplot_unpruned.pdf")
layout(matrix(1:3,1,3),widths=c(0.35,0.3,0.35))
par(cex=1.1) ## make sure the correct font size is used in subplots


plot(map1,fsize=c(0,0.4),ftype=c("off","reg"),sig=1,legend=5,
     xlim=c(-0,1)*max(nodeHeights(pruned.tree)),mar=c(1.1,0.1,4.1,0.1))



ylim<-c(1-0.12*(length(tree$tip.label)-1),length(tree$tip.label))
plot.new()
plot.window(xlim=c(-0.1,0.1),ylim=ylim)
text(rep(0,length(tree$tip.label)), 1:length(tree$tip.label),rev(tree$tip.label),
     font=3, cex = 0.3)

ylim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim
plot(map2,fsize=c(0,0.4),ftype=c("off","reg"),
     direction="leftwards", sig=1,legend=F,
     xlim=c(-0,1)*max(nodeHeights(tree)),mar=c(1.1,0.1,4.1,0.1), ylim = ylim)


dev.off()



