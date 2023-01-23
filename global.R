library(rtracklayer) 
library(Rsubread)
library(Rsamtools)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(shiny)
library(DT)
library(gdata)
library(rstatix)
library(multcomp)
library(tidyverse)
library(tools)
library(ggpubr)
library(ggrepel)
library(ggdendro)
library(ggplotify)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(ggnewscale)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(msigdbr)
library(genefilter)
library(ComplexHeatmap)
library(shinyBS, verbose = FALSE)
library(plotly,verbose=FALSE)
library('shinyjs', verbose = FALSE)
library(BiocManager)
library(clusterProfiler.dplyr)
library(dorothea)
library(GenomicRanges)
library(soGGi) ##devtools::install_github("ColeWunderlich/soGGi")
library(ChIPseeker)
library(ChIPpeakAnno)
library(rGREAT)
library(FindIT2)
library(limma)
library(ggseqlogo)
library(marge)
library(plotmics)
library(colorspace)
library(ggcorrplot)
library(RColorBrewer)
library(bedtorch)
options('homer_path' = "/usr/local/homer")
check_homer()
options(rsconnect.max.bundle.size=31457280000000000000)
species_list <- c("not selected", "Mus musculus (mm10)","Homo sapiens (hg19)","Homo sapiens (hg38)")
gene_set_list <- c("MSigDB Hallmark", "KEGG", "Reactome", "PID (Pathway Interaction Database)",
                   "BioCarta","WikiPathways", "GO biological process", 
                   "GO cellular component","GO molecular function", "Human phenotype ontology", 
                   "DoRothEA regulon (activator)", "DoRothEA regulon (repressor)",
                   "Transcription factor targets", "miRNA target")
gene_set_list_genome <- c("MSigDB Hallmark", "KEGG", "Reactome", "PID (Pathway Interaction Database)",
                   "BioCarta","WikiPathways", "GO biological process", 
                   "GO cellular component","GO molecular function", "Human phenotype ontology", 
                   "Transcription factor targets", "miRNA target")
org <- function(Species){
  if(Species != "not selected"){
    switch (Species,
            "Homo sapiens (hg38)" = org <- org.Hs.eg.db,
            "Homo sapiens (hg19)" = org <- org.Hs.eg.db,
            "Mus musculus (mm10)" = org <- org.Mm.eg.db
            )
    return(org)
  }
}

txdb_function <- function(Species){
  switch (Species,
          "Mus musculus (mm10)" = txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene,
          "Homo sapiens (hg19)" = txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene,
          "Homo sapiens (hg38)" = txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene)
  return(txdb)
}

promoter <- function(txdb, upstream, downstream,input_type = "Promoter",files = NULL, bam=F,RPM){
  if(input_type == "Promoter"){
  return(promoters(genes(txdb),upstream = upstream, downstream = downstream))
  }else{
    consensusToCount <- soGGi:::runConsensusRegions(GRangesList(files), "none")
    occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% 
      rowSums
    consensusToCount <- consensusToCount[occurrences >= 2, ]
    return(consensusToCount)
  }
}
promoter_clustering <- function(txdb, upstream, downstream,input_type = "Promoter",files = NULL, bam=F,RPM){
  if(input_type == "Promoter"){
    return(promoters(genes(txdb),upstream = upstream, downstream = downstream))
  }else{
    if(length(names(files)) > 1){
    consensusToCount <- soGGi:::runConsensusRegions(GRangesList(files), "none")
    }else consensusToCount <- files[[1]]
    return(consensusToCount)
  }
}

Bigwig2count <- function(bw, promoter, Species, input_type = "Promoter"){
bed1<-promoter
write.table(bed1,file = "bed.bed",sep = "\t")
bed1 <- read.table("bed.bed",header = T)
bed1 <- dplyr::arrange(bed1, seqnames)
data <- as.data.frame(bed1)[,1:3]
colnames(data)<-c("chr","start","end")
bed<-with(data,GRanges(chr,IRanges(start,end)))
counts <- matrix(NA, nrow = length(bed), ncol = length(bw))
colnames(counts) <- names(bw)
chromnames=levels(seqnames(bed))
perc <- 0
withProgress(message = "converting BigWig to gene count data",{
for(i in seq_len(length(bw))) {
    perc <- perc + 1
  last=1 	 
  coverage <- import(bw[[i]], as = 'RleList')
  for(j in chromnames){
    range_vals=ranges(bed[seqnames(bed)==j])
    cur_coverage=coverage[[j]] 
    if(is.null(cur_coverage)){
      counts[last:(last+length(range_vals)-1),i]=matrix(0,nrow=length(range_vals),ncol=1)
    }else{
      newvals=sum(Views(cur_coverage, ranges(bed[seqnames(bed)==j])))
      counts[last:(last+length(newvals)-1), i] <-newvals 
    }   
    last=last+length(range_vals) 
  }
  incProgress(1/length(colnames(counts)), message = paste("Finish '", colnames(counts)[i], "', ", 
                                            perc, "/", length(colnames(counts)),sep = ""))
}
})
if(input_type == "Promoter"){
  switch (Species,
          "Mus musculus (mm10)" = org <- org.Mm.eg.db,
          "Homo sapiens (hg19)" = org <- org.Hs.eg.db,
          "Homo sapiens (hg38)" = org <- org.Hs.eg.db)
  rownames(counts) <- bed1$gene_id
my.symbols <- rownames(counts)
gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                keytype = "ENTREZID",
                                columns = c("ENTREZID","SYMBOL"))
colnames(gene_IDs) <- c("Row.names","SYMBOL")
gene_IDs <- gene_IDs %>% distinct(SYMBOL, .keep_all = T) %>% na.omit()
gene_IDs <- data.frame(SYMBOL = gene_IDs$SYMBOL, row.names = gene_IDs$Row.names)
data2 <- merge(gene_IDs,counts, by=0)
rownames(data2) <- data2$SYMBOL
counts <- data2[,-1:-2]
}else{
  a <- as.data.frame(bed)
  Row.name <- paste0(a$seqnames,":",a$start,"-",a$end)
  rownames(counts) <- Row.name
}
return(counts)
}



PCAplot <- function(data, plot){
  if(length(grep("SYMBOL", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "SYMBOL")]
  }
  if(length(grep("Unique_ID", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "Unique_ID")]
  }
  pca <- prcomp(data, scale. = T)
  label<- colnames(data)
  label<- gsub("\\_.+$", "", label)
  lab_x <- paste(summary(pca)$importance[2,1]*100,
                 "% of variance)", sep = "")
  lab_x <- paste("PC1 (", lab_x, sep = "")
  lab_y <- paste(summary(pca)$importance[2,2]*100,
                 "% of variance)", sep = "")
  lab_y <- paste("PC2 (", lab_y, sep = "")
  pca$rotation <- as.data.frame(pca$rotation)
  if(plot==TRUE){
  g1 <- ggplot(pca$rotation,aes(x=pca$rotation[,1],
                                y=pca$rotation[,2],
                                col=label, label = label)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA)) +
    xlab(lab_x) + ylab(lab_y) + geom_text_repel()  +
    theme(legend.position="none", aspect.ratio=1)
  rho <- cor(data,method="spearman")
  d <- dist(1-rho)
  mds <- as.data.frame(cmdscale(d))
  label<-colnames(data)
  label<-gsub("\\_.+$", "", label)
  g2 <- ggplot(mds, aes(x = mds[,1], y = mds[,2],
                        col = label, label = label)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA)) +
    xlab("dim 1") + ylab("dim 2") +
    geom_text_repel() + theme(legend.position="none", aspect.ratio=1)
  x <- NULL
  y <- NULL
  xend <- NULL
  yend <- NULL
  data.t <- t(data)
  hc <- hclust(dist(data.t), "ward.D2")
  dendr <- dendro_data(hc, type="rectangle")
  g3 <- ggplot() +
    geom_segment(data=segment(dendr),
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=label(dendr),
              aes(x, y, label=label, hjust=0), size=3) +
    theme(legend.position = "none",
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),axis.text.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          axis.title.y=element_blank(),panel.background=element_rect(fill="white"))+
    coord_flip()+ scale_y_reverse(expand=c(0.6, 0))
  p2 <- plot_grid(g1, g2, g3, nrow = 1)
  return(p2)
  }else return(pca$rotation)
}
umap_plot <- function(data, n_neighbors){
  umap <- umap::umap(t(data),n_neighbors = n_neighbors, random_state = 123)
  data2 <- umap$layout %>% as.data.frame()
  label<- colnames(data)
  label<- gsub("\\_.+$", "", label)
  p<- ggplot(data2, mapping = aes(V1,V2, color = label, label = colnames(data)))+
    geom_point()+geom_text_repel()+ xlab("UMAP_1") + ylab("UMAP_2")+
    theme(panel.background =element_rect(fill=NA,color=NA),panel.border = element_rect(fill = NA),
          aspect.ratio=1)
  return(p)
}

GOIheatmap <- function(data.z, show_row_names = TRUE){
  ht <- Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                clustering_method_columns = 'ward.D2',
                show_row_names = show_row_names, show_row_dend = F,column_names_side = "top",
                row_names_gp = gpar(fontface = "italic"))
  return(ht)
}

pdfSize_for_GOI <- paste(strong("Heatmap:"),"height = 10, width = 7 <br>", 
                         strong("Boxplot:"),"<br>",
                         "Gene number = 1,","height = 3, width = 3 <br>",
                         "Gene number = 2,","height = 3, width = 6 <br>",
                         "Gene number = 3,","height = 3, width = 9 <br>",
                         "Gene number = 4,","height = 6, width = 6 <br>",
                         "Gene number = 5 ~ 6,","height = 6, width = 9 <br>",
                         "Gene number = 7 ~ 9,","height = 7.5, width = 6.75 <br>",
                         "Gene number = 10 ~ 12,","height = 7.5, width = 9 <br>",
                         "Gene number = 13 ~ 16,","height = 9, width = 9 <br>",
                         "Gene number = 17 ~ 25,","height = 11.5, width = 11.5 <br>",
                         "Gene number = 26 ~ 36,","height = 13.5, width = 13.5 <br>",
                         "Gene number = 37 ~ 49,","height = 15.75, width = 15.75 <br>",
                         "Gene number = 50 ~ 64,","height = 18, width = 18 <br>",
                         "Gene number = 65 ~ 81,","height = 20.5, width = 20.5 <br>",
                         "Gene number = 82 ~ 200,","height = 22.5, width = 22.5 <br>",
                         "Gene number > 200,", "height = 30, width = 30 <br>")

dotplot_for_output <- function(data, plot_genelist, Gene_set, Species){
  if(!is.null(Gene_set)){
    if(is.null(data)){
      return(NULL)
    }else{
      withProgress(message = "Plot results",{
        if(Species != "not selected"){
          print(plot_genelist)
        }
        incProgress(1)
      })
    }
  }
}
GeneList_for_enrichment <- function(Species, Gene_set, org, Custom_gene_list){
  if(Species != "not selected" || is.null(Gene_set) || is.null(org)){
    switch (Species, 
            "Mus musculus (mm10)" = species <- "Mus musculus",
            "Homo sapiens (hg19)" = species <- "Homo sapiens",
            "Homo sapiens (hg38)" = species <- "Homo sapiens")
    if(Gene_set == "MSigDB Hallmark"){
      H_t2g <- msigdbr(species = species, category = "H") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description) 
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HALLMARK_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="P53", replacement = "p53")
    }
    if(Gene_set == "KEGG"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="KEGG_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "Transcription factor targets"){
      H_t2g <- msigdbr(species = species, category = "C3")
      H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
    if(Gene_set == "DoRothEA regulon (activator)"){
      H_t2g <- as_tibble(dorothea(species = Species,  type = "DoRothEA regulon (activator)")) %>%
        dplyr::select(gs_name, entrez_gene, confidence)
      H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
    }
    if(Gene_set == "DoRothEA regulon (repressor)"){
      H_t2g <- as_tibble(dorothea(species = Species,  type = "DoRothEA regulon (repressor)")) %>%
        dplyr::select(gs_name, entrez_gene, confidence)
      H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
    }
    if(Gene_set == "Reactome"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="REACTOME_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "miRNA target"){
      H_t2g <- msigdbr(species = species, category = "C3")
      H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "MIR:MIRDB" | gs_subcat == "MIR:MIR_Legacy") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
    if(Gene_set == "GO biological process"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "GO:BP") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOBP_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "GO cellular component"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "GO:CC") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOCC_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "GO molecular function"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "GO:MF") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOMF_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "Human phenotype ontology"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "HPO") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HP_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "WikiPathways"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="WP_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "PID (Pathway Interaction Database)"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:PID") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="PID_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "BioCarta"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:BIOCARTA") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="BIOCARTA_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "Custom gene set"){
      if(!is.null(Custom_gene_list)){
        H_t2g <- gene_list_convert_for_enrichment(data= Custom_gene_list, Species = Species)
        H_t2g <- data.frame(gs_name = H_t2g$Group, entrez_gene = H_t2g$ENTREZID)
        H_t2g$gs_name <- gsub(":", "_", H_t2g$gs_name)
      }else H_t2g <- NULL
    }
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Tnf", replacement = "TNF")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Tgf", replacement = "TGF")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Pi3k_akt_mtor", replacement = "PI3K_Akt_mTOR")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Il6_", replacement = "IL6_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Il2_", replacement = "IL2_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Kras", replacement = "KRas")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Uv_r", replacement = "UV_r")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Dna_", replacement = "DNA_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Rna_", replacement = "RNA_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Mrna_", replacement = "mRNA_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="E2f", replacement = "E2F")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="G2m", replacement = "G2M")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Mtorc", replacement = "mTORC")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Ecm_", replacement = "ECM_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Abc_", replacement = "ABC_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="No1_", replacement = "NO1_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_mirna", replacement = "_miRNA")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="gtpase_", replacement = "GTPase_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="rho_", replacement = "Rho_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="jnk_", replacement = "JNK_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_jak", replacement = "_JAK")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_stat", replacement = "_STAT")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_nfkb", replacement = "_NFkB")
    return(H_t2g)
  }else return(NULL)
}
gene_list_for_enrichment_genome <- function(H_t2g){
  df <- list()
  set <- unique(H_t2g$gs_name)
  for(name in set){
    data <- dplyr::filter(H_t2g, gs_name == name)
    df[[name]] <- data$entrez_gene
  }
  return(df)
}
dorothea <- function(species, confidence = "recommend",type){
  if(species == "Mus musculus (mm10)"){
    net <- dorothea::dorothea_mm
  }else{
    net <- dorothea::dorothea_hs
  }
  if(confidence == "recommend"){
    net2 <- net %>% filter(confidence != "D") %>% filter(confidence != "E")
  }else net2 <- net
  if(type == "DoRothEA regulon (activator)") net2 <- net2%>% filter(mor == 1)
  if(type == "DoRothEA regulon (repressor)") net2 <- net2%>% filter(mor == -1)
  my.symbols <- gsub("\\..*","", net2$target)
  gene_IDs<-AnnotationDbi::select(org(species),keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("SYMBOL", "ENTREZID"))
  colnames(gene_IDs) <- c("target", "ENTREZID")
  gene_IDs <- gene_IDs %>% distinct(target, .keep_all = T)
  gene_IDs <- na.omit(gene_IDs)
  net2 <- merge(net2, gene_IDs, by="target")
  net3 <- data.frame(gs_name = net2$tf, entrez_gene = net2$ENTREZID, target = net2$target, confidence = net2$confidence)
  net3 <- dplyr::arrange(net3, gs_name)
  return(net3)
}
cnet_for_output <- function(data, plot_data, Gene_set, Species){
  if(!is.null(Gene_set)){
    if(is.null(data)){
      return(NULL)
    }else{
      if(Species != "not selected"){
        withProgress(message = "cnet plot",{
          p <- plot_data
          print(p)
          incProgress(1)
        })
      }else return(NULL)
    }
  }
}
range_changer <- function(data){
  chr<-gsub("\\:.+$", "", rownames(data))
  start <- gsub(".+\\:","",gsub("\\-.+$", "", rownames(data)) )
  end <- gsub(".+\\-","",rownames(data))
  data$chr <- chr
  data$start <- as.numeric(start)
  data$end <- as.numeric(end)
  return(data)
}
enrich_for_table <- function(data, H_t2g, Gene_set){
  if(length(as.data.frame(data)$Description) == 0 || is.null(H_t2g)){
    return(NULL)
  }else{
    colnames(data)[1] <- "gs_name"
    H_t2g <- H_t2g %>% distinct(gs_name, .keep_all = T)
    data2 <- left_join(data, H_t2g, by="gs_name")  %>% as.data.frame()
    if(Gene_set == "DoRothEA regulon (activator)" || Gene_set == "DoRothEA regulon (repressor)"){
      data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name, Confidence = data2$confidence,
                          Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                          p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
    }else{
      if(Gene_set == "Custom gene set"){
        data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name,
                            Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                            p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
      }else{
        data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name, ID = data2$gs_id, Description = data2$gs_description,
                            Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                            p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
        
      }
      return(data3) 
    }
  }
}

read_known_results<-function (path, homer_dir = TRUE) {
  if (homer_dir == TRUE) {
    path <- paste0(path, "/knownResults.txt")
  }
  if (!file.exists(path)) {
    warning(paste("File", path, "does not exist"))
    return(NULL)
  }
  col_spec <- readr::cols("c", "c", "d", "-", "d", "d", "c", 
                          "d", "c")
  raw <- readr::read_tsv(path, col_types = col_spec)
  colnames(raw) <- c("motif_name", "consensus", "log_p_value", 
                     "fdr", "tgt_num", "tgt_pct", "bgd_num", "bgd_pct")
  tmp <- raw %>% tidyr::separate_("motif_name", c("motif_name", 
                                                  "experiment", "database"), "/", extra = "drop")
  parsed <- .parse_homer_subfields(tmp) %>% dplyr::mutate_at(vars(contains("pct")), 
                                                             .parse_pcts)
  known <- .append_known_pwm(parsed)
  names(known$motif_pwm) <- known$motif_name
  return(known)
}
.parse_homer_subfields <- function(motif_tbl) {
  cond <- stringr::str_detect(motif_tbl$motif_name, "/") %>%
    sum(., na.rm = TRUE) > 0
  if (cond == TRUE) {
    motif_tbl <- motif_tbl %>%
      tidyr::separate_('motif_name',
                       c('motif_name', 'experiment', 'database'),
                       '/', extra = "drop", fill = "right")
  }
  
  ## Detect if parentheses are present in motif_name
  ## to break apart into motif_name vs. motif_family
  cond <- stringr::str_detect(motif_tbl$motif_name, '\\(') %>%
    sum(., na.rm = TRUE) > 0
  if (cond == TRUE) {
    motif_tbl <- motif_tbl %>%
      tidyr::separate_('motif_name',
                       c('motif_name', 'motif_family'),
                       '\\(', extra = "drop", fill = "right")
    motif_tbl$motif_family <- stringr::str_replace(motif_tbl$motif_family, '\\)', '')
  }
  
  ## Detect If parentheses are present in experiment
  ## to break apart into experiment vs. accession
  if ("experiment" %in% colnames(motif_tbl)) {
    cond <- stringr::str_detect(motif_tbl$experiment, '\\(') %>%
      sum(., na.rm = TRUE) > 0
    if (cond == TRUE) {
      motif_tbl <- motif_tbl %>%
        tidyr::separate_('experiment',
                         c('experiment', 'accession'),
                         '\\(', extra = "drop", fill = "right")
      motif_tbl$accession <- stringr::str_replace(motif_tbl$accession, '\\)', '')
    }
  }
  
  return(motif_tbl)
}
.parse_pcts <- function(x) {
  stringr::str_replace(x, "%", "") %>%
    as.numeric() * 0.01
}
.append_known_pwm <- function(known_results) {
  data(HOMER_motifs, envir = environment())
  hm <- HOMER_motifs %>%
    ##        dplyr::filter(rlang::UQ(rlang::sym('log_odds_detection')) > 0) %>%
    dplyr::select("motif_name", "motif_family", "experiment", "accession",
                  "motif_pwm", "log_odds_detection")
  dplyr::inner_join(known_results, hm, 
                    by = c("motif_name", "motif_family", "experiment", "accession"))
}   
findMotif <- function(df,anno_data = NULL,Species,type = "Genome-wide",motif,size,back="random",bw_count=NULL,other_data=NULL){
    switch(Species,
           "Mus musculus (mm10)" = ref <- "mm10",
           "Homo sapiens (hg19)" = ref <- "hg19",
           "Homo sapiens (hg38)" = ref <- "hg38")
    switch(motif,
           "known motif" = time <- "10 ~ 20",
           "known and de novo motifs" = time <- "20 ~ 30")
  withProgress(message = paste0("Motif analysis takes about ",time," min per group"),{
    if(type == "Genome-wide" || type == "Other") {
      group_name <- names(df)
      group_file <- length(names(df))
    }else{
      group_name <- unique(df$group)
      group_file <- length(unique(df$group))
    }
    perc <- 0
    df2 <- list()
    path <- format(Sys.time(), "%Y%m%d_%H%M_homer")
    dir.create(path = path)
    print(group_name)
    for(name in group_name){
      group_dir <- paste0(path, "/",name)
      dir.create(path = group_dir)
      perc <- perc + 1
      if(!is.null(anno_data)){
        if(type == "Genome-wide"){
          data <- df[[name]]
          data2 <- anno_data %>% dplyr::filter(locus %in% rownames(data))
          y <- with(data2, GRanges(seqnames = seqnames, 
                                   ranges = IRanges(start,end)))
        }else{
          data <- dplyr::filter(df, group == name)
          y <- subset(anno_data, gene_id %in% data$ENTREZID)
          y <- as.data.frame(y)
          y <- with(y, GRanges(seqnames = seqnames, 
                                   ranges = IRanges(start,end)))
        }
      }else y<- df[[name]]
      y <- as.data.frame(y)
      print(head(y))
      if(dim(y)[1] != 0){
        if(type=="Genome-wide"){
        if(back == "random"){
          bg <-'automatic'
        }else{
          data <- anno_data %>% dplyr::filter(! locus %in% rownames(data))
          data2 <- range_changer(data)
          bg <- data.frame(seqnames = data2$chr,start=data2$start,end=data2$end)
        }
          }
        if(type== "Promoter") {
          bg <- subset(anno_data, ! gene_id %in% data$ENTREZID) 
          bg <- as.data.frame(bg)
          bg <- with(bg, GRanges(seqnames = seqnames, 
                              ranges = IRanges(start,end)))
          bg <- as.data.frame(bg)
        }
        if(type=="Other"){
          if(back == "random"){
            bg <-'automatic'
          }else{
            z <- with(y, GRanges(seqnames = seqnames, 
                                 ranges = IRanges(start,end)))
            bg <- exclude_bed(other_data,z)
            bg <- as.data.frame(bg)
          }
            }
        print(head(bg))
      switch(motif,
             "known motif" = motif_type <- TRUE,
             "known and de novo motifs" = motif_type <- FALSE)
      find_motifs_genome(
        y,
        path = group_dir,
        genome = ref, 
        motif_length = c(8,10,12),
        scan_size = size,
        optimize_count = 25,
        background = bg,
        local_background = FALSE,
        only_known = motif_type, only_denovo = FALSE,
        cores = 2, cache = 500,
        overwrite = TRUE, keep_minimal = FALSE
      )
      df2[[name]] <- group_dir
      }
      incProgress(1/group_file, message = paste("Finish motif analysis of Group '", name, "', ", perc, "/", group_file,sep = ""))
    }
    return(df2)
  })
}

known_motif <- function(df){
  df2 <- data.frame(matrix(rep(NA, 14), nrow=1))[numeric(0), ]
  colnames(df2) <- c("motif_name", "motif_family", "experiment", "accession", "database", "consensus", "p_value",
                     "fdr", "tgt_num", "tgt_pct", "bgd_num", "bgd_pct", "log_odds_detection","Group")
  for(name in names(df)){
    if(file.exists(paste0(df[[name]],"/knownResults.txt")) == TRUE){
    known <-  as.data.frame(read_known_results(path = df[[name]]))[,-13]
    known$Group <- name
    df2 <- rbind(df2,known)
    }
  }
  colnames(df2) <- c("motif_name", "motif_family", "experiment", "accession", "database", "consensus", "p_value",
                     "fdr", "tgt_num", "tgt_pct", "bgd_num", "bgd_pct", "log_odds_detection","Group")
  return(df2)
}

denovo_motif <- function(df){
  df2 <- data.frame(matrix(rep(NA, 18), nrow=1))[numeric(0), ]
  colnames(df2) <- c("consensus","motif_name", "log_odds_detection","motif_id",  "log_p_value_detection", 
                     "tgt_num", "tgt_pct", "bgd_num", "bgd_pct"," log_p_value","fdr","tgt_pos","tgt_std","bgd_pos","bgd_std","strand_bias","multiplicity", "Group")
  for(name in names(df)){
    if(file.exists(paste0(df[[name]],"/homerResults.html")) == TRUE){
      known <-  as.data.frame(read_denovo_results(path = df[[name]]))[,-5]
      known$Group <- name
      df2 <- rbind(df2,known)
    }
  }
  return(df2)
}


homer_Motifplot <- function(df, showCategory=5){
  df2 <- data.frame(matrix(rep(NA, 15), nrow=1))[numeric(0), ]

  for(name in names(df)){
    print(file.exists(paste0(df[[name]],"/knownResults.txt")))
    
    if(file.exists(paste0(df[[name]],"/knownResults.txt")) == TRUE){
    res <- as.data.frame(as.data.frame(read_known_results(path = df[[name]])))
    if(length(res$motif_name) != 0){
    res$Group <- name
    if(length(rownames(res)) > showCategory){
      res <- res[1:showCategory,]
    }
    df2 <- rbind(df2, res)
    }
    }
  }
  if(length(df2$motif_name) == 0){
    return(NULL)
  }else{
    colnames(df2) <- c("motif_name", "motif_family", "experiment", "accession", "database", "consensus", "p_value",
                       "fdr", "tgt_num", "tgt_pct", "bgd_num", "bgd_pct","motif_pwm","log_odds_detection","Group")
    df2 <- dplyr::mutate(df2, x = paste0(Group, 1/-log10(eval(parse(text = "p_value")))))
    df2$x <- gsub(":","", df2$x)
    df2 <- dplyr::arrange(df2, x)
    idx <- order(df2[["x"]], decreasing = FALSE)
    df2$motif_name <- factor(df2$motif_name,
                             levels=rev(unique(df2$motif_name[idx])))
    d <- ggplot(df2, aes(x = Group,y= motif_name,color=p_value,size=log_odds_detection))+
      geom_point() +
      scale_color_continuous(low="red", high="blue",
                             guide=guide_colorbar(reverse=TRUE)) +
      scale_size(range=c(1, 6))+ theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) +
      scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top")+
      theme(plot.margin=margin(l=-0.75,unit="cm"))
    df2 <- df2 %>% distinct(motif_name, .keep_all = T)
    pfm <- df2$motif_pwm
    pfm1 <- list()
    for(name in names(pfm)){
      pfm1[[name]] <- t(pfm[[name]])
    }
    Seqlogo <- as.grob(ggseqlogo(pfm1,ncol = 1)+
                         theme_void()) 
    p <- plot_grid(plot_grid(NULL, Seqlogo, ncol = 1, rel_heights = c(0.05:10)),as.grob(d))
    
    return(p)
  }
}

files2GRangelist <- function(files){
  Glist <- GRangesList()
  for(name in names(files)){
    extension <- gsub(".+\\.","",files$name)
    file_name <- gsub(".+\\/","",files$name)
    file_name<- gsub("\\..+$", "", file_name)
    print(file_name)
    print(files$name)
    input <- toGRanges(files$name, format=extension, header=FALSE) 
    Glist[[file_name]] <- input 
  }
  return(Glist)
}

integrate_ChIP_RNA <- function (result_geneRP, result_geneDiff, lfc_threshold = 1, 
                                padj_threshold = 0.05) {
  if ("GRanges" %in% class(result_geneRP)) {
    stop("sorry, please use the the simplify result or metadata(fullRP_hit)$peakRP_gene", 
         call. = FALSE)
  }
  merge_result <- dplyr::left_join(result_geneRP, result_geneDiff, 
                                   by = "gene_id")
  allGenes_N <- as.double(nrow(merge_result))
  merge_result <- merge_result %>% dplyr::mutate(diff_rank = rank(padj, 
                                                                  na.last = "keep"), diff_rank = dplyr::case_when(is.na(diff_rank) ~ 
                                                                                                                    allGenes_N, TRUE ~ diff_rank), rankProduct = RP_rank * 
                                                   diff_rank, rankOf_rankProduct = rank(rankProduct)) %>% 
    dplyr::arrange(rankOf_rankProduct) %>% dplyr::mutate(gene_category = dplyr::case_when(log2FoldChange > 
                                                                                            lfc_threshold & padj < padj_threshold ~ "up", log2FoldChange < 
                                                                                            -lfc_threshold & padj < padj_threshold ~ "down", TRUE ~ 
                                                                                            "NS"), gene_category = factor(gene_category, levels = c("up", 
                                                                                                                                                        "down", "NS")))
  upGenes_rank <- filter(merge_result, gene_category == "up")$RP_rank
  downGenes_rank <- filter(merge_result, gene_category == "down")$RP_rank
  staticGenes_rank <- filter(merge_result, gene_category == 
                               "NS")$RP_rank
  if (length(upGenes_rank) == 0 & length(downGenes_rank) == 
      0) {
    warning("no significant genes, just returing rank product result", 
            call. = FALSE)
    return(merge_result)
  }
  else if (length(upGenes_rank) == 0) {
    warning("no significant up genes, just returing rank product result", 
            call. = FALSE)
    return(merge_result)
  }
  else if (length(downGenes_rank) == 0) {
    warning("no significant down genes, just returing rank product result", 
            call. = FALSE)
    return(merge_result)
  }
  up_static_pvalue <- suppressWarnings(ks.test(upGenes_rank, 
                                               staticGenes_rank)$p.value)
  down_static_pvalue <- suppressWarnings(ks.test(downGenes_rank, 
                                                 staticGenes_rank)$p.value)
  ks_test <- paste0("\n Kolmogorov-Smirnov Tests ", "\n pvalue of up vs NS: ", 
                    format(up_static_pvalue, digits = 3, scientific = TRUE), 
                    "\n pvalue of down vs NS: ", format(down_static_pvalue, 
                                                            digits = 3, scientific = TRUE))
  annotate_df <- data.frame(xpos = -Inf, ypos = Inf, annotateText = ks_test, 
                            hjustvar = 0, vjustvar = 1)
  p <- merge_result %>% ggplot2::ggplot(aes(x = RP_rank)) + 
    ggplot2::stat_ecdf(aes(color = gene_category), geom = "line") + 
    ggplot2::geom_text(data = annotate_df, aes(x = xpos, 
                                               y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText)) + 
    ggplot2::xlab("Regulatory potential rank") + ggplot2::ylab("Cumulative Probability")+
    ggplot2::scale_x_continuous(expand = c(0,0)) + ggplot2::theme_bw(base_size = 15) + 
    ggplot2::scale_color_manual(breaks = c("up","down","NS"),values = c("#F8766D","#00BFC4","grey"))
  return(p)
}

enrich_viewer_forMulti2 <- function(data3, Species, Gene_set, org, org_code, H_t2g){
  if(!is.null(Gene_set)){
    if(is.null(data3)){
      return(NULL)
    }else{
      withProgress(message = "enrichment analysis",{
        if(is.null(H_t2g)){
          df <- NULL
        }else{
          H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in unique(data3$Group)) {
            sum <- length(data3$ENTREZID[data3$Group == name])
            em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(setReadable(em, org, 'ENTREZID'))
                cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
                df <- rbind(df, cnet1)
              }
            }
          }
        }
        if(length(df$ID) !=0){
          df$GeneRatio <- parse_ratio(df$GeneRatio)
          return(df)
        }else return(NULL)
      })
    }
  } 
}

enrich_gene_list <- function(data, Gene_set, H_t2g, org){
  if(!is.null(Gene_set)){
    if(is.null(data)){
      return(NULL)
    }else{
      if(is.null(H_t2g)){
        df <- NULL
      }else{
        H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
        df <- list()
        for (name in unique(data$Group)) {
          sum <- length(data$ENTREZID[data$Group == name])
          em <- enricher(data$ENTREZID[data$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
          if (length(as.data.frame(em)$ID) != 0) {
            if(length(colnames(as.data.frame(em))) == 9){
              cnet1 <- setReadable(em, org, 'ENTREZID')
              df[[name]] <- cnet1
            }
          }
        }
      }
      return(df)
    }
  }
}

enrich_genelist <- function(data, enrich_gene_list, showCategory=5){
  if(is.null(data) || is.null(enrich_gene_list)){
    return(NULL)
  }else{
    df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
    colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
    for (name in names(enrich_gene_list)) {
      sum <- length(data$ENTREZID[data$Group == name])
      em <- enrich_gene_list[[name]]
      if (length(as.data.frame(em)$ID) != 0) {
        if(length(colnames(as.data.frame(em))) == 9){
          cnet1 <- as.data.frame(em)
          cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
          cnet1 <- cnet1[sort(cnet1$pvalue, decreasing = F, index=T)$ix,]
          if (length(cnet1$pvalue) > showCategory){
            cnet1 <- cnet1[1:showCategory,]
          }
          df <- rbind(df, cnet1)
        }
      }
    }
    if ((length(df$Description) == 0) || length(which(!is.na(unique(df$qvalue)))) == 0) {
      p1 <- NULL
    } else{
      df$GeneRatio <- parse_ratio(df$GeneRatio)
      df <- dplyr::filter(df, !is.na(qvalue))
      df$Description <- gsub("_", " ", df$Description)
      df <- dplyr::mutate(df, x = paste0(Group, 1/(-log10(eval(parse(text = "qvalue"))))))
      df$x <- gsub(":","", df$x)
      df <- dplyr::arrange(df, x)
      idx <- order(df[["x"]], decreasing = FALSE)
      df$Description <- factor(df$Description,
                               levels=rev(unique(df$Description[idx])))
      p1 <- as.grob(ggplot(df, aes(x = Group,y= Description,color=qvalue,size=GeneRatio))+
                      geom_point() +
                      scale_color_continuous(low="red", high="blue",
                                             guide=guide_colorbar(reverse=TRUE)) +
                      scale_size(range=c(1, 6))+ theme_dose(font.size=12)+ylab(NULL)+xlab(NULL)+
                      scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
      p <- plot_grid(p1)
      return(p)
    }
  }
}
enrich_for_table <- function(data, H_t2g, Gene_set){
  if(length(as.data.frame(data)$Description) == 0 || is.null(H_t2g)){
    return(NULL)
  }else{
    colnames(data)[1] <- "gs_name"
    H_t2g <- H_t2g %>% distinct(gs_name, .keep_all = T)
    data2 <- left_join(data, H_t2g, by="gs_name")  %>% as.data.frame()
    if(Gene_set == "DoRothEA regulon (activator)" || Gene_set == "DoRothEA regulon (repressor)"){
      data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name, Confidence = data2$confidence,
                          Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                          p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
    }else{
      if(Gene_set == "Custom gene set"){
        data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name,
                            Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                            p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
      }else{
        data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name, ID = data2$gs_id, Description = data2$gs_description,
                            Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                            p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
        
      }
      return(data3) 
    }
  }
}
symbol2gene_id <- function(data,org){
  my.symbols <- rownames(data)
  gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("SYMBOL","ENTREZID"))
  colnames(gene_IDs) <- c("SYMBOL","gene_id")
  gene_IDs <- gene_IDs %>% distinct(SYMBOL, .keep_all = T) %>% na.omit()
  gene_IDs <- data.frame(gene_id = gene_IDs$gene_id, row.names = gene_IDs$SYMBOL)
  data <- merge(gene_IDs,data,by=0)
  data <- data[,-1]
  return(data)
}

data_trac <- function(y,gene_position,gen,txdb,org,filetype=NULL,bw_files,
                      bam_files,track_additional_files){
  chr <- gene_position$seqnames
  grtrack <- GeneRegionTrack(txdb,
                             chromosome = chr, name = "UCSC known genes",geneSymbol = TRUE,
                             transcriptAnnotation = "symbol",genome = gen,
                             background.title = "grey",cex = 1.25)
  symbols <- unlist(mapIds(org, gene(grtrack), "SYMBOL", "ENTREZID", multiVals = "first"))
  symbol(grtrack) <- symbols[gene(grtrack)]
  displayPars(grtrack) <- list(fontsize = 15)
  if(!is.null(filetype)){
  switch(filetype,
         "Row1" = bw_files <- bw_files,
         "Row2" = bw_files <- bam_files)
  }
  if(!is.null(track_additional_files)) bw_files <- c(bw_files, track_additional_files)
  df <- list()
  c <- 0
  unique <- unique(gsub("\\_.+$", "", names(bw_files)))
  num <- length(unique)
  col <- list()
  col_count = num + 1
  for(name in unique){
    col_count <- col_count - 1
    col[[name]] <- rainbow_hcl(num,c = 100)[col_count]
  }
  for(name in names(bw_files)){
    c <- c + 1
    if(!is.null(filetype)){
    if(filetype == "Row2"){ 
      bai <- file.exists(paste0(bw_files[[name]],".bai"))
      if(bai == FALSE) a <- indexBam(bw_files[[name]])
    }}
    name2 <- gsub("\\_.+$", "", name)
    df[[name]] <- DataTrack(range = bw_files[[name]], type = "l",genome = gen,
                            name = gsub("\\..+$", "", name), window = -1,
                            chromosome = chr, background.title = col[[name2]],cex = 0.8,
                            col.histogram = col[[name2]], cex.axis=0.8,cex.main=0.8, cex.title = 0.8,
                            fill.histogram = col[[name2]])
  }
  df[["grtrack"]] <- grtrack
  return(df)
}

RNAseqDEGimport <- function(tmp,exampleButton){
  withProgress(message = "Importing an RNA-seq DEG regult file, please wait",{
    if(is.null(tmp) && exampleButton > 0 )  tmp = "data/RNAseq.txt"
    if(is.null(tmp)) {
      return(NULL)
    }else{
      if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1,quote = "")
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = "")
      rownames(df) = gsub("\"", "", rownames(df))
      if(length(colnames(df)) != 0){
        if(str_detect(colnames(df)[1], "^X\\.")){
          colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
        }
      }
      return(df)
    }
  })
}
RNAseqDEG_ann <- function(RNAdata,org){
  RNAdata$log2FoldChange <- -RNAdata$log2FoldChange
  if(str_detect(rownames(RNAdata)[1], "ENS")){
    my.symbols <- gsub("\\..*","", rownames(RNAdata))
    gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                    keytype = "ENSEMBL",
                                    columns = c("ENSEMBL","SYMBOL","ENTREZID"))
    colnames(gene_IDs) <- c("EnsemblID","Symbol","gene_id")
    RNAdata$EnsemblID <- gsub("\\..*","", rownames(RNAdata))
    gene_IDs <- gene_IDs %>% distinct(EnsemblID, .keep_all = T)
  }else{
    my.symbols <- rownames(RNAdata)
    gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                    keytype = "SYMBOL",
                                    columns = c("SYMBOL", "ENTREZID"))
    colnames(gene_IDs) <- c("Symbol", "gene_id")
    gene_IDs <- gene_IDs %>% distinct(Symbol, .keep_all = T)
    RNAdata$Symbol <- rownames(RNAdata) 
  }
  data <- merge(RNAdata, gene_IDs, by="Symbol")
  return(data)
}
mmAnno <- function(peak,genomic_region=NULL,txdb,peak_distance){
  if(!is.null(peak)){
    if(!is.null(genomic_region)){
    if(genomic_region == "Promoter"){ 
      peak <- peak  %>% as.data.frame() %>% distinct(start, .keep_all = T)
      peak <- with(peak, GRanges(seqnames = seqnames,ranges = IRanges(start=start,end=end)))
    }}
    range <- peak_distance * 1000
    mcols(peak) <- DataFrame(feature_id = paste0("peak_", seq_len(length(peak))))
    if(length(as.data.frame(peak)$seqnames) != 0){
      mmAnno_up <- mm_geneScan(peak, txdb,upstream = range,downstream = range)
    }else mmAnno_up <- NULL
    return(mmAnno_up)
  }else return(NULL)
}
RP_f <- function(mmAnno,txdb){
  if(!is.null(mmAnno)) {
    result_geneRP_up <- calcRP_TFHit(mmAnno = mmAnno,Txdb = txdb)
  }else result_geneRP_up <- NULL
  if(!is.null(result_geneRP_up)) {
    result_geneRP <-result_geneRP_up
    result_geneRP <- result_geneRP %>% dplyr::arrange(-sumRP)
    result_geneRP$RP_rank <- rownames(result_geneRP) %>% as.numeric()
    return(result_geneRP)
  }else return(NULL)
}
regulatory_potential_f <- function(species,data,result_geneRP,DEG_fc,DEG_fdr){
  if(species != "not selected"){
    merge_data <- integrate_ChIP_RNA(
      result_geneRP = result_geneRP,
      result_geneDiff = data,lfc_threshold = log(DEG_fc,2),padj_threshold = DEG_fdr
    )
    return(merge_data)
  }
}

GetGRanges <- function (LoadFile, simple = FALSE,sepr = "\t", simplify = FALSE) {
    if (sum(class(LoadFile) == "character")) {
      RangesTable <- read.delim(LoadFile, sep = sepr, header = FALSE, 
                                comment.char = "#")
      if(str_detect(RangesTable[1,2], "tart") == TRUE) {
        RangesTable <- RangesTable[-1]
      }else RangesTable <- RangesTable
    }
    Chromosomes <- as.vector(RangesTable[, 1])
    Start <- as.numeric(as.vector(RangesTable[, 2]))
    End <- as.numeric(as.vector(RangesTable[, 3]))
    RegionRanges <- GRanges(seqnames = Chromosomes, ranges = IRanges(start = Start, 
                                                                     end = End))
    if (simple == FALSE) {
      if (ncol(RangesTable) > 4) {
        ID <- as.vector(RangesTable[, 4])
        Score <- as.vector(RangesTable[, 5])
        if (ncol(RangesTable) > 6) {
          Strand <- rep("*", nrow(RangesTable))
          RemainderColumn <- as.data.frame(RangesTable[, 
                                                       -c(1:6)])
          elementMetadata(RegionRanges) <- cbind(ID, 
                                                 Score, Strand, RemainderColumn)
        }
        else {
          elementMetadata(RegionRanges) <- cbind(ID, 
                                                 Score)
        }
      }
    }
  return(RegionRanges)
}
