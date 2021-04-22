#*********************************
#*Script for RPPA data analysis
#*
#*MECHANICAL STRESS IN PANCREATIC CANCER: 
#*SIGNALING PATHWAY ADAPTATION ACTIVATES CYTOSKELETAL REMODELING 
#*AND ENHANCES CELL MIGRATION
#*
#*2021-03-30
#*********************************

###Clear environment and set options
rm(list=ls())
options(stringsAsFactors = FALSE)

{
  # library(d3heatmap)
  library(gplots)
  library(dplyr)
  library(ggplot2)
  library(hrbrthemes)
  library(Hmisc)
  library(rafalib)
  library(RColorBrewer)
  
}

###Load RPPA data
{
  data_unprocessed <- read.table("Kalli_compressed_uncompressed_monoculture.xls",sep="\t",header=T) #
  data_processed <- data_unprocessed[,2:ncol(data_unprocessed)] %>% as.matrix %>% na.omit()
  data_processed<-t(data_processed)
  protein_list <- data_unprocessed[1:ncol(data_processed),1]
  colnames(data_processed)=protein_list
  rownames(data_processed) <- gsub('_MIA.PaCa2','',rownames(data_processed))
}

###Identify samples in dataset
{
  compressed_16<-c(2,4,6)
  uncompressed_16<-c(1,3,5)
}

###Heatmap for Figure 1
{
  #setup heatmap palette
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)
  
  #calculate fold changes of each protein
  #16h
  fc_16h<-c()
  for(i in 1:ncol(data_processed)){
    fc_16h[i]<-(sum(data_processed[compressed_16,i])-sum(data_processed[uncompressed_16,i]))/sum(data_processed[uncompressed_16,i])
  }
  names(fc_16h)<-colnames(data_processed)
  
  ###adjust fold change to include all proteins in experiments
  experiment_proteinlist<-read.csv("experiment_proteinlist.csv",header = TRUE)
  experiment_proteinlist_idx<-which(colnames(data_processed) %in% experiment_proteinlist$RPPA_protein_name)
  
  #set fold change threshold
  foldchange_threshold<-0.2 #e.g. set 20% as 0.2
  
  idx_fc_up<-which(fc_16h>foldchange_threshold)
  idx_fc_down<-which(fc_16h<(-foldchange_threshold))
  
  names_up<-names(fc_16h)[idx_fc_up]
  names_down<-names(fc_16h)[idx_fc_down]
  
  proteinlist_fc_heatmap_16<-c(names_up,names_down)
  data_fc_heatmap_16<-data_processed[c(uncompressed_16,compressed_16),proteinlist_fc_heatmap_16]
  
  fc_appendix<-data.frame(proteinlist_fc_heatmap_16,fc_16h[c(names_up,names_down)])
  #write fc>20% gene names
  write.csv(colnames(data_fc_heatmap_16),"names_foldchange_appdixI.csv")
  write.csv(fc_appendix,"foldchange_appendixI.csv")
  
  #Heatmap with 3 rows of fold change, compressed/control average
  hm_fc_order_16h<-order(fc_16h[proteinlist_fc_heatmap_16],decreasing = TRUE)
  avg_control<-colSums(data_fc_heatmap_16[1:3,hm_fc_order_16h])/3
  
  avg_fc_data<-cbind(
    (data_fc_heatmap_16[4,hm_fc_order_16h]/avg_control),
    (data_fc_heatmap_16[5,hm_fc_order_16h]/avg_control),
    (data_fc_heatmap_16[6,hm_fc_order_16h]/avg_control)
  )
  colnames(avg_fc_data)<-c("Comp.#3","Comp.#2","Comp.#1")
  ylim1=0;
  ylim2=2;
  heatmap21<-heatmap.2(t(avg_fc_data),Rowv = FALSE,Colv=FALSE, 
                       main="16h 20% fold change", trace="none", margins= c(14,13), col=my_palette, breaks=seq(ylim1,ylim2,length.out=101),key.par=list(cex=0.5),keysize =0.5,cexRow=1.5,cexCol=1.2,
                       lhei=c(2,4), lwid=c(2,12),offsetRow =0,offsetCol =0,
                       adjCol = c(1,0.5))
  
}


###Volcano plot for appendix
###use fold changes and p values for volcano plots
{# Paired t-test
  #16h
    pVal_16_paired<-c()
    for(i in 1:ncol(data_processed)){
      pVal_16_paired[i]<-t.test(data_processed[compressed_16,i],data_processed[uncompressed_16,i],paired=TRUE)$p.value
      
    }
    
    names(pVal_16_paired)=protein_list[1:425]
    
    genelist_16_paired<-c(names(pVal_16_paired[pVal_16_paired<0.05]))
    
  volcanodata_16<-cbind(fc_16h,pVal_16_paired)%>%as.data.frame
  colnames(volcanodata_16)<-c("fc","pval")
  
  volcanodata_16$delabel <- NA
  volcanodata_16$delabel[abs(volcanodata_16$fc)>0.2 & volcanodata_16$pval<0.05] <- rownames(volcanodata_16)[abs(volcanodata_16$fc)>0.2 & volcanodata_16$pval<0.05]
  
  volcanodata_16$trend_16<-replicate(length(fc_16h),"no")
  volcanodata_16$trend_16[volcanodata_16$fc>0.2 & volcanodata_16$pval<0.05]<-"up"
  volcanodata_16$trend_16[volcanodata_16$fc<(-0.2) & volcanodata_16$pval<0.05]<-"down"
  
  colnames(volcanodata_16)<-c("fc","pval","lab","trend")
  
  library(ggrepel)
  
  ggplot(volcanodata_16,
         aes(x=log2(fc+1), y=-log10(pval),label=lab,color=trend))+
    geom_point(cex=1)+
    theme_minimal()+
    geom_text_repel(cex=3) +
    scale_color_manual(values=c("blue", "black", "red"))+
    geom_vline(xintercept=c(log2(-0.2+1), log2(0.2+1)), col="grey60") +
    geom_hline(yintercept=-log10(0.05), col="grey60")+
    ggtitle(label="16h paired")
}

###Pathway heatmap for appendix
{ #data for heatmap
  data_unprocessed <- read.table("Kalli_compressed_uncompressed_monoculture.xls",sep="\t",header=T) #
  data_processed <- data_unprocessed[,2:ncol(data_unprocessed)] %>% as.matrix %>% na.omit()
  data_processed<-t(data_processed)
  colnames(data_processed)=protein_list
  rownames(data_processed) <- gsub('_MIA.PaCa2','',rownames(data_processed))
  
  compressed_16<-c(2,4,6)
  uncompressed_16<-c(1,3,5)
  alldata_16h_ordered_by_fc<-data_processed[c(uncompressed_16,compressed_16),order(fc_16h,decreasing = TRUE)]
}

{#3 row of foldchange pairs for Ras/MAPK pathway members
  proteinlist_rasmapk <- c("c-Jun_pS73","C-Raf_pS338","JNK_pT183_Y185","MAPK_pT202_Y204", "MEK1_p_S217_S221","p90RSK_pT573","Shc_pY317","YB1_pS102") #"p38_pT180_Y182",
  data_rasmapk<-data_processed[c(uncompressed_16,compressed_16),proteinlist_rasmapk]
  rasmapk_fc_order_16h<-order(fc_16h[proteinlist_rasmapk],decreasing = TRUE)
  rasmapk_fc_data<-cbind(
    (data_rasmapk[4,rasmapk_fc_order_16h]/data_rasmapk[1,rasmapk_fc_order_16h]),
    (data_rasmapk[5,rasmapk_fc_order_16h]/data_rasmapk[2,rasmapk_fc_order_16h]),
    (data_rasmapk[6,rasmapk_fc_order_16h]/data_rasmapk[3,rasmapk_fc_order_16h])
  )
  colnames(rasmapk_fc_data)<-c("rep#1","rep#2","rep#3")
  
  ylim1=0;
  ylim2=2;
  par(cex.main=0.90)
  heatmap31<-heatmap.2(t(rasmapk_fc_data),Rowv = FALSE,Colv=FALSE, 
                       main="Ras/MAPK pathway protein fold change", 
                       trace="none", margins= c(14,13), 
                       col=my_palette, breaks=seq(ylim1,ylim2,length.out=101),
                       key.par=list(cex=0.5),keysize =0.5,
                       cexRow=1.5,cexCol=1.2,
                       lhei=c(2,6), lwid=c(2,10),#colorkey size
                       adjCol = c(1,0.5))
  
}

{#3 row of foldchange pairs for tsc/mtor pathway members
  proteinlist_tscmtor<-c("4E-BP1_pS65","4E-BP1_pT37_T46","p70-S6K_pT389","mTOR_pS2448","Rictor_pT1135","S6_pS235_S236","S6_pS240_S244")
  data_tscmtor<-data_processed[c(uncompressed_16,compressed_16),proteinlist_tscmtor]
  tscmtor_fc_order_16h<-order(fc_16h[proteinlist_tscmtor],decreasing = TRUE)
  tscmtor_fc_data<-cbind(
    (data_tscmtor[4,tscmtor_fc_order_16h]/data_tscmtor[1,tscmtor_fc_order_16h]),
    (data_tscmtor[5,tscmtor_fc_order_16h]/data_tscmtor[2,tscmtor_fc_order_16h]),
    (data_tscmtor[6,tscmtor_fc_order_16h]/data_tscmtor[3,tscmtor_fc_order_16h])
  )
  colnames(tscmtor_fc_data)<-c("rep#1","rep#2","rep#3")
  
  ylim1=0;
  ylim2=2;
  par(cex.main=0.90)
  heatmap32<-heatmap.2(t(tscmtor_fc_data),Rowv = FALSE,Colv=FALSE, 
                       main="TSC/mTOR pathway protein fold change", 
                       trace="none", margins= c(14,13), 
                       col=my_palette, breaks=seq(ylim1,ylim2,length.out=101),
                       key.par=list(cex=0.5),keysize =0.5,
                       cexRow=1.5,cexCol=1.2,
                       lhei=c(2,6), lwid=c(2,10),#colorkey size
                       adjCol = c(1,0.5))
  
}


###Pathway analysis for Figure 1
###Create functions that will be used later for error bars
{
  ###Create function to calculate standard error because R doesn't have one
  se <- function(x, na.rm=FALSE) {
    if (na.rm) x <- na.omit(x)
    sqrt(var(x)/length(x))
  }
  
  ###Function to summarize data for multiple replicates using mean, SD, and SE
  ###data = data table to summarize, varname = variable to summarize, groupnames = variables used to group replicates
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE), 
        se = se(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
    return(data_sum)
  }
  
}

###Select proteins for each pathway
{
  pi3kakt_pathway_members <- c("INPP4b","PTEN","Akt_pS473","Akt_pT308","GSK-3a-b_pS21_S9","p27_pT157","p27_pT198","PRAS40_pT246") 
  #pi3kakt_pathway_members_select <- c("Akt_pS473","Akt_pT308","GSK-3a-b_pS21_S9","p27_pT157","p27_pT198","PRAS40_pT246") 
  rtk_pathway_members <- c("EGFR_pY1173","HER2_pY1248","HER3_pY1289","Shc_pY317","Src_pY416","Src_pY527") # EGFR_pY1068 missing
  tscmtor_pathway_members <- c("4E-BP1_pS65","4E-BP1_pT37_T46","p70-S6K_pT389","mTOR_pS2448","Rictor_pT1135","S6_pS235_S236","S6_pS240_S244") # 
  rasmapk_pathway_members <- c("c-Jun_pS73","C-Raf_pS338","JNK_pT183_Y185","MAPK_pT202_Y204", "MEK1_p_S217_S221","p90RSK_pT573","Shc_pY317","YB1_pS102") #"p38_pT180_Y182",
  #rasmapk_pathway_select_members <- c("C-Raf_pS338","MAPK_pT202_Y204", "MEK1_p_S217_S221","p90RSK_pT573","YB1_pS102")
  
  cellcycle_pathway_members <- c("Cyclin-B1","Cyclin-D1","Cyclin-E1","p27_pT157","PCNA")
  # cellcycle_pathway_members <- c("Rb","CDKN2A","Cyclin-B1","Cyclin-D1","Cyclin-E1","p27_pT157","PCNA","Cdc42","CDK1_pT14","cdc2_pY15")
  apo_pathway_members <- c("Bad_pS112","Bcl-xL","Bcl2","Mcl-1","XIAP","Bak","Bax","Bid","Bim","Puma","Smac","Caspase-3-cleaved","Caspase-8","Caspase-7-cleaved-")
  dnadamage_pathway_members <- c("53BP1","ATM","Chk1_pS345","Chk2_pT68","p53","Rad50","Rad51","XRCC1")
  ##hormoneA_pathway_members <- c("ER","PR") #- this was for oct18
  hormoneA_pathway_members <- c("ER-a","PR") #- this is for oct19
  hormoneB_pathway_members <- c("AR","INPP4b","GATA3","Bcl2")
  
}

###Get IDs for select proteins
{
  id_apo_pathway_match <- match(apo_pathway_members,colnames(data_processed), nomatch=0) 
  id_cellcycle_pathway_match <- match(cellcycle_pathway_members,colnames(data_processed), nomatch=0) 
  id_dnadamage_pathway_match <- match(dnadamage_pathway_members,colnames(data_processed), nomatch=0) 
  id_hormoneA_pathway_match <- match(hormoneA_pathway_members,colnames(data_processed), nomatch=0) 
  id_hormoneB_pathway_match <- match(hormoneB_pathway_members,colnames(data_processed), nomatch=0) 
  id_pi3kakt_pathway_match <- match(pi3kakt_pathway_members,colnames(data_processed), nomatch=0) 
  #id_pi3kakt_pathway_select_match <- match(pi3kakt_pathway_members_select,colnames(data_processed), nomatch=0) 
  id_rtk_pathway_match <- match(rtk_pathway_members,colnames(data_processed), nomatch=0) 
  id_tscmtor_pathway_match <- match(tscmtor_pathway_members,colnames(data_processed), nomatch=0) 
  id_rasmapk_pathway_match <- match(rasmapk_pathway_members,colnames(data_processed), nomatch=0) 
  #id_rasmapk_pathway_select_match <- match(rasmapk_pathway_select_members,colnames(data_processed), nomatch=0) 
  
  id_pathways_to_plot <- c(id_apo_pathway_match, id_cellcycle_pathway_match, id_dnadamage_pathway_match, id_hormoneA_pathway_match, id_hormoneB_pathway_match, id_pi3kakt_pathway_match, id_rasmapk_pathway_match, id_rtk_pathway_match, id_tscmtor_pathway_match)
  
}

###Set coefficients for proteins in pathways (positive is positive regulator of pathway, negative is negative regulator of pathway)
{
  apo_pathway_coeff <- c(-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1)
  cellcycle_pathway_coeff <- c(1,1,1,1,1)
  # cellcycle_pathway_coeff <- c(-1,-1,1,1,1,1,1,1,1,1)
  dnadamage_pathway_coeff <- c(1,1,1,1,1,1,1,1)
  hormoneA_pathway_coeff <- c(1,1)
  hormoneB_pathway_coeff <- c(1,1,1,1)
  pi3kakt_pathway_coeff <- c(-1,-1,1,1,1,1,1,1)
  rtk_pathway_coeff <- c(1,1,1,1,1,1)
  rasmapk_pathway_coeff <- c(1,1,1,1,1,1,1,1)
  tscmtor_pathway_coeff <- c(1,1,1,1,1,1,1)
  
}

#Compressed 16h
###Calculate pathway scores for MIAPaCa2
{
  ###Select data for each pathway
  apo_pathway_data_mono_lap_192 <- data_processed[compressed_16,id_apo_pathway_match]
  cellcycle_pathway_data_mono_lap_192 <- data_processed[compressed_16,id_cellcycle_pathway_match]
  dnadamage_pathway_data_mono_lap_192 <- data_processed[compressed_16,id_dnadamage_pathway_match]
  hormoneA_pathway_data_mono_lap_192 <- data_processed[compressed_16,id_hormoneA_pathway_match]
  hormoneB_pathway_data_mono_lap_192 <- data_processed[compressed_16,id_hormoneB_pathway_match]
  pi3kakt_pathway_data_mono_lap_192 <- data_processed[compressed_16,id_pi3kakt_pathway_match]
  rasmapk_pathway_data_mono_lap_192 <- data_processed[compressed_16,id_rasmapk_pathway_match]
  rtk_pathway_data_mono_lap_192 <- data_processed[compressed_16,id_rtk_pathway_match]
  tscmtor_pathway_data_mono_lap_192 <- data_processed[compressed_16,id_tscmtor_pathway_match]
  
  ###Multiply each protein's expression by its coefficient and sum all proteins in pathway
  apo_pathway_score_mono_lap_192 <- apply(t(apo_pathway_data_mono_lap_192)*apo_pathway_coeff,2,sum) 
  cellcycle_pathway_score_mono_lap_192 <- apply(t(cellcycle_pathway_data_mono_lap_192)*cellcycle_pathway_coeff,2,sum) 
  dnadamage_pathway_score_mono_lap_192 <- apply(t(dnadamage_pathway_data_mono_lap_192)*dnadamage_pathway_coeff,2,sum) 
  hormoneA_pathway_score_mono_lap_192 <- apply(t(hormoneA_pathway_data_mono_lap_192)*hormoneA_pathway_coeff,2,sum) 
  hormoneB_pathway_score_mono_lap_192 <- apply(t(hormoneB_pathway_data_mono_lap_192)*hormoneB_pathway_coeff,2,sum) 
  pi3kakt_pathway_score_mono_lap_192 <- apply(t(pi3kakt_pathway_data_mono_lap_192)*pi3kakt_pathway_coeff,2,sum) 
  rasmapk_pathway_score_mono_lap_192 <- apply(t(rasmapk_pathway_data_mono_lap_192)*rasmapk_pathway_coeff,2,sum) 
  rtk_pathway_score_mono_lap_192 <- apply(t(rtk_pathway_data_mono_lap_192)*rtk_pathway_coeff,2,sum)                 #sum columns of transposed matrix
  tscmtor_pathway_score_mono_lap_192 <- apply(t(tscmtor_pathway_data_mono_lap_192)*tscmtor_pathway_coeff,2,sum) 
  
  pathway_score_mono_lap_192 <- c(apo_pathway_score_mono_lap_192, cellcycle_pathway_score_mono_lap_192, dnadamage_pathway_score_mono_lap_192, hormoneA_pathway_score_mono_lap_192, hormoneB_pathway_score_mono_lap_192, pi3kakt_pathway_score_mono_lap_192, rasmapk_pathway_score_mono_lap_192, rtk_pathway_score_mono_lap_192, tscmtor_pathway_score_mono_lap_192)
  
}


#Uncompressed 16h
###Calculate pathway scores for MIAPaCa2 uncomp16
{
  ###Select data for each pathway
  apo_pathway_data_transar22_lap_192 <- data_processed[uncompressed_16,id_apo_pathway_match]
  cellcycle_pathway_data_transar22_lap_192 <- data_processed[uncompressed_16,id_cellcycle_pathway_match]
  dnadamage_pathway_data_transar22_lap_192 <- data_processed[uncompressed_16,id_dnadamage_pathway_match]
  hormoneA_pathway_data_transar22_lap_192 <- data_processed[uncompressed_16,id_hormoneA_pathway_match]
  hormoneB_pathway_data_transar22_lap_192 <- data_processed[uncompressed_16,id_hormoneB_pathway_match]
  pi3kakt_pathway_data_transar22_lap_192 <- data_processed[uncompressed_16,id_pi3kakt_pathway_match]
  rasmapk_pathway_data_transar22_lap_192 <- data_processed[uncompressed_16,id_rasmapk_pathway_match]
  rtk_pathway_data_transar22_lap_192 <- data_processed[uncompressed_16,id_rtk_pathway_match]
  tscmtor_pathway_data_transar22_lap_192 <- data_processed[uncompressed_16,id_tscmtor_pathway_match]
  
  ###Multiply each protein's expression by its coefficient and sum all proteins in pathway
  apo_pathway_score_transar22_lap_192 <- apply(t(apo_pathway_data_transar22_lap_192)*apo_pathway_coeff,2,sum) 
  cellcycle_pathway_score_transar22_lap_192 <- apply(t(cellcycle_pathway_data_transar22_lap_192)*cellcycle_pathway_coeff,2,sum) 
  dnadamage_pathway_score_transar22_lap_192 <- apply(t(dnadamage_pathway_data_transar22_lap_192)*dnadamage_pathway_coeff,2,sum) 
  hormoneA_pathway_score_transar22_lap_192 <- apply(t(hormoneA_pathway_data_transar22_lap_192)*hormoneA_pathway_coeff,2,sum) 
  hormoneB_pathway_score_transar22_lap_192 <- apply(t(hormoneB_pathway_data_transar22_lap_192)*hormoneB_pathway_coeff,2,sum) 
  pi3kakt_pathway_score_transar22_lap_192 <- apply(t(pi3kakt_pathway_data_transar22_lap_192)*pi3kakt_pathway_coeff,2,sum) 
  rasmapk_pathway_score_transar22_lap_192 <- apply(t(rasmapk_pathway_data_transar22_lap_192)*rasmapk_pathway_coeff,2,sum) 
  rtk_pathway_score_transar22_lap_192 <- apply(t(rtk_pathway_data_transar22_lap_192)*rtk_pathway_coeff,2,sum)                 #sum columns of transposed matrix
  tscmtor_pathway_score_transar22_lap_192 <- apply(t(tscmtor_pathway_data_transar22_lap_192)*tscmtor_pathway_coeff,2,sum) 
  
  pathway_score_transar22_lap_192 <- c(apo_pathway_score_transar22_lap_192, cellcycle_pathway_score_transar22_lap_192, dnadamage_pathway_score_transar22_lap_192, hormoneA_pathway_score_transar22_lap_192, hormoneB_pathway_score_transar22_lap_192, pi3kakt_pathway_score_transar22_lap_192, rasmapk_pathway_score_transar22_lap_192, rtk_pathway_score_transar22_lap_192, tscmtor_pathway_score_transar22_lap_192)
  
}
### Pathway score for only 16h samples (2 groups)
###Organize mia-paca data into data frame
{
  pathway <- c(replicate(3, "Apoptosis"), replicate(3, "Cell Cycle"), replicate(3, "DNA Damage"), replicate(3, "Hormone A"), replicate(3, "Hormone B"), replicate(3, "PI3KAKT"), replicate(3, "RASMAPK"), replicate(3, "RTK"), replicate(3, "TSCMTOR"))
  scores <- c( pathway_score_mono_lap_192, pathway_score_transar22_lap_192)
  cell_line <- c(replicate(54, "MIA PaCa2"))
  treatment <- c(replicate(27, "Compressed 16h"), replicate(27, "Uncompressed 16h"))
  data <- cbind.data.frame(cell_line, treatment, pathway, scores)
  
}

###Prepare mono_16h data for plotting
{
  ###Summarize the biological replicates in the MIAPaCa2 data using previously created function
  expression_stat <- data_summary(data, varname="scores", groupnames=c("cell_line", "pathway", "treatment"))
  colnames(expression_stat) <- c("cell_line", "pathway", "treatment", "scores", "sd", "se")
  
  ###Make data numeric
  class(expression_stat$se) <- "numeric"
  class(expression_stat$scores) <- "numeric"
  
  ###Set levels of treatment data (this sets the order of the bars in the plot)
  expression_stat$treatment <- factor(expression_stat$treatment, levels = c("Uncompressed 16h", "Compressed 16h"))
  
  
  ###Normalize to 1
  {
    new_stat<-c()
    new_error<-c()
    new_pval<-c()
    for(i in 1:9){
      new_stat[2*i]<-1
      new_stat[2*i-1]<-expression_stat$scores[2*i-1]/expression_stat$scores[2*i]
      
      new_error[2*i]<-expression_stat$se[2*i]/expression_stat$scores[2*i]
      new_error[2*i-1]<-expression_stat$se[2*i-1]/expression_stat$scores[2*i]
      
      new_pval[2*i]<-0
      
    }
    expression_stat$scores<-new_stat
    expression_stat$se<-new_error
  }
  
}

###Create mono_16h plot--this plot is not showed in the paper
{
  ggplot(expression_stat, aes(x=pathway, y=scores)) + 
    geom_col(mapping = aes(x = pathway, y = scores, fill = treatment), position = position_dodge()) +
    geom_errorbar(aes(ymin=scores-se, ymax=scores+se, color = treatment), width=.5, position=position_dodge(width = 1), show.legend = FALSE) +
    scale_fill_manual(values = c("skyblue", "palevioletred1")) +
    scale_color_manual(values = c("gray40",  "gray40")) +
    labs(x = "", 
         y = "Normalized score", fill = "Sample") +
    scale_x_discrete(labels=c("Apoptosis", "Cell Cycle", "DNA Damage", "Hormone A", "Hormone B", "PI3K/AKT", "Ras/MAPK","RTK", "TSC-mTOR")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5,size=12),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5, size = 13), axis.text = element_text(size=15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
    ggtitle("Pathway")
  
}

###barplot
{
  data_compressed<-data[1:27,]
  data_uncompressed<-data[28:54,]
  
  data_bar<-c()
  
  for (i in 1:9){
    data_bar[i]<-1+log2(mean(data_compressed$scores[c((3*i-2),(3*i-1),3*i)])/mean(data_uncompressed$scores[c((3*i-2),(3*i-1),3*i)]))
  }
  
  names(data_bar)<-c("Apoptosis", "Cell Cycle", "DNA Damage", "Hormone A", "Hormone B", "PI3K/AKT", "Ras/MAPK","RTK", "TSC-mTOR")
  
  data_bar<-data_bar[order(data_bar,decreasing = FALSE)]
  
  par(mar = c(7,2,2,2))
  barplot(data_bar,ylim=c(0,1.6),main="Pathway score difference",
          las=2)
  abline(h=0.8)
  abline(h=1.2)
}


