```R
## Load required packages
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(pROC)

## Load data
nreads <- seq(10000,50000,5000)
soupercell_result.ls <- lapply(1:length(nreads),function(i){
	df <- read.table(str_c("~/MitoSort/simulation/souporcell/",nreads[i],"/clusters.tsv"),header=TRUE,sep="\t")
	df
	})
Freemuxlet_result.ls <- lapply(1:length(nreads),function(i){
	df <- read.table(gzfile(str_c("~/MitoSort/simulation/Freemuxlet/",nreads[i],".clust1.samples.gz")),header=TRUE,sep="\t")
	df
	})
Vireo_result.ls <- lapply(1:length(nreads),function(i){
	df <- read.table(str_c("~/MitoSort/simulation/Vireo/",nreads[i],"/donor_ids.tsv"),header=TRUE,sep="\t")
	df
	})
MitoSort_result.ls <- lapply(1:length(nreads),function(i){
	df <- read.table(str_c("~/MitoSort/simulation/MitoSort/",nreads[i],"/MitoSort/Demultiplex_output/result_pvalue.txt"),header=TRUE,sep="\t")
	df
	})
true_labels.ls <- lapply(1:length(nreads),function(i){
	df <- read.table(str_c("~/MitoSort/simulation/bam_data/simulate_depth/100_",nreads[i],"_1_6_simulated.bam_barcodes.txt"),header=TRUE,sep="\t",comment.char="")
	df$type <- ifelse(grepl("\\+",df$old_read),"Doublet","Singlet")
	df
	})

## Figure2a, percentage of cell assignment
combined_df.ls <- lapply(1:length(nreads),function(i){
	soupercell_df <- soupercell_result.ls[[i]][,c("barcode","status","assignment")]
	soupercell_count_df <- soupercell_df %>%
		dplyr::group_by(status) %>%
		dplyr::summarise(counts=n()) %>%
		dplyr::mutate(percentage=counts/sum(counts)*100) %>%
		dplyr::rename(type=status)
	soupercell_count_df$method <- "Souporcell"
	soupercell_count_df$type[which(soupercell_count_df$type=="singlet")] <- "Singlet"
	soupercell_count_df$type[which(soupercell_count_df$type=="doublet")] <- "Doublet"
	soupercell_count_df$type[which(soupercell_count_df$type=="unassigned")] <- "Unassign"

	Freemuxlet_df <- Freemuxlet_result.ls[[i]][,c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS","DBL.BEST.GUESS")]
	Freemuxlet_count_df <- Freemuxlet_df %>% 
		dplyr::group_by(DROPLET.TYPE) %>% 
		dplyr::summarise(counts=n()) %>%
		dplyr::mutate(percentage=counts/sum(counts)*100) %>%
		dplyr::rename(type=DROPLET.TYPE)
	Freemuxlet_count_df$method <- "Freemuxlet"
	Freemuxlet_count_df$type[which(Freemuxlet_count_df$type=="SNG")] <- "Singlet"
	Freemuxlet_count_df$type[which(Freemuxlet_count_df$type=="DBL")] <- "Doublet"
	Freemuxlet_count_df$type[which(Freemuxlet_count_df$type=="AMB")] <- "Unassign"

	Vireo_df <- Vireo_result.ls[[i]][,c("cell","donor_id")]
	Vireo_df$type <- ifelse(Vireo_df$donor_id %in% str_c("donor",0:5),"Singlet",Vireo_df$donor_id)
	Vireo_count_df <- Vireo_df %>% 
		dplyr::group_by(type) %>% 
		dplyr::summarise(counts=n()) %>%
		dplyr::mutate(percentage=counts/sum(counts)*100)
	Vireo_count_df$method <- "Vireo"
	Vireo_count_df$type[which(Vireo_count_df$type=="doublet")] <- "Doublet"
	Vireo_count_df$type[which(Vireo_count_df$type=="unassigned")] <- "Unassign"


	MitoSort_df <- MitoSort_result.ls[[i]][,c("Barcode","Demultiplex")]
	MitoSort_df$type <- ifelse(MitoSort_df$Demultiplex=="Doublet","Doublet","Singlet")
	MitoSort_count_df <- MitoSort_df %>%
		dplyr::group_by(type) %>%
		dplyr::summarise(counts=n()) %>%
		dplyr::mutate(percentage=counts/sum(counts)*100)
	MitoSort_count_df$method <- "MitoSort"

	
	true_label_df <- true_labels.ls[[i]][,c("barcode","CB_barcode","type")]
	true_label_count_df <- true_label_df %>%
		dplyr::group_by(type) %>%
		dplyr::summarise(counts=n()) %>%
		dplyr::mutate(percentage=counts/sum(counts)*100)
	true_label_count_df$method <- "True label"

	df <- do.call(rbind,list(soupercell_count_df,Freemuxlet_count_df,Vireo_count_df,MitoSort_count_df,true_label_count_df))
	df$nreads <- nreads[i]
	df
	})
combined_df <- do.call(rbind,combined_df.ls)
combined_df$type <- factor(combined_df$type,levels=c("Singlet","Doublet","Unassign"))
combined_df$nreads <- factor(combined_df$nreads,levels=nreads)
combined_df$method <- factor(combined_df$method,levels=c("True label","MitoSort","Freemuxlet","Souporcell","Vireo"))

newpalette <- c(brewer.pal(9,"Blues")[6],brewer.pal(9,"Reds")[6],"lightgrey")
pdf("~/MitoSort/output/figure2/different_methods_cell_assignment_percentage.pdf",height=13,width=8)
ggplot(data=combined_df,aes(x=method,y=percentage,fill=type))+
	geom_bar(stat="identity",position="stack",width=0.7)+
	facet_wrap(. ~ nreads,ncol = 3)+
	theme_classic()+
	scale_fill_manual(values=newpalette,guide = guide_legend())+
	labs(x="",y="Percentage ( % )",fill="Classification")+
	theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.5)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5),angle = 45,hjust = 1),axis.title.x=element_text(color="black",size=rel(1.5)),axis.title.y=element_text(color="black",size=rel(1.3)),axis.line = element_line(colour="black",size = 0.8),plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm"),legend.title=element_text(color="black",size=rel(1.3)),legend.text=element_text(color="black",size=rel(1.2)),strip.background=element_rect(colour="black",fill="white",linetype="solid"),
        strip.text = element_text(face="bold",size=rel(1.1)),legend.position="bottom")
dev.off()


## Figure2b - runtime
MitoSort_runtime_df <- c()
for (i in 1:length(nreads)){
	df <- read.table(str_c("~/MitoSort/simulation/MitoSort/runtime/",nreads[i],"_runtime.txt"),header=FALSE,sep="\t")
	generate_matrix_time <- as.numeric(unlist(strsplit(df$V1[1],":"))[4])
	demultiplex_time <- as.numeric(unlist(strsplit(df$V1[2],":"))[4])
	total_time <- generate_matrix_time + demultiplex_time
	MitoSort_runtime_df <- data.frame(nreads=nreads[i],method="MitoSort",total_seconds=total_time) %>% rbind(MitoSort_runtime_df,.)
}

soupercell_runtime_df <- data.frame(nreads=nreads,method="Souporcell",total_seconds=c(509,578,606,674,720,725,710,737,771))
 
#2023/05/29_13:00:29-13:10-13:18:29
#2023/05/29_13:02:00-13:12-13:21:38
#2023/05/29_13:03:38-13:15-13:25:06
#2023/05/29_13:04:26-13:16-13:27:14
#2023/05/29_13:05:29-13:17-13:29:00
#2023/05/29_13:05:41-13:19-13:31:05
#2023/05/29_13:06:34-13:22-13:33:50
#2023/05/29_13:06:46-13:22-13:34:17
#2023/05/29_13:07:41-13:20-13:32:51

Freemuxlet_runtime_df <- c()
for (i in 1:length(nreads)){
	df <- read.table(str_c("~/MitoSort/simulation/Freemuxlet/runtime/",nreads[i],"_runtime.txt"),header=FALSE,sep="\t")
	Pileup_hours <- as.numeric(unlist(strsplit(df$V1[1],":"))[2])
	Pileup_minutes <- as.numeric(unlist(strsplit(df$V1[1],":"))[3])
	Pileup_seconds <- as.numeric(unlist(strsplit(df$V1[1],":"))[4])
	Pileup_time <- Pileup_hours*60*60+Pileup_minutes*60+Pileup_seconds
	demultiplex_hours <- as.numeric(unlist(strsplit(df$V1[2],":"))[2])
	demultiplex_minutes <- as.numeric(unlist(strsplit(df$V1[2],":"))[3])
	demultiplex_seconds <- as.numeric(unlist(strsplit(df$V1[2],":"))[4])
	demultiplex_time <- demultiplex_hours*60*60+demultiplex_minutes*60+demultiplex_seconds
	total_time <- Pileup_time + demultiplex_time
	Freemuxlet_runtime_df <- data.frame(nreads=nreads[i],method="Freemuxlet",total_seconds=total_time) %>% rbind(Freemuxlet_runtime_df,.)
}

Vireo_runtime_df <- c()
for (i in 1:length(nreads)){
	df <- read.table(str_c("~/MitoSort/simulation/Vireo/runtime/",nreads[i],"_runtime.txt"),header=FALSE,sep="\t")
	genotype_hours <- as.numeric(unlist(strsplit(df$V1[1],":"))[2])
	genotype_minutes <- as.numeric(unlist(strsplit(df$V1[1],":"))[3])
	genotype_seconds <- as.numeric(unlist(strsplit(df$V1[1],":"))[4])
	genotype_time <- genotype_hours*60*60+genotype_minutes*60+genotype_seconds
	demultiplex_hours <- as.numeric(unlist(strsplit(df$V1[2],":"))[2])
	demultiplex_minutes <- as.numeric(unlist(strsplit(df$V1[2],":"))[3])
	demultiplex_seconds <- as.numeric(unlist(strsplit(df$V1[2],":"))[4])
	demultiplex_time <- demultiplex_hours*60*60+demultiplex_minutes*60+demultiplex_seconds
	total_time <- genotype_time + demultiplex_time
	Vireo_runtime_df <- data.frame(nreads=nreads[i],method="Vireo",total_seconds=total_time) %>% rbind(Vireo_runtime_df,.)
}


runtime_df <- do.call(rbind,list(MitoSort_runtime_df,soupercell_runtime_df,Freemuxlet_runtime_df,Vireo_runtime_df))
runtime_df$total_minutes <- runtime_df$total_seconds/60
runtime_df$method <- factor(runtime_df$method,levels=c("MitoSort","Souporcell","Vireo","Freemuxlet"))
pdf("~/MitoSort/output/figure2/different_methods_runtime.pdf",height=4,width=4)
ggplot(data=runtime_df,aes(x=nreads,y=total_minutes,color=method))+
	geom_point() +
	geom_line() + 
	labs(x="Number of reads",y="Minites",title="Runtime",color="Approach")+
	scale_y_continuous(limits=c(0,85),breaks = seq(0,85,10))+
	scale_color_manual(values=brewer.pal(8,"Dark2")[c(1:3,7)])+
	theme_classic() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.4)), axis.text.x = element_text(color="black",size=rel(1.4)),axis.title.x=element_text(color="black",size=rel(1.5)),axis.title.y=element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position=c(0.25,0.8),plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm"),legend.text=element_text(color="black",size=rel(1.1)),legend.title=element_blank())
dev.off()


## Singlet/doublet classification - ROC plot and AUC calculation
assignment_df.ls <- lapply(1:length(nreads),function(i){
	soupercell_df <- soupercell_result.ls[[i]]
	colnames(soupercell_df) <- str_c("souporcell_",colnames(soupercell_df))

	Freemuxlet_df <- Freemuxlet_result.ls[[i]]
	colnames(Freemuxlet_df) <- str_c("Freemuxlet_",colnames(Freemuxlet_df))

	Vireo_df <- Vireo_result.ls[[i]]
	colnames(Vireo_df) <- str_c("Vireo_",colnames(Vireo_df))

	MitoSort_df <- MitoSort_result.ls[[i]]
	colnames(MitoSort_df) <- str_c("MitoSort_",colnames(MitoSort_df))

	true_label_df <- true_labels.ls[[i]][,c("barcode","CB_barcode","type")]
	true_label_df$CR_barcode <- sapply(1:nrow(true_label_df),function(j){
		unlist(strsplit(true_label_df$barcode[j],"_"))[2]
		})
	true_label_df$CB_barcode <- sapply(1:nrow(true_label_df),function(j){
		unlist(strsplit(true_label_df$CB_barcode[j],"_"))[2]
		})
	true_label_df$sample <- sapply(1:nrow(true_label_df),function(j){
		unlist(strsplit(true_label_df$barcode[j],"_"))[1]
		})

	assignment_df <- merge(true_label_df,soupercell_df,by.x="CB_barcode",by.y="souporcell_barcode")
	assignment_df <- merge(assignment_df,Freemuxlet_df,by.x="CB_barcode",by.y="Freemuxlet_BARCODE")
	assignment_df <- merge(assignment_df,Vireo_df,by.x="CB_barcode",by.y="Vireo_cell")
	assignment_df <- merge(assignment_df,MitoSort_df,by.x="CR_barcode",by.y="MitoSort_Barcode")
	assignment_df
	})

# FigureS1,ROC plot (doublet classification)
for (i in 1:length(nreads)){
	labels <- as.vector(factor(assignment_df.ls[[i]]$type,levels=c("Singlet","Doublet"),labels=c(1,0)))
	MitoSort_scores <- assignment_df.ls[[i]]$MitoSort_P_value_1
	Souporcell_scores <- exp(assignment_df.ls[[i]]$souporcell_log_prob_singleton)
	Freemuxlet_scores <- assignment_df.ls[[i]]$Freemuxlet_SNG.POSTERIOR
	Vireo_scores <- assignment_df.ls[[i]]$Vireo_prob_max

	# Calculate ROC curve
	MitoSort_roc_obj <- roc(labels, MitoSort_scores)
	Souporcell_roc_obj <- roc(labels, Souporcell_scores)
	Freemuxlet_roc_obj <- roc(labels, Freemuxlet_scores)
	Vireo_roc_obj <- roc(labels, Vireo_scores)

	# Calculate AUC
	MitoSort_auc_val <- auc(MitoSort_roc_obj)
	Souporcell_auc_val <- auc(Souporcell_roc_obj)
	Freemuxlet_auc_val <- auc(Freemuxlet_roc_obj)
	Vireo_auc_val <- auc(Vireo_roc_obj)

	pdf(str_c("~MitoSort/output/figure2/",nreads[i],"_reads_singlet_doublet_assignment_ROC_plot.pdf"))
	plot(MitoSort_roc_obj, col = brewer.pal(8,"Dark2")[1], main = str_c("ROC Curve - ", nreads[i]," reads per cell"),lwd=3,cex.lab=1.5, cex.axis=1.5,cex.main=1.5)
	plot(Souporcell_roc_obj, col = brewer.pal(8,"Dark2")[2], add = TRUE,lwd=3)
	plot(Vireo_roc_obj,col = brewer.pal(8,"Dark2")[3], add = TRUE,lwd=3)
	plot(Freemuxlet_roc_obj,col = brewer.pal(8,"Dark2")[7], add = TRUE,lwd=3)
	text(x = 0.24, y = 0.4, str_c("MitoSort AUC : ",round(as.vector(MitoSort_auc_val),2)),col=brewer.pal(8,"Dark2")[1])
	text(x = 0.2, y = 0.35, str_c("Souporcell AUC : ",round(as.vector(Souporcell_auc_val),2)),col=brewer.pal(8,"Dark2")[2])
	text(x = 0.24, y = 0.3, str_c("Vireo AUC : ",round(as.vector(Vireo_auc_val),2)),col=brewer.pal(8,"Dark2")[3])
	text(x = 0.2, y = 0.25, str_c("Freemuxlet AUC : ",round(as.vector(Freemuxlet_auc_val),2)),col=brewer.pal(8,"Dark2")[7])
	legend("bottomright", legend = c("MitoSort", "Souporcell","Vireo","Freemuxlet"), col =brewer.pal(8,"Dark2")[c(1:3,7)] , lwd = 2)
	dev.off()
}

auc_df <- c()
for (i in 1:length(nreads)){
	labels <- as.vector(factor(assignment_df.ls[[i]]$type,levels=c("Singlet","Doublet"),labels=c(1,0)))
	MitoSort_scores <- assignment_df.ls[[i]]$MitoSort_P_value_1
	Souporcell_scores <- exp(assignment_df.ls[[i]]$souporcell_log_prob_singleton)
	Freemuxlet_scores <- assignment_df.ls[[i]]$Freemuxlet_SNG.POSTERIOR
	Vireo_scores <- assignment_df.ls[[i]]$Vireo_prob_max

	# Calculate ROC curve
	MitoSort_roc_obj <- roc(labels, MitoSort_scores)
	Souporcell_roc_obj <- roc(labels, Souporcell_scores)
	Freemuxlet_roc_obj <- roc(labels, Freemuxlet_scores)
	Vireo_roc_obj <- roc(labels, Vireo_scores)

	# Calculate AUC
	MitoSort_auc_val <- auc(MitoSort_roc_obj)
	Souporcell_auc_val <- auc(Souporcell_roc_obj)
	Freemuxlet_auc_val <- auc(Freemuxlet_roc_obj)
	Vireo_auc_val <- auc(Vireo_roc_obj)

	auc_df <- data.frame(nreads=nreads[i],method=c("MitoSort","Souporcell","Freemuxlet","Vireo"),AUC=c(as.vector(MitoSort_auc_val),as.vector(Souporcell_auc_val),as.vector(Freemuxlet_auc_val),as.vector(Vireo_auc_val))) %>% rbind(auc_df,.)
}

# Figure2c - AUC for singlet/doublet detection
auc_df$method <- factor(auc_df$method,levels=c("MitoSort","Souporcell","Vireo","Freemuxlet"))
pdf("~/MitoSort/output/figure2/different_methods_singlet_doublet_classification_AUC.pdf",height=4,width=4)
ggplot(data=auc_df,aes(x=nreads,y=AUC,color=method))+
	geom_point() +
	geom_line() + 
	scale_y_continuous(limits=c(0,1.005),breaks = seq(0,1,0.2))+
	labs(x="Number of reads",y="Area under the curve (AUC)",title="Doublet detection",color="Approach")+
	scale_color_manual(values=brewer.pal(8,"Dark2")[c(1:3,7)])+
	theme_classic() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.4)), axis.text.x = element_text(color="black",size=rel(1.4)),axis.title.x=element_text(color="black",size=rel(1.5)),axis.title.y=element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position=c(0.8,0.21),plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm"),legend.text=element_text(color="black",size=rel(1.1)),legend.title=element_blank())
dev.off()

## Figure2d, singlet assignment - the adjusted Rand index 
library(fpc)
ari_df <- c()

for (i in 1:length(nreads)){
	Souporcell_singlet_df <- assignment_df.ls[[i]] %>%
		filter(souporcell_status=="singlet",type=="Singlet")
	ari_df <- data.frame(nreads=nreads[i],method="Souporcell",ARI=mclust::adjustedRandIndex(Souporcell_singlet_df$sample, Souporcell_singlet_df$souporcell_assignment)) %>% rbind(ari_df,.)

	Freemuxlet_singlet_df <- assignment_df.ls[[i]] %>%
		filter(Freemuxlet_DROPLET.TYPE=="SNG",type=="Singlet")
	ari_df <- data.frame(nreads=nreads[i],method="Freemuxlet",ARI=mclust::adjustedRandIndex(Freemuxlet_singlet_df$sample, Freemuxlet_singlet_df$Freemuxlet_SNG.BEST.GUESS)) %>% rbind(ari_df,.)

	Vireo_singlet_df <- assignment_df.ls[[i]] %>%
		filter(Vireo_donor_id %in% str_c("donor",0:5),type=="Singlet")
	ari_df <- data.frame(nreads=nreads[i],method="Vireo",ARI=mclust::adjustedRandIndex(Vireo_singlet_df$sample, Vireo_singlet_df$Vireo_donor_id)) %>% rbind(ari_df,.)

	
	MitoSort_singlet_df <- assignment_df.ls[[i]] %>%
		filter(MitoSort_Demultiplex!="Doublet",type=="Singlet")
	ari_df <- data.frame(nreads=nreads[i],method="MitoSort",ARI=mclust::adjustedRandIndex(MitoSort_singlet_df$sample,MitoSort_singlet_df$MitoSort_Demultiplex)) %>% rbind(ari_df,.)
}


ari_df$method <- factor(ari_df$method,levels=c("MitoSort","Souporcell","Vireo","Freemuxlet"))
pdf("~/MitoSort/output/figure2/different_methods_singlet_assignment_ARI.pdf",height=4,width=4)
ggplot(data=ari_df,aes(x=nreads,y=ARI,color=method))+
	geom_point() +
	geom_line() + 
	scale_y_continuous(limits=c(0,1.005),breaks = seq(0,1,0.2))+
	labs(x="Number of reads",y="Adjusted Rand index",title="Singlet assignment",color="Approach")+
	scale_color_manual(values=brewer.pal(8,"Dark2")[c(1:3,7)])+
	theme_classic() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.4)), axis.text.x = element_text(color="black",size=rel(1.4)),axis.title.x=element_text(color="black",size=rel(1.5)),axis.title.y=element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position=c(0.8,0.4),plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm"),legend.text=element_text(color="black",size=rel(1.1)),legend.title=element_blank())
dev.off()

## Heatmap showing VAF for sample-specific variants
# Figure2e
# MitoSort
MitoSort_ref_mtx <- read.csv("~/MitoSort/simulation/MitoSort/50000/MitoSort/SNP_matrix/ref.csv",row.names=1)
MitoSort_alt_mtx <- read.csv("~/MitoSort/simulation/MitoSort/50000/MitoSort/SNP_matrix/alt.csv",row.names=1)
MitoSort_freq_mtx <- read.csv("~/MitoSort/simulation/MitoSort/50000/MitoSort/SNP_matrix/frequency.csv",row.names=1)
MitoSort_freq_mtx[is.na(MitoSort_freq_mtx)] <- 0
MitoSort_specific_germline <- read.table("~/MitoSort/simulation/MitoSort/50000/MitoSort/Demultiplex_output/specific_germline.txt",header=TRUE,sep="\t")
MitoSort_specific_germline <- MitoSort_specific_germline %>%
	arrange(Sample)
MitoSort_specific_germline <- MitoSort_specific_germline[-132,]

MitoSort_cluster_df <- read.table("~/MitoSort/simulation/MitoSort/50000/MitoSort/Demultiplex_output/result_pvalue.txt",header=TRUE,sep="\t")
MitoSort_singlet_df <- MitoSort_cluster_df %>%
	filter(Demultiplex %in% str_c("Sample",0:5)) %>%
	arrange(Demultiplex)

MitoSort_specfic_germline_freq_mtx <- MitoSort_freq_mtx[MitoSort_specific_germline$Specific.germline,MitoSort_singlet_df$Barcode]
anno_df <- data.frame(row.names=MitoSort_singlet_df$Barcode,MitoSort_assignment=MitoSort_singlet_df$Demultiplex)
ha <- HeatmapAnnotation(
  df = anno_df,
  col = list(
    MitoSort_assignment=c("Sample0"=brewer.pal(8,"Dark2")[1],"Sample1"=brewer.pal(8,"Dark2")[2],"Sample2"=brewer.pal(8,"Dark2")[3],"Sample3"=brewer.pal(8,"Dark2")[4],"Sample4"=brewer.pal(8,"Dark2")[5],"Sample5"=brewer.pal(8,"Dark2")[6])
    )
  )
col <- brewer.pal(9,"Reds")
pdf("~/MitoSort/simulation/MitoSort/50000/MitoSort/Demultiplex_output/singlet_specific_germline_variant_heatmap.pdf",width=10.5,height=9)
h1 <- Heatmap(MitoSort_specfic_germline_freq_mtx,
  name="Frequency",
  col=col,
  cluster_row_slices=FALSE,
  cluster_column_slices=FALSE,
  row_split = c(rep("",22),rep("",7),rep("",55),rep("",21),rep("",16),rep("",16)),
  column_split=c(rep("Sample0",98),rep("Sample1",97),rep("Sample2",99),rep("Sample3",99),rep("Sample4",96),rep("Sample5",95)),
  show_row_names=FALSE,
  column_names_gp = gpar(fontsize = 9),
  cluster_rows=FALSE,
  cluster_columns=TRUE,
  show_column_names=FALSE,
  top_annotation=ha,
  row_names_gp = gpar(fontsize = 7.5)) 
h1 <- draw(h1)
dev.off()

# Figure2g, tsne for MitoSort sample-specific variants
library(Rtsne)
tSNE_fit <- as.data.frame(t(MitoSort_specfic_germline_freq_mtx)) %>%
  dplyr::select(where(is.numeric)) %>%
  Rtsne::Rtsne(check_duplicates = FALSE,initial_dims=50,perplexity=60)
tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1="V1",tSNE2="V2") 
rownames(tSNE_df) <- colnames(MitoSort_specfic_germline_freq_mtx) 
tSNE_df <- merge(tSNE_df,assignment_df.ls[[9]][,c("CR_barcode","sample","souporcell_assignment","MitoSort_Demultiplex")],by.x="row.names",by.y="CR_barcode")
tSNE_df$sample <- factor(tSNE_df$sample,levels=c("sample1","CCL1","CD34","CRC","BMMC","lib2"),labels=c(str_c("Sample",0:5)))
newpalette <- brewer.pal(8,"Dark2")[1:6]
pdf("~/MitoSort/simulation/MitoSort/50000/MitoSort/Demultiplex_output/singlet_specific_germline_variant_tSNE_with_sample_label.pdf",width=8)
ggplot(data=tSNE_df,aes(x=tSNE1,y=tSNE2))+
	geom_point(aes(color=sample),size=2)+
	scale_color_manual(values=newpalette)+
	labs(x="tSNE_01",y="tSNE_02",title="mtDNA variants",color="Sample")+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',size=1.3),plot.title = element_text(hjust = 0.5,size=rel(2)),axis.text=element_text(color="black",size=rel(1.2)),axis.ticks=element_blank(),axis.title=element_text(color="black",size=rel(1.3)),plot.margin=unit(c(0.2, 0.2, 0.7, 0.5), "cm"),legend.title = element_text(size=15),legend.text=element_text(size=10))
dev.off()

# Figure2i
# Compute pairwise Euclidean distances between points in the t-SNE space
# Function to calculate Euclidean distance between two points
euclidean_distance <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
}

# Calculate distances among points in t-SNE space
library(proxy)
distances <- proxy::dist(tSNE_fit$Y, method = euclidean_distance)
distance_matrix <- as.matrix(distances)
rownames(distance_matrix) <- colnames(MitoSort_specfic_germline_freq_mtx) 
colnames(distance_matrix) <- colnames(MitoSort_specfic_germline_freq_mtx) 

n <- nrow(distance_matrix)
row_indices <- c()
col_indices <- c()
dist_values <- vector("numeric")

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    row_indices <- c(row_indices, rownames(distance_matrix)[i])
    col_indices <- c(col_indices, colnames(distance_matrix)[j])
    dist_values <- c(dist_values, distance_matrix[i, j])
  }
}

distance_df <- data.frame(
  cell_1 = row_indices,
  cell_2 = col_indices,
  distance = dist_values
)
distance_df <- merge(distance_df,assignment_df.ls[[9]][,c("CR_barcode","sample","souporcell_assignment","MitoSort_Demultiplex")],by.x="cell_1",by.y="CR_barcode",all.x=TRUE,sort=FALSE)
distance_df <- distance_df %>%
	rename(cell_1_group=sample)
distance_df <- merge(distance_df,assignment_df.ls[[9]][,c("CR_barcode","sample","souporcell_assignment","MitoSort_Demultiplex")],by.x="cell_2",by.y="CR_barcode",all.x=TRUE,sort=FALSE)
distance_df <- distance_df %>%
	rename(cell_2_group=sample)
distance_df$type <- ifelse(distance_df$cell_1_group==distance_df$cell_2_group,"within sample","between sample")
distance_df$type <- factor(distance_df$type,levels=c("within sample","between sample"))
#newpalette <- brewer.pal(9,"Pastel1")[1:2]
newpalette <- c("firebrick3","dodgerblue3")
pdf("~/MitoSort/simulation/MitoSort/50000/MitoSort/Demultiplex_output/singlet_specific_germline_variant_tSNE_within_vs_between_sample_distance_violin_plot.pdf")
ggplot(data=distance_df,aes(x=type,y=distance,fill=type))+
	#geom_violin(scale="width",trim=FALSE)+
    geom_boxplot(width=0.3,outlier.size=0.8)+
    labs(x="",title="mtDNA variants",y="Euclidean distance")+
    scale_fill_manual(values=newpalette)+
	theme_classic()+
	theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.4)), axis.text.x = element_text(color="black",size=rel(1.4)),axis.title.x=element_text(color="black",size=rel(1.5)),axis.title.y=element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position="non",plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm"),legend.text=element_text(color="black",size=rel(1.1)),legend.title=element_blank())
dev.off()

## Figure2f
# souporcell 
ref_mtx <- Matrix::readMM('~/MitoSort/simulation/souporcell/50000/ref.mtx')
alt_mtx <- Matrix::readMM('~/MitoSort/simulation/souporcell/50000/alt.mtx')
freq_mtx <- alt_mtx / (ref_mtx+alt_mtx)
freq_mtx[is.na(freq_mtx)] <- 0

variants_df <- read.table("~/MitoSort/simulation/souporcell/50000/souporcell_merged_sorted_vcf.vcf",sep="\t")
variants <- str_c(variants_df$V1,"-",variants_df$V2,"-",variants_df$V4,"-",variants_df$V5)
rownames(freq_mtx) <- variants

clusters_df <- read.table("~/MitoSort/simulation/souporcell/50000/clusters.tsv",header=TRUE)
colnames(freq_mtx) <- clusters_df$barcode
singlet_clusters_df <- clusters_df %>%
	dplyr::filter(status=="singlet")

# sleect sample-specific variants
cluster_genotypes_df <- read.table("~/MitoSort/simulation/souporcell/50000/cluster_genotypes.vcf",sep="\t")
cluster_genotypes_df$cluster0_GT <- sapply(1:nrow(cluster_genotypes_df),function(i){
	ifelse(unlist(strsplit(cluster_genotypes_df$V10[i],":"))[1]=="1/1",1,0)
	})
cluster_genotypes_df$cluster1_GT <- sapply(1:nrow(cluster_genotypes_df),function(i){
	ifelse(unlist(strsplit(cluster_genotypes_df$V11[i],":"))[1]=="1/1",1,0)
	})
cluster_genotypes_df$cluster2_GT <- sapply(1:nrow(cluster_genotypes_df),function(i){
	ifelse(unlist(strsplit(cluster_genotypes_df$V12[i],":"))[1]=="1/1",1,0)
	})
cluster_genotypes_df$cluster3_GT <- sapply(1:nrow(cluster_genotypes_df),function(i){
	ifelse(unlist(strsplit(cluster_genotypes_df$V13[i],":"))[1]=="1/1",1,0)
	})
cluster_genotypes_df$cluster4_GT <- sapply(1:nrow(cluster_genotypes_df),function(i){
	ifelse(unlist(strsplit(cluster_genotypes_df$V14[i],":"))[1]=="1/1",1,0)
	})
cluster_genotypes_df$cluster5_GT <- sapply(1:nrow(cluster_genotypes_df),function(i){
	ifelse(unlist(strsplit(cluster_genotypes_df$V15[i],":"))[1]=="1/1",1,0)
	})
cluster_genotypes_df$all_GT_sum <- cluster_genotypes_df$cluster0_GT+cluster_genotypes_df$cluster1_GT+cluster_genotypes_df$cluster2_GT+cluster_genotypes_df$cluster3_GT+cluster_genotypes_df$cluster4_GT+cluster_genotypes_df$cluster5_GT

specfic_germline_df <- c()
for (i in 0:5){
	sample_specfic_df <- cluster_genotypes_df[which(cluster_genotypes_df$all_GT_sum==1 & cluster_genotypes_df[,str_c("cluster",i,"_GT")]==1),]
	sample_specfic_variants <- str_c(sample_specfic_df$V1,"-",sample_specfic_df$V2,"-",sample_specfic_df$V4,"-",sample_specfic_df$V5)
	specfic_germline_df <- data.frame(variant=sample_specfic_variants,sample=str_c("cluster",i)) %>% rbind(specfic_germline_df,.)
}


souporcell_specfic_germline_freq_mtx <- as.matrix(freq_mtx[which(rownames(freq_mtx) %in% specfic_germline_df$variant),which(colnames(freq_mtx) %in% singlet_clusters_df$barcode)])

souporcell_specfic_germline_df <- c()
for (i in 0:5){
	sample_specfic_df <- cluster_genotypes_df[which(cluster_genotypes_df$all_GT_sum==1 & cluster_genotypes_df[,str_c("cluster",i,"_GT")]==1),]
	sample_specfic_variants <- str_c(sample_specfic_df$V1,"-",sample_specfic_df$V2,"-",sample_specfic_df$V4,"-",sample_specfic_df$V5)
	souporcell_specfic_germline_df <- data.frame(variant=sample_specfic_variants,sample=str_c("cluster",i)) %>% rbind(souporcell_specfic_germline_df,.)
}

souporcell_clusters_df <- read.table("~/MitoSort/simulation/souporcell/50000/clusters.tsv",header=TRUE)
souporcell_singlet_clusters_df <- souporcell_clusters_df %>%
	filter(status=="singlet") %>%
	arrange(assignment)
souporcell_singlet_clusters_df <- merge(souporcell_singlet_clusters_df,)	
table(assignment_df.ls[[9]]$souporcell_assignment,assignment_df.ls[[9]]$MitoSort_Demultiplex)
table(assignment_df.ls[[9]]$souporcell_assignment,assignment_df.ls[[9]]$sample)
table(assignment_df.ls[[9]]$MitoSort_Demultiplex,assignment_df.ls[[9]]$sample)


souporcell_specfic_germline_freq_mtx <- as.matrix(freq_mtx[souporcell_specfic_germline_df$variant,souporcell_singlet_clusters_df$barcode])
library(ComplexHeatmap)
anno_df <- data.frame(row.names=souporcell_singlet_clusters_df$barcode,souporcell_assignment=as.vector(factor(souporcell_singlet_clusters_df$assignment,levels=0:5,labels=str_c("Cluster",0:5))))
ha <- HeatmapAnnotation(
  df = anno_df,
  col = list(
    souporcell_assignment=c("Cluster0"=brewer.pal(8,"Dark2")[1],"Cluster1"=brewer.pal(8,"Dark2")[2],"Cluster2"=brewer.pal(8,"Dark2")[3],"Cluster3"=brewer.pal(8,"Dark2")[4],"Cluster4"=brewer.pal(8,"Dark2")[5],"Cluster5"=brewer.pal(8,"Dark2")[6])
    )
  )
col <- brewer.pal(9,"Reds")
pdf("~/MitoSort/simulation/souporcell/50000/singlet_specific_germline_variant_heatmap.pdf",width=10.5,height=9)
h1 <- Heatmap(souporcell_specfic_germline_freq_mtx,
  name="Frequency",
  col=col,
  cluster_row_slices=FALSE,
  cluster_column_slices=FALSE,
  row_split = c(rep("",2470),rep("",856),rep("",1090),rep("",778),rep("",1125),rep("",862)),
  column_split=c(rep("Cluster0",100),rep("Cluster1",47),rep("Cluster2",145),rep("Cluster3",96),rep("Cluster4",110),rep("Cluster5",54)),
  show_row_names=FALSE,
  column_names_gp = gpar(fontsize = 9),
  cluster_rows=FALSE,
  cluster_columns=TRUE,
  show_column_names=FALSE,
  top_annotation=ha,
  row_names_gp = gpar(fontsize = 7.5)) 
h1 <- draw(h1)
dev.off()


## Figure2h
# tsne for souporcell  sample-specific variants
souporcell_specfic_germline_freq_mtx <- as.matrix(freq_mtx[souporcell_specfic_germline_df$variant,souporcell_singlet_clusters_df$barcode])
library(Rtsne)
tSNE_fit <- as.data.frame(t(souporcell_specfic_germline_freq_mtx)) %>%
  dplyr::select(where(is.numeric)) %>%
  Rtsne::Rtsne(check_duplicates = FALSE,initial_dims=50,perplexity=60)
tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1="V1",tSNE2="V2") 
rownames(tSNE_df) <- colnames(souporcell_specfic_germline_freq_mtx) 
tSNE_df <- merge(tSNE_df,assignment_df.ls[[9]][,c("CB_barcode","sample","souporcell_assignment","MitoSort_Demultiplex")],by.x="row.names",by.y="CB_barcode")
tSNE_df$sample <- factor(tSNE_df$sample,levels=c("sample1","CCL1","CD34","CRC","BMMC","lib2"),labels=c(str_c("Sample",0:5)))
newpalette <- brewer.pal(8,"Dark2")[1:6]
pdf("/md01/shipy3/tzj/SMILE/simulation/souporcell/50000/singlet_specific_germline_variant_tSNE_with_sample_label.pdf",width=8)
ggplot(data=tSNE_df,aes(x=tSNE1,y=tSNE2))+
	geom_point(aes(color=sample),size=2)+
	scale_color_manual(values=newpalette)+
	labs(x="tSNE_01",y="tSNE_02",title="Nuclear variants",color="Sample")+
	#geom_segment(aes(x = -20, y = -30, xend = -10, yend = -30),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  	#geom_segment(aes(x = -20, y = -30, xend = -20, yend = -20),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  	#annotate(geom = "text", x = -17, y = -32, label = "tSNE_1", color = "black") +
  	#annotate(geom = "text", x = -20, y = -25, label = "tSNE_2", color = "black",angle = 90)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',size=1.3),plot.title = element_text(hjust = 0.5,size=rel(2)),axis.text=element_text(color="black",size=rel(1.2)),axis.ticks=element_blank(),axis.title=element_text(color="black",size=rel(1.3)),plot.margin=unit(c(0.2, 0.2, 0.7, 0.5), "cm"),legend.title = element_text(size=15),legend.text=element_text(size=10))
dev.off()


## Figure2j
# Compute pairwise Euclidean distances between points in the t-SNE space
# Function to calculate Euclidean distance between two points
euclidean_distance <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
}

# Calculate distances among points in t-SNE space
library(proxy)
distances <- proxy::dist(tSNE_fit$Y, method = euclidean_distance)
distance_matrix <- as.matrix(distances)
rownames(distance_matrix) <- colnames(souporcell_specfic_germline_freq_mtx) 
colnames(distance_matrix) <- colnames(souporcell_specfic_germline_freq_mtx) 

n <- nrow(distance_matrix)
row_indices <- c()
col_indices <- c()
dist_values <- vector("numeric")

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    row_indices <- c(row_indices, rownames(distance_matrix)[i])
    col_indices <- c(col_indices, colnames(distance_matrix)[j])
    dist_values <- c(dist_values, distance_matrix[i, j])
  }
}

distance_df <- data.frame(
  cell_1 = row_indices,
  cell_2 = col_indices,
  distance = dist_values
)
distance_df <- merge(distance_df,assignment_df.ls[[9]][,c("CB_barcode","sample","souporcell_assignment","MitoSort_Demultiplex")],by.x="cell_1",by.y="CB_barcode",all.x=TRUE,sort=FALSE)
distance_df <- distance_df %>%
	rename(cell_1_group=sample)
distance_df <- merge(distance_df,assignment_df.ls[[9]][,c("CB_barcode","sample","souporcell_assignment","MitoSort_Demultiplex")],by.x="cell_2",by.y="CB_barcode",all.x=TRUE,sort=FALSE)
distance_df <- distance_df %>%
	rename(cell_2_group=sample)
distance_df$type <- ifelse(distance_df$cell_1_group==distance_df$cell_2_group,"within sample","between sample")
distance_df$type <- factor(distance_df$type,levels=c("within sample","between sample"))
newpalette <- c("firebrick3","dodgerblue3")
pdf("/md01/shipy3/tzj/SMILE/simulation/souporcell/50000/singlet_specific_germline_variant_tSNE_within_vs_between_sample_distance_violin_plot.pdf")
ggplot(data=distance_df,aes(x=type,y=distance,fill=type))+
	#geom_violin(scale="width",trim=FALSE)+
    geom_boxplot(width=0.3,outlier.size=0.8)+
    labs(x="",title="Nuclear variants",y="Euclidean distance")+
    scale_fill_manual(values=newpalette)+
	theme_classic()+
	theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.4)), axis.text.x = element_text(color="black",size=rel(1.4)),axis.title.x=element_text(color="black",size=rel(1.5)),axis.title.y=element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position="non",plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm"),legend.text=element_text(color="black",size=rel(1.1)),legend.title=element_blank())
dev.off()
```