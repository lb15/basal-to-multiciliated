### proportion testing
 library(Seurat)
 library("scProportionTest")
 library(colorspace)

 seur=readRDS("v1/early_int_v1.rds")
 

seur$merge_mcc <-seur$ord_clus
seur$merge_mcc[seur$ord_clus ==4] <- 3
seur$merge_mcc[seur$ord_clus ==5] <- 3

DimPlot(seur, group.by="merge_mcc",label=T)
prop_test <- sc_utils(seur)
prop_test <- permutation_test(
         prop_test, cluster_identity = "merge_mcc",
         sample_1 = "NT1_KO1_soupx", sample_2 = "Gmnc_soupx",
         sample_identity = "orig.ident"
 ) 

 permutation_plot(prop_test) 
 
 prop_test2 <- permutation_test(
         prop_test, cluster_identity = "merge_mcc",
         sample_1 = "NT1_KO1_soupx", sample_2 = "Mcidas_soupx",
         sample_identity = "orig.ident"
 ) 
 
 log2FD_threshold = log2(1.5)
 FDR_threshold=0.05
 p=permutation_plot(prop_test2,log2FD_threshold = log2(1.5),order_clusters = T) 
 g=permutation_plot(prop_test,log2FD_threshold = log2(1.5),order_clusters = T) 
 pdf("early_int_v1_proportion_permtest_Mcidas.pdf",height=4,width=4,useDingbats = F)
p + labs(color="Significance") + scale_color_manual(values=c("salmon","grey"),breaks=c(paste("FDR <",FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), "n.s."), labels = c("FDR < 0.05 & FC > 1.5", "n.s"))+ylab("log2(FC)")+xlab("Clusters")+theme_classic()+ggtitle("Mcidas KO vs. NT KO")
dev.off()

pdf("early_int_v1_proportion_permtest_Gmnc.pdf",height=4,width=4,useDingbats = F)
g + labs(color="Significance") + scale_color_manual(values=c("salmon","grey"),breaks=c(paste("FDR <",FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), "n.s."), labels = c("FDR < 0.05 & FC > 1.5", "n.s"))+ylab("log2(FC)")+xlab("Clusters")+theme_classic()+ggtitle("Gmnc KO vs. NT KO")
dev.off()

## make volcano plot

data=p$data
data$comparison <- rep("Mcidas_vs_NT",length(rownames(data)))
colnames(data)[2:3] <- c("sample1","sample2")
data_2=g$data
data_2$comparison <- rep("Gmnc_vs_NT",length(rownames(data_2)))
colnames(data_2)[2:3] <- c("sample1","sample2")
full_dat = rbind(data, data_2)

full_dat$sig_clus <- factor(paste(full_dat$clusters,full_dat$significance, sep="_"),levels=c("1_n.s.","2_n.s.","2_FDR < 0.05 & abs(Log2FD) > 0.58","3_FDR < 0.05 & abs(Log2FD) > 0.58","4_FDR < 0.05 & abs(Log2FD) > 0.58","5_FDR < 0.05 & abs(Log2FD) > 0.58","6_n.s.","7_n.s.","7_FDR < 0.05 & abs(Log2FD) > 0.58","8_n.s.","9_n.s.","10_n.s.","11_n.s."))

mcc_color <- colorRampPalette(c("lightskyblue2", "royalblue4"))

fig1_col = c("#F7A12D","#009245",mcc_color(3),"palevioletred2","palegreen1","darkgoldenrod","#998675","#996699","tomato3")

desat_col <- lighten(fig1_col, amount=0)

## dont like this plot - too busy
ggplot(full_dat, aes(x=obs_log2FD,y=-log(FDR,10),shape=comparison,fill=sig_clus,colour=significance))+
        geom_point(size=4, alpha=1,stroke=1.5)+
        scale_fill_manual(values=c(desat_col[1:2],fig1_col[2:5],desat_col[6:7],fig1_col[7],desat_col[8:11]))+
        scale_shape_manual(values=c(21,24),name = "", labels = c("Gmnc KO vs. NT KO", "Mcidas KO vs. NT KO"))+
        scale_color_manual(values=c("black","grey"))+
        theme_classic()+
        ylab("-log(FDR)")+
        xlab("log2(FC)")

## combine permutation plots
full_dat$clusters <- factor(full_dat$clusters, levels=c(1:11))

pdf("early_int_v1_permuationtest_full.pdf",height=6,width=8,useDingbats = F)
ggplot(full_dat, aes(x = clusters, y = obs_log2FD,shape=comparison)) + 
        geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, 
                            color = significance),size=1,position = position_dodge(width = 0.7)) + 
        theme_bw() + 
        geom_hline(yintercept = log2FD_threshold, lty = 2) + 
        geom_hline(yintercept = -log2FD_threshold,  lty = 2) + 
        geom_hline(yintercept = 0,size=0.5) + 
        scale_color_manual(values = c("salmon",  "grey")) +
        coord_flip()+
        labs(color="Significance",shape="Comparison") + 
        scale_color_manual(values=c("salmon","grey"),breaks=c(paste("FDR <",FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), "n.s."), labels = c("FDR < 0.05 & FC > 1.5", "n.s"))+ylab("log2(FC)")+
        xlab("Clusters")+
        theme(axis.text = element_text(size=12),axis.title = element_text(size=14),panel.grid = element_blank())+
        geom_vline(xintercept = seq(0.5, length(full_dat$clusters), by = 1), color="gray", size=.5, alpha=.5)
dev.off()

sub_dat = filter(full_dat, clusters %in% c(1:6))

pdf("Figure/early_int_v1_permuationtest_clusters_mergeMCC.pdf",height=4.5,width=4,useDingbats = F)
ggplot(sub_dat, aes(x = clusters, y = obs_log2FD,shape=comparison,fill=comparison)) + 
        geom_bar(stat="identity",position="dodge")+
        geom_errorbar(aes(ymin=boot_CI_2.5,
                          ymax=boot_CI_97.5), position = "dodge")+
        theme_bw() + 
        geom_hline(yintercept = log2FD_threshold, lty = 2) + 
        geom_hline(yintercept = -log2FD_threshold,  lty = 2) + 
        geom_hline(yintercept = 0,size=0.5) + 
        scale_fill_manual(values = c("#205c2e",  "#f06324")) +
        labs(color="Significance",shape="Comparison") + 
        scale_color_manual(breaks=c(paste("FDR <",FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), "n.s."), labels = c("FDR < 0.05 & FC > 1.5", "n.s"))+ylab("log2(FC)")+
        xlab("Clusters")+
        theme(axis.text = element_text(size=12),axis.title = element_text(size=14),panel.grid = element_blank())+
        geom_vline(xintercept = seq(0.5, length(full_dat$clusters), by = 1), color="gray", size=.5, alpha=.5)
dev.off()

