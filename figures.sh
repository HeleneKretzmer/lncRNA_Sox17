### Lineage tree heatmap showing SOX17 (green) and LNCNSOX17 (orange) expression
```R
require(reshape2)
require(gplots)

n <- c('Homo_sapiens_heart_tissue_male_adult_(34_years)', 'Homo_sapiens_heart_left_ventricle_tissue_female_adult_(51_year)', 'Homo_sapiens_bladder_microvascular_endothelial_primary_cell_male_adult_(46_yea
rs)', 'Homo_sapiens_cardiac_atrium_fibroblast_primary_cell_male_child_(2_years)', 'Homo_sapiens_cardiac_ventricle_fibroblast_primary_cell_male_adult_(18_years)', 'Homo_sapiens_dermis_blood_vessel_endothel
ial_primary_cell_female_child_(16_years)', 'Homo_sapiens_endothelial_of_coronary_artery_primary_cell_female_adult_(41_year)', 'Homo_sapiens_fibroblast_of_lung_primary_cell_male_adult_(23_years)', 'Homo_sa
piens_lung_microvascular_endothelial_primary_cell_female_adult_(55_years)', 'Homo_sapiens_pulmonary_artery_endothelial_primary_cell_male_adult_(23_years)', 'Homo_sapiens_tracheal_epithelial_primary_cell_m
ale_adult_(21_year)', 'Homo_sapiens_thoracic_aorta_endothelial_primary_cell_female_adult_(22_years)', 'Homo_sapiens_vein_endothelial_primary_cell_male_adult_(48_years)')

# load TPM from Encode (stringtie output)
data1 <- read.table('Encode_TPM.txt', header=F, col.names=c('sample','gene','TPM'))

# load TPM from publication (stringtie output)
data2 <- read.table('BingLim_EndodemInduction_CellStemCell_2014_TPM.txt', header=F, col.names=c('sample','gene','TPM'))

# combine data sets
df <- dcast(rbind(data1,data2), gene~sample)

# plot
heatmap.2(as.matrix(log2(df[,-which(colnames(df) %in% n)][c(1,1),-1]+1)), labRow=df$gene[1], cexRow=1, trace='none', margins=c(25,5), col=colorpanel(100, 'snow', 'chartreuse3'))
heatmap.2(as.matrix(log2(df[,-which(colnames(df) %in% n)][c(2,2),-1]+1)), labRow=df$gene[2], cexRow=1, trace='none', margins=c(25,5), col=colorpanel(100, 'snow', 'orange1'))
```

### Heatmap showing SOX17 binding distribution genome-wide in sgCtrl and sgLNCSOX17 EN
deeptools/bin/computeMatrix scale-regions \
-S signal/sgCtrl_day5_IP_SOX17_FE.bw signal/sgLNCSOX17_day5_IP_FE.bw \
-R peaks/sgLNCSOX17_day5_IP_narrowpeak_q3_FE3.bed \
-o sgLNCSOX17_peaks.tab.gz \
-a 2500 -b 2500 \
--binSize 50

deeptools/bin/plotHeatmap \
-m sgLNCSOX17_peaks.tab.gz \
-out sgLNCSOX17_peaks.pdf \
--startLabel peak \
--endLabel peak \
--missingDataColor 'white' \
--colorMap 'Greys' \
--heatmapHeight 14


### Scatter plot highlighting differentially expressed genes between sgLNCSOX17 and sgCtrl EN cells
```R
require(ggplot2)
require(ggrepel)

data <- read.table('differential_expression/D9_log2FC.txt', header=T)
TPM <- read.table('avg_TPMs.bed', header=T, comment.char='')
colnames(TPM)[5] <- 'ENSG'

df <- merge(TPM[,c(5,9,10)], data, by='ENSG')
df$avg_expr <- rowMeans(df[,c(2,3)])
df$change <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1, ifelse(df$log2FoldChange>0,'Upregulated','Downregulated'), 'not-significant')
df$change[is.na(df$change)] <- 'not-significant'
df$change <- factor(df$change, levels=c('not-significant','Downregulated','Upregulated'))

ggplot(subset(df,avg_expr>0), aes(x=log10(avg_expr), y=log2FoldChange)) + theme_classic() + geom_hline(yintercept=0, color='grey') + scale_color_manual(values=c('lightgrey','cornflowerblue','magenta')) + geom_point(size=0.5, aes(color=change))
```


### Scatter plots displaying DNA methylation levels, ATAC signal, endoderm TFs occupancy as measured by ChIP-seq and TFs binding motifs abundance at the LNCSOX17 locus in EN cells
cat bins.bed
chr8    55139461        55140801        r1
chr8    55138192        55139461        r2
chr8    55136923        55138192        eSOX17
chr8    55135654        55136923        r3
chr8    55134385        55135654        r4
chr8    55133116        55134385        r5
chr8    55131847        55133116        r6
chr8    55130578        55131847        r7
chr8    55129309        55130578        r8
chr8    55128040        55129309        r9
chr8    55126771        55128040        r10
chr8    55125502        55126771        r11
chr8    55124233        55125502        r12
chr8    55122964        55124233        r13
chr8    55121695        55122964        r14
chr8    55120426        55121695        r15
chr8    55119157        55120426        r16
chr8    55117888        55119157        r17

#### 2) Calculate for each region the reads density/mean methyl
#### GATA4, GATA6, SOX17, FOXA2, ATAC, DNAme
##### ATAC
UCSCtools/bigWigAverageOverBed hDE_ATAC.pooled.fc.bigWig bins.bed tmp
cut -f1,5 tmp >hDE_ATAC_bins_avg.tsv
## FOXA2
UCSCtools/bigWigAverageOverBed hDE_FOXA2.pooled.fc.bigWig bins.bed tmp
cut -f1,5 tmp >hDE_FOXA2_bins_avg.tsv
## GATA4
UCSCtools/bigWigAverageOverBed hDE_GATA4.pooled.fc.bigWig bins.bed tmp
cut -f1,5 tmp >hDE_GATA4_bins_avg.tsv
## GATA6
UCSCtools/bigWigAverageOverBed hDE_GATA6.pooled.fc.bigWig bins.bed tmp
cut -f1,5 tmp >hDE_GATA6_bins_avg.tsv
## SOX17
UCSCtools/bigWigAverageOverBed sgCtrl_day5_IP_SOX17_FE.bw bins.bed tmp
cut -f1,5 tmp >hDE_SOX17_bins_avg.tsv

##### DNAme
bigWigAverageOverBed HUES_EN.bw bins.bed tmp
cut -f1,6 tmp >hDE_DNAme_bins_avg.tsv

#### 3) get TF motif file
cat homer.KnownMotifs.hg19.170917.bed/output.hg19.bed | grep Gata6 >TF_motifs.bed
cat homer.KnownMotifs.hg19.170917.bed/output.hg19.bed | grep Gata4 >>TF_motifs.bed
cat homer.KnownMotifs.hg19.170917.bed/output.hg19.bed | grep Sox17 >>TF_motifs.bed
cat homer.KnownMotifs.hg19.170917.bed/output.hg19.bed | grep Foxa2 >>TF_motifs.bed

cat homer.KnownMotifs.hg19.170917.bed/output.hg19.bed | grep Gata >TFfamily_motifs.bed
cat homer.KnownMotifs.hg19.170917.bed/output.hg19.bed | grep Sox >>TFfamily_motifs.bed
cat homer.KnownMotifs.hg19.170917.bed/output.hg19.bed | grep Foxa >>TFfamily_motifs.bed

#### count TF motifs per bin
bedtools intersect -loj -a bins.bed -b TF_motifs.bed | sort -k4,4 -k8,8 | grep -v Oct | sed 's/\t\t/\t/' | cut -f1-4,8 | bedtools groupby -g 1,2,3,4,5 -c 5 -o count >bins_TF_motifs_counts.txt
bedtools intersect -loj -a bins.bed -b TFfamily_motifs.bed | sort -k4,4 -k8,8 | grep -v Oct | sed 's/Gata.*/Gata/' | sed 's/Foxa.*/Foxa/' | sed 's/Sox.*/Sox/' | sed 's/\t\t/\t/' | cut -f1-4,8 | bedtools g
roupby -g 1,2,3,4,5 -c 5 -o count >bins_TFfamily_motifs_counts.txt


### Scatter plot showing the expression of a set of endoderm lncRNAs and the corresponding TFs
```R
require(ggplot2)
require(reshape2)
library(ggExtra)

n <- c('Homo_sapiens_heart_tissue_male_adult_(34_years)', 'Homo_sapiens_heart_left_ventricle_tissue_female_adult_(51_year)', 'Homo_sapiens_bladder_microvascular_endothelial_primary_cell_male_adult_(46_yea
rs)', 'Homo_sapiens_cardiac_atrium_fibroblast_primary_cell_male_child_(2_years)', 'Homo_sapiens_cardiac_ventricle_fibroblast_primary_cell_male_adult_(18_years)', 'Homo_sapiens_dermis_blood_vessel_endothel
ial_primary_cell_female_child_(16_years)', 'Homo_sapiens_endothelial_of_coronary_artery_primary_cell_female_adult_(41_year)', 'Homo_sapiens_fibroblast_of_lung_primary_cell_male_adult_(23_years)', 'Homo_sa
piens_lung_microvascular_endothelial_primary_cell_female_adult_(55_years)', 'Homo_sapiens_pulmonary_artery_endothelial_primary_cell_male_adult_(23_years)', 'Homo_sapiens_tracheal_epithelial_primary_cell_m
ale_adult_(21_year)', 'Homo_sapiens_thoracic_aorta_endothelial_primary_cell_female_adult_(22_years)', 'Homo_sapiens_vein_endothelial_primary_cell_male_adult_(48_years)','Homo_sapiens_ovary_tissue_female_a
dult_(53_years)','Homo_sapiens_uterus_tissue_female_adult_(51_year)')

data1 <- read.table('Encode_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
data2 <- read.table('BingLim_EndodemInduction_CellStemCell_2014_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
Sox17 <- dcast(rbind(data1,data2), sample~gene)
colnames(Sox17) <- c('sample', 'gene', 'lnc')
Sox17$ID <- 'Sox17 - lncSox17'

data1 <- read.table('GATA6_Encode_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
data2 <- read.table('GATA6_BingLim_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
GATA6 <- dcast(rbind(data1,data2), sample~gene)
colnames(GATA6) <- c('sample', 'gene', 'lnc')
GATA6$ID <- 'Gata6 - Gata6AS'

data1 <- read.table('FOXA2_Encode_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
data2 <- read.table('FOXA2_BingLim_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
FOXA2 <- dcast(rbind(data1,data2), sample~gene)
colnames(FOXA2) <- c('sample', 'gene', 'lnc')
FOXA2$ID <- 'FOXA2 - DEANR1'

data1 <- read.table('NKX2_Encode_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
data2 <- read.table('NKX2_BingLim_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
NKX2 <- dcast(rbind(data1,data2), sample~gene)
colnames(NKX2) <- c('sample', 'gene', 'lnc')
NKX2$ID <- 'NKX2-1 - NKX2-1-AS1'

data1 <- read.table('GATA3_Encode_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
data2 <- read.table('GATA3_BingLim_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
GATA3 <- dcast(rbind(data1,data2), sample~gene)
colnames(GATA3) <- c('sample', 'gene', 'lnc')
GATA3$ID <- 'GATA3 - GATA3-AS1'

data1 <- read.table('LHX1_Encode_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
data2 <- read.table('LHX1_BingLim_TPM.txt', header=F, col.names=c('sample','gene','TPM'))
LHX1 <- dcast(rbind(data1,data2), sample~gene)
colnames(LHX1) <- c('sample', 'gene', 'lnc')
LHX1$ID <- 'LHX1 - LHX1-DT'

df <- rbind(Sox17, GATA6, FOXA2, NKX2, GATA3, LHX1)
df <- df[! df$sample %in% n, ]

p <- ggplot(df, aes(x=log1p(gene), y=log1p(lnc), color=ID)) + geom_point(size=2) + theme_classic() + geom_smooth(data=subset(df, sample!='Definitive_Endoderm_(DE)'), aes(x=log1p(gene), y=log1p(lnc), color
=ID), method='lm', formula= y~x) + geom_text(data=subset(df, sample=='Definitive_Endoderm_(DE)'), aes(x=log1p(gene), y=log1p(lnc)), label='DE', color='black', size=3, hjust = 0, nudge_x = 0.1) + ggtitle('
Linear model excluding DE')
ggMarginal(p1, type="boxplot", groupFill = TRUE)
```
