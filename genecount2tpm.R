##############################
### Gene length extraction ###
##############################

#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)

all<-read.table(file="EVE/all_comb_htseq_ordered.txt", header = T)

##########
# bosTau #
##########
txdb <- makeTxDbFromGFF("EVE/gtf/bosTau.gtf", format="gtf")

exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths <- as.data.frame(sum(width(reduce(exons_by_gene))))

head(gene_lengths)
gene_lengths$Gene<-rownames(gene_lengths)
colnames(gene_lengths)<-c("Length","Gene")

bT<-as.data.frame(cbind(all$OrthoGroup, all$bosTau_GeneID,
                        all$bosTau_SRR17288213, all$bosTau_SRR17288214, all$bosTau_SRR17288216))
colnames(bT)<-c("OrthoGroup","Gene_ID","bosTau_SRR17288213","bosTau_SRR17288214","bosTau_SRR17288216")
bT2<-merge(bT, gene_lengths, by.x = "Gene_ID", by.y = "Gene")

#S1
bT2$bosTau_SRR17288213_rpk<-as.numeric(bT2$bosTau_SRR17288213)/as.numeric(bT2$Length/1000)
s1_sf<-sum(bT2$bosTau_SRR17288213_rpk)
bT2$bosTau_SRR17288213_tpm <- (bT2$bosTau_SRR17288213_rpk / s1_sf) * 1e6

#S2
bT2$bosTau_SRR17288214_rpk<-as.numeric(bT2$bosTau_SRR17288214)/as.numeric(bT2$Length/1000)
s1_sf<-sum(bT2$bosTau_SRR17288214_rpk)
bT2$bosTau_SRR17288214_tpm <- (bT2$bosTau_SRR17288214_rpk / s1_sf) * 1e6

#S3
bT2$bosTau_SRR17288216_rpk<-as.numeric(bT2$bosTau_SRR17288216)/as.numeric(bT2$Length/1000)
s1_sf<-sum(bT2$bosTau_SRR17288216_rpk)
bT2$bosTau_SRR17288216_tpm <- (bT2$bosTau_SRR17288216_rpk / s1_sf) * 1e6

bT_tpm<-as.data.frame(cbind(bT2$Gene_ID, bT2$OrthoGroup,
                            bT2$bosTau_SRR17288213_tpm, bT2$bosTau_SRR17288214_tpm, bT2$bosTau_SRR17288216_tpm))
colnames(bT_tpm)<-c("Gene_ID","OrthoGroup","bosTau_SRR17288213_TPM","bosTau_SRR17288214_TPM","bosTau_SRR17288216_TPM")

##########
# musMus #
##########
txdb <- makeTxDbFromGFF("EVE/gtf/musMus.gtf", format="gtf")

exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths <- as.data.frame(sum(width(reduce(exons_by_gene))))

head(gene_lengths)
gene_lengths$Gene<-rownames(gene_lengths)
colnames(gene_lengths)<-c("Length","Gene")

mM<-as.data.frame(cbind(all$OrthoGroup, all$musMus_GeneID,
                        all$musMus_SRR24530105, all$musMus_SRR24530106, all$musMus_SRR24530107))
colnames(mM)<-c("OrthoGroup","Gene_ID","musMus_SRR24530105","musMus_SRR24530106","musMus_SRR24530107")
mM2<-merge(mM, gene_lengths, by.x = "Gene_ID", by.y = "Gene")

#S1
mM2$musMus_SRR24530105_rpk<-as.numeric(mM2$musMus_SRR24530105)/as.numeric(mM2$Length/1000)
s1_sf<-sum(mM2$musMus_SRR24530105_rpk)
mM2$musMus_SRR24530105_tpm <- (mM2$musMus_SRR24530105_rpk / s1_sf) * 1e6

#S2
mM2$musMus_SRR24530106_rpk<-as.numeric(mM2$musMus_SRR24530106)/as.numeric(mM2$Length/1000)
s1_sf<-sum(mM2$musMus_SRR24530106_rpk)
mM2$musMus_SRR24530106_tpm <- (mM2$musMus_SRR24530106_rpk / s1_sf) * 1e6

#S3
mM2$musMus_SRR24530107_rpk<-as.numeric(mM2$musMus_SRR24530107)/as.numeric(mM2$Length/1000)
s1_sf<-sum(mM2$musMus_SRR24530107_rpk)
mM2$musMus_SRR24530107_tpm <- (mM2$musMus_SRR24530107_rpk / s1_sf) * 1e6

mM_tpm<-as.data.frame(cbind(mM2$Gene_ID, mM2$OrthoGroup,
                            mM2$musMus_SRR24530105_tpm, mM2$musMus_SRR24530106_tpm, mM2$musMus_SRR24530107_tpm))
colnames(mM_tpm)<-c("Gene_ID","OrthoGroup","musMus_SRR24530105_TPM","musMus_SRR24530106_TPM","musMus_SRR24530107_TPM")


##########
# myoMyo #
##########
txdb <- makeTxDbFromGFF("EVE/gtf/myoMyo.gtf", format="gtf")

exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths <- as.data.frame(sum(width(reduce(exons_by_gene))))

head(gene_lengths)
gene_lengths$Gene<-rownames(gene_lengths)
colnames(gene_lengths)<-c("Length","Gene")

mY<-as.data.frame(cbind(all$OrthoGroup, all$myoMyo_GeneID,
                        all$myoMyo_fibro_MMY6321_1.1_S16, all$myoMyo_fibro_MMY6321_1.2_S17, all$myoMyo_fibro_MMY6321_1.3_S18))
colnames(mY)<-c("OrthoGroup","Gene_ID","myoMyo_fibro_MMY6321_1.1_S16","myoMyo_fibro_MMY6321_1.2_S17","myoMyo_fibro_MMY6321_1.3_S18")
mY2<-merge(mY, gene_lengths, by.x = "Gene_ID", by.y = "Gene")

#S1
mY2$myoMyo_fibro_MMY6321_1.1_S16_rpk<-as.numeric(mY2$`myoMyo_fibro_MMY6321_1.1_S16`)/as.numeric(mY2$Length/1000)
s1_sf<-sum(mY2$myoMyo_fibro_MMY6321_1.1_S16_rpk)
mY2$myoMyo_fibro_MMY6321_1.1_S16_tpm <- (mY2$myoMyo_fibro_MMY6321_1.1_S16_rpk / s1_sf) * 1e6

#S2
mY2$myoMyo_fibro_MMY6321_1.2_S17_rpk<-as.numeric(mY2$myoMyo_fibro_MMY6321_1.2_S17)/as.numeric(mY2$Length/1000)
s1_sf<-sum(mY2$myoMyo_fibro_MMY6321_1.2_S17_rpk)
mY2$myoMyo_fibro_MMY6321_1.2_S17_tpm <- (mY2$myoMyo_fibro_MMY6321_1.2_S17_rpk / s1_sf) * 1e6

#S3
mY2$myoMyo_fibro_MMY6321_1.3_S18_rpk<-as.numeric(mY2$myoMyo_fibro_MMY6321_1.3_S18)/as.numeric(mY2$Length/1000)
s1_sf<-sum(mY2$myoMyo_fibro_MMY6321_1.3_S18_rpk)
mY2$myoMyo_fibro_MMY6321_1.3_S18_tpm <- (mY2$myoMyo_fibro_MMY6321_1.3_S18_rpk / s1_sf) * 1e6

mY_tpm<-as.data.frame(cbind(mY2$Gene_ID, mY2$OrthoGroup,
                            mY2$myoMyo_fibro_MMY6321_1.1_S16_tpm, mY2$myoMyo_fibro_MMY6321_1.2_S17_tpm, mY2$myoMyo_fibro_MMY6321_1.3_S18_tpm))
colnames(mY_tpm)<-c("Gene_ID","OrthoGroup","myoMyo_fibro_MMY6321_1.1_S16_TPM","myoMyo_fibro_MMY6321_1.2_S17_TPM","myoMyo_fibro_MMY6321_1.3_S18_TPM")

##########
# susScr #
##########
txdb <- makeTxDbFromGFF("EVE/gtf/susScr.gtf", format="gtf")

exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths <- as.data.frame(sum(width(reduce(exons_by_gene))))

head(gene_lengths)
gene_lengths$Gene<-rownames(gene_lengths)
colnames(gene_lengths)<-c("Length","Gene")

sS<-as.data.frame(cbind(all$OrthoGroup, all$susScr_GeneID,
                        all$susScr_SRR13844546, all$susScr_SRR13844547))
colnames(sS)<-c("OrthoGroup","Gene_ID","susScr_SRR13844546","susScr_SRR13844547")
sS2<-merge(sS, gene_lengths, by.x = "Gene_ID", by.y = "Gene")

#S1
sS2$susScr_SRR13844546_rpk<-as.numeric(sS2$susScr_SRR13844546)/as.numeric(sS2$Length/1000)
s1_sf<-sum(sS2$susScr_SRR13844546_rpk)
sS2$susScr_SRR13844546_tpm <- (sS2$susScr_SRR13844546_rpk / s1_sf) * 1e6

#S2
sS2$susScr_SRR13844547_rpk<-as.numeric(sS2$susScr_SRR13844547)/as.numeric(sS2$Length/1000)
s1_sf<-sum(sS2$susScr_SRR13844547_rpk)
sS2$susScr_SRR13844547_tpm <- (sS2$susScr_SRR13844547_rpk / s1_sf) * 1e6

sS_tpm<-as.data.frame(cbind(sS2$Gene_ID, sS2$OrthoGroup,
                            sS2$susScr_SRR13844546_tpm, sS2$susScr_SRR13844547_tpm))
colnames(sS_tpm)<-c("Gene_ID","OrthoGroup","susScr_SRR13844546_TPM","susScr_SRR13844547_TPM")

all_tpm<- bT_tpm %>% left_join(mM_tpm, by = "OrthoGroup") %>%
  left_join(mY_tpm, by = "OrthoGroup") %>%
  left_join(sS_tpm, by = "OrthoGroup")
write.table(file="EVE/all_TPM.txt", all_tpm, quote = F, row.names = F, col.names = T, sep = "\t")
