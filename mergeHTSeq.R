library(purrr)
library(dplyr)
library(GenomicFeatures)

spp<-read.table(file="EVE/spp", header = T)
keys<-read.table(file="EVE/combined_keys.txt", header = T)

for(q in 1:nrow(spp)){
  spec=spp[q,]
  spec_key<-keys[,c(1,q+1)]
  
  file_list<-paste0("EVE/htseq/",spec,"/",list.files(path = paste0("EVE/htseq/",spec)))
  samples<-paste0(spec,"_",gsub("_GeneCount.txt","",list.files(path = paste0("EVE/htseq/",spec))))
  file_df<-lapply(file_list, read.table)
  merged_df <- file_df %>%
    purrr::reduce(inner_join, by = "V1")
  colnames(merged_df)<-c("Gene",samples)
  merged_df$Gene<-paste0(spec, "_", merged_df$Gene)
  
  OG_DF<-merge(spec_key, merged_df, by.x = spec, by.y = "Gene")
  
  gtf_file<-paste0("EVE/refs/",spec,"/",list.files(path = paste0("EVE/refs/",spec)))
  txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
  
  exons_by_gene <- exonsBy(txdb, by = "gene")
  gene_lengths <- as.data.frame(sum(width(reduce(exons_by_gene))))
  
  head(gene_lengths)
  
  gene_lengths$Gene<-rownames(gene_lengths)
  gene_lengths <- gene_lengths %>%
    mutate_all(~ if(is.character(.)) toupper(.) else .)
  gene_lengths$Gene<-paste0(spec,"_",gene_lengths$Gene)
  colnames(gene_lengths)<-c("Length","Gene")
  
  OG_GL<-merge(OG_DF, gene_lengths, by.x = spec, by.y = "Gene")
  
  for(x in 1:length(samples)){
    OG_GL[,x+2]<-as.numeric(OG_GL[,x+2])/as.numeric(OG_GL$Length/1000)
    sf<-sum(as.numeric(OG_GL[,x+2]))
    OG_GL[,x+2] <- (OG_GL[,x+2] / sf) * 1e6
  }
  OG_GL$Length<-NULL
  OG_GL[,1]<-NULL
  write.table(file=paste0("EVE/tpm/",spec,"_TPM.txt"), OG_GL, quote = F, row.names = F, sep = "\t")
}


all_list<-paste0("EVE/TPM/",list.files(path="EVE/TPM"))
all_df<-lapply(all_list, function(file) read.table(file, header = TRUE))
all_spp_merged <- all_df %>%
  purrr::reduce(inner_join, by = "OG")

homSap_key<-as.data.frame(cbind(keys$OG, keys$homSap))
homSap_key$V2<-gsub("homSap_","",homSap_key$V2)

final_df<-merge(homSap_key, all_spp_merged, by.x = "V1", by.y = "OG")
colnames(final_df)[c(1,2)]<-c("OG", "Gene_ID")

write.table(file="EVE/TPM/all_species_TPM.txt", final_df, row.names = F, quote = F, sep = "\t")

eve_input<-final_df
eve_input$OG<-NULL
eve_input$bosTau_SRR17288215<-NULL
write.table(file="EVE/TPM/EVE_input_2024.10.15.txt", eve_input, row.names = F, col.names = F, quote = F, sep = "\t")

# Re-structure to fit with the updated tree
df<-as.data.frame(cbind(eve_input$Gene_ID, eve_input$musMus_SRR24530105, eve_input$musMus_SRR24530106, eve_input$musMus_SRR24530107,
                        eve_input$homSap_SRR23630190, eve_input$homSap_SRR23630193, eve_input$homSap_SRR23630200, eve_input$homSap_SRR23630208,
                        eve_input$rouAeg_SAMEA9699320, eve_input$rouAeg_SAMEA9699322, eve_input$rouAeg_SAMEA9699328, eve_input$rouAeg_SAMEA9699331,
                        eve_input$rhiFer_SRR22938632, eve_input$rhiFer_SRR22938633, eve_input$rhiFer_SRR22938634,
                        eve_input$pipKuh_SAMEA9699324, eve_input$pipKuh_SAMEA9699326, eve_input$pipKuh_SAMEA9699330, eve_input$pipKuh_SAMEA9699342,
                        eve_input$myoMyo_fibro_MMY6321_1.1_S16, eve_input$myoMyo_fibro_MMY6321_1.2_S17, eve_input$myoMyo_fibro_MMY6321_1.3_S18,
                        eve_input$canFam_SRR17870687, eve_input$canFam_SRR17870688, eve_input$canFam_SRR17870691, eve_input$canFam_SRR17870692,
                        eve_input$felCat_SRR24448362, eve_input$felCat_SRR24448363, eve_input$felCat_SRR24448364,
                        eve_input$equCab_SRR18906497, eve_input$equCab_SRR18906498, eve_input$equCab_SRR18906499, eve_input$equCab_SRR18906500,
                        eve_input$bosTau_SRR17288213, eve_input$bosTau_SRR17288214, eve_input$bosTau_SRR17288216,
                        eve_input$susScr_SRR25461976, eve_input$susScr_SRR25461977, eve_input$susScr_SRR25461978))

row.names(df)<-df$V1
df$V1<-NULL
df_filtered <- df %>%
  mutate_all(as.numeric) %>%
  select(where(~ sum(.) != 0))

write.table(file="EVE/TPM/EVE_input_restructured.txt", df_filtered, row.names = T, col.names = F, quote = F, sep = "\t")
