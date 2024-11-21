#setwd("EVE/refs/")

spp<-read.table(file="spp", header = F)

for(q in 1:nrow(spp)){
  
  library(dplyr)
  spec=spp[q,]
  prot<-list.files(path=spec, pattern = "protein.faa")
  
  df1<-read.table(file=paste0(spec,"/",spec,"_genes.tsv"), header = T, sep = "\t", comment.char = "", quote = "", na.strings = "")
  df2<-read.table(file=paste0(spec,"/",prot), header = F, sep = "\t", comment.char = "", quote = "")
  df3<-df1[!is.na(df1$Protein.accession),]
  
  df3 <- df3 %>%
    group_by(Protein.accession) %>%
    summarise(
      Begin = first(Begin),
      End = first(End),
      Chromosome = first(Chromosome),
      Orientation = first(Orientation),
      Name = first(Name),
      Symbol = first(Symbol),
      Gene.ID = first(Gene.ID),
      Gene.Type = first(Gene.Type),
      Transcripts.accession = first(Transcripts.accession),
      Protein.length = first(Protein.length),
      Locus.tag = first(Locus.tag),
      header = first(header)
    )
  
  df2$V1<-gsub(">","",df2$V1)
  df2$V1<-gsub("\\ .*","",df2$V1)
  
  df3$header<-paste0(">",spec,"_",df3$Protein.accession,"; ",df3$Name,"; gene:",df3$Symbol)
  
  df3<-df3 %>% group_by(header) %>% filter(duplicated(header) | n()==1)
  
  df2_updated <- df2 %>%
    left_join(df3, by = c("V1" = "Protein.accession")) %>%
    mutate(V1 = ifelse(!is.na(header), header, V1))
  
  df2_updated<-as.data.frame(df2_updated$V1)
  
  genes<-as.data.frame(df2_updated[grepl(">",df2_updated$`df2_updated$V1`),])
  genes[duplicated(genes$`df2_updated[grepl(">", df2_updated$\`df2_updated$V1\`), ]`),]
  
  write.table(file=paste0(spec,"/",spec,"_protein.fasta"),
              df2_updated, quote = F, row.names = F, col.names = F)
  }