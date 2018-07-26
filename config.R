#setwd("/srv/shiny-server/")
virus_id     <- paste0(getwd(),"/database/virus.id")
virus_db     <- paste0(getwd(),"/database/virus_db.fasta")
virus_msa_db <- paste0(getwd(),"/database/virus_msa.fasta")
orf_db       <- paste0(getwd(),"/database/orf_db.fasta")
gb_id        <- paste0(getwd(),"/database/gb.txt") 

 mafft_path   <- "mafft"
 blastdb_path <- "/ncbi-blast-2.7.1+/bin/makeblastdb"
 blastp_path  <- "/ncbi-blast-2.7.1+/bin/blastp"
 orf_path     <- "ORF/ORFfinder"

#mafft_path   <- "/Users/mayu/anaconda3/bin/mafft"
#blastdb_path <- "/Users/mayu/anaconda3/bin/makeblastdb"
#blastp_path  <- "/Users/mayu/anaconda3/bin/blastp"
#orf_path     <- "ORF/ORFfinder"


orfViewer <- function(x){
  orf_asn <- read.table(paste0("ORF/",x,".asn"),sep = "\n",stringsAsFactors = F)
  orf_label <- gsub(",","",gsub(".*str ","",orf_asn[grep("product",orf_asn$V1),]))
  orf_start  <- na.omit(as.numeric(gsub(",","",gsub(".*from ","",orf_asn[grep("from",orf_asn$V1),]))))
  orf_end    <- na.omit(as.numeric(gsub(",","",gsub(".*to ","",orf_asn[grep("to",orf_asn$V1),]))))
  orf_strand <- gsub(",","",gsub(".*strand ","",orf_asn[grep("strand",orf_asn$V1),]))
  orf_length <- na.omit(as.numeric(gsub(",","",gsub(".*length ","",orf_asn[grep("length",orf_asn$V1),]))))
  orf_sequence <- gsub("\n","",gsub(".*iupacaa ","",orf_asn[grep("iupacaa",orf_asn$V1),]))
  
  orf_table <- data.frame(orf_label,orf_strand,orf_start,orf_end,orf_length,orf_sequence,stringsAsFactors = F)
  orf_table$orf_strand[orf_table$orf_strand == "plus"] <- "+"
  orf_table$orf_strand[orf_table$orf_strand == "minus"] <- "-"
  orf_table$nucleo_length <- abs(orf_table$orf_start - orf_table$orf_end) + 1
  write.table(orf_table,paste0("ORF/",x,".table"),sep = "\t",row.names = F,quote = F)
  orf_model <- orf_table[,c(1,3,4,7,2)]
  orf_model$chromosome <- "virus"
  orf_model$orf_label <- paste0("ORF",1:nrow(orf_table))
  orf_model <- orf_model[,c(6,2,3,4,5,1)]
  colnames(orf_model) <- c("chromosome","start","end","width","strand","symbol")
  orf_model$transcript <- orf_model$symbol
  orf_model$gene <- orf_model$symbol
  orf_model$start <- as.numeric(orf_model$start)
  orf_model$end <- as.numeric(orf_model$end)
  orf_model$width <- as.numeric(orf_model$width)
  orf_model$feature <- "protein_coding"
  write.table(orf_model,paste0("ORF/",x,".viewer"),sep = "\t",row.names = F,quote = F)
}
