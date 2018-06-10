#setwd("/srv/shiny-server/")
virus_id     <- paste0(getwd(),"/database/virus.id")
virus_db     <- paste0(getwd(),"/database/virus_db.fasta")
virus_msa_db <- paste0(getwd(),"/database/virus_msa.fasta")
orf_db       <- paste0(getwd(),"/database/orf_db.fasta")

mafft_path   <- "mafft"
blastdb_path <- "/ncbi-blast-2.7.1+/bin/makeblastdb"
blastp_path  <- "/ncbi-blast-2.7.1+/bin/blastp"
orf_path     <- "ORF/ORFfinder"

# mafft_path   <- "/Users/mayu/anaconda3/bin/mafft"
# blastdb_path <- "/Users/mayu/anaconda3/bin/makeblastdb"
# blastp_path  <- "/Users/mayu/anaconda3/bin/blastp"
# orf_path     <- "ORF/ORFfinder"