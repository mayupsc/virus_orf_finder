library(shiny)
library(data.table)
library(jsonlite)
library(shinydashboard)
library(seqinr)
library(DECIPHER)
source("config.R")

server <-shinyServer(function(input, output) {

  

  output$uploadfiles <- renderTable({

    inFile <- input$file
    if (is.null(inFile)) return(NULL)
    files <- unzip(inFile$datapath)
    gsub("./.*/","",files)

  })

  output$intervirus <- renderUI({
    inFile <- input$file
    files <- unzip(inFile$datapath)
    selectInput("intervirus", "Interest Virus:",
                choices = gsub(".fasta","",gsub("./.*/","",files)),
                selected = gsub(".fasta","",gsub("./.*/","",files))[7]
    )
  })

  output$interorf <- renderUI({
    blast_json <- fromJSON(paste0("blastp/",input$intervirus,".blast"))
    selectInput("interorf", "Interest ORF:",
                choices = paste0("ORF",1:length(blast_json$BlastOutput2$report$program),"_",input$intervirus),
                selected =paste0("ORF",1:length(blast_json$BlastOutput2$report$program),"_",input$intervirus)[1]
    )
  })

  output$plot1 <- renderPlot({
      Sys.sleep(2);
      inFile <- input$file
      if (is.null(inFile)) return(NULL)
      files <- unzip(inFile$datapath)
      system(paste("cat",paste0(files,collapse = " "),">",virus_db))
      system(paste0("grep '>' ",virus_db," | sed 's/>//g' > ",virus_id))
      system(paste0(mafft_path," ",virus_db," > ",virus_msa_db))
      #########dend###########
      dbConn <-dbConnect(SQLite(),":memory:")
      Seqs2DB(virus_msa_db,type = "FASTA",dbFile = dbConn,identifier = "")
      query_desc<- dbGetQuery(dbConn,"select description from Seqs")$description
      Add2DB(myData = data.frame(identifier=query_desc,stringsAsFactors = FALSE),dbFile = dbConn)
      query_desc<-unlist(lapply(strsplit(query_desc,"curl"),'[',2L))
      query_desc<-paste0("TYLC",query_desc)
      Add2DB(myData = data.frame(identifier=query_desc,stringsAsFactors = FALSE),dbConn)
      cons <- IdConsensus(dbConn,threshold=0.3,minInformation=0.1)
      d <- DistanceMatrix(cons,correction = "Jukes-Cantor")
      dend<-IdClusters(d,method="ML",type = "dendrogram",myXStringSet=cons)
      p<- par(mar=c(1,1,1,30),xpd=TRUE)
      plot(dend,yaxt="n",horiz = TRUE)
      par(p)

  })
  output$orfdb <- renderTable({
    inFile <- input$file
    files <- unzip(inFile$datapath)
    virus_gb <- c()
    orf_number <- c()
    withProgress("Calculating Progress",value = 0,min=0,max = length(files),{
    for (i in 1:length(files)){
      input$refresh
      incProgress(1, detail = paste("Analyzing virus ", i))
      gb <-gsub(".fasta","",gsub("./.*/","",files[i]))
      orf_cmd <- paste0(orf_path," -id ",gb," -ml ",input$minorflen," -s 0 -out ",paste0("ORF/",gb,".orfs"))
      system(orf_cmd)
      if (file.info(paste0("ORF/",gb,".orfs"))$size != 0){
        virus <- read.fasta(paste0("ORF/",gb,".orfs"))
        virus_gb <-c(gb,virus_gb)
        orf_number <- c(length(virus),orf_number)
      }

    }
    })
    system(paste0("cat ORF/*.orfs > ",orf_db))
    system(paste0(blastdb_path," -in ",orf_db," -parse_seqids -dbtype prot"))

    orfdb <- data.frame(virus_gb,orf_number)

  })

  output$blastp <- renderTable({
    system(paste0(blastp_path," -query ORF/",input$intervirus,".orfs -db ",orf_db," -outfmt 15"," -evalue ",input$evalue," > blastp/",input$intervirus,".blast"))
    blast_json <- fromJSON(paste0("blastp/",input$intervirus,".blast"))
    hits_number <- c()
    for (i in 1:length(blast_json$BlastOutput2$report$program)) {
      hits_number<- c(dim(blast_json$BlastOutput2$report$results$search$hits[[i]])[1],hits_number)
    }
    blast_result <- data.frame("query_title"=blast_json[["BlastOutput2"]][["report"]][["results"]][["search"]][["query_title"]],
                               "query_length"=blast_json[["BlastOutput2"]][["report"]][["results"]][["search"]][["query_len"]],
                               "hits_number"= hits_number
                               )

  })

  output$cal <- renderTable({

    virus_id <- fread(virus_id)
    virus_id <- data.frame("gi"= virus_id$V2,"gb"=virus_id$V4,"title"=virus_id$V5)
    blast_json <- fromJSON(paste0("blastp/",input$intervirus,".blast"))
    dbConn <-dbConnect(SQLite(),":memory:")
    Seqs2DB(virus_msa_db,type = "FASTA",dbFile = dbConn,identifier = "")
    query_desc<- dbGetQuery(dbConn,"select description from Seqs")$description
    Add2DB(myData = data.frame(identifier=query_desc,stringsAsFactors = FALSE),dbFile = dbConn)
    query_desc<-unlist(lapply(strsplit(query_desc,"curl"),'[',2L))
    query_desc<-paste0("TYLC",query_desc)
    Add2DB(myData = data.frame(identifier=query_desc,stringsAsFactors = FALSE),dbConn)
    cons <- IdConsensus(dbConn,threshold=0.3,minInformation=0.1)
    d <- DistanceMatrix(cons,correction = "Jukes-Cantor")
    dend<-IdClusters(d,method="ML",type = "dendrogram",myXStringSet=cons)
    withProgress("Calculating Progress",value = 0,min=0,max = dim(blast_json$BlastOutput2)[1],{
    
    for (i in 1:dim(blast_json$BlastOutput2)[1]) {
      json_abb <- blast_json$BlastOutput2$report$results$search$hits[[i]]
      incProgress(1, detail = paste("Doing part", i))
      pdf(paste0("www/ORF",i,"_",input$intervirus,".pdf"),width = 11,height = 8.5)

      if (length(json_abb)>1){
        file_name <- paste0("ORF",i,"_",input$intervirus)
        orf_hit_desc<-c()
        system(paste0("touch output/msa/",file_name,".fasta"))

        for (j in 1:dim(json_abb)[1]){
          id <- json_abb$description[[j]]$id
          n_gap_fill <- paste0(rep("-",json_abb$hsps[[j]]$query_from-1),collapse = "")
          query <- paste0(n_gap_fill,json_abb$hsps[[j]]$hseq)
          query_title <- paste0(">",id)
          msa_query_title_cmd <- paste0("echo '",query_title,"' >> ",paste0("output/msa/",file_name,".fasta"))
          msa_query_cmd <- paste0("echo '",query,"' >> ",paste0("output/msa/",file_name,".fasta"))
          system(msa_query_title_cmd)
          system(msa_query_cmd)
          orf_hit_desc_pre <- substring(strsplit(id,":")[[1]][1],gregexpr("_",id)[[1]][1]+1)
          orf_hit_trans <- as.character(virus_id$title[which(orf_hit_desc_pre == virus_id$gb)])
          orf_hit_trans_split <- paste0("TYLC",strsplit(orf_hit_trans,"curl")[[1]][2])
          orf_hit_desc <- c(orf_hit_desc,orf_hit_trans_split)

        }
        local({
          colLab <<- function(n){
            if(is.leaf(n)){
              a <- attributes(n)
              i <<- i+1
              if(a$label %in% orf_hit_desc){
                attr(n,"nodePar") <- c(a$nodePar, list(lab.col = "firebrick3"))
              }

            }
            n
          }
          i <-0
        }
        )
        dL <- dendrapply(dend, colLab)
        p<- par(mar=c(1,1,1,30),xpd=TRUE)
        plot(dL,yaxt="n",horiz = TRUE,main =file_name)
        par(p)
        dev.off()
        dbConn <-dbConnect(SQLite(),":memory:")
        Seqs2DB(paste0("output/msa/",file_name,".fasta"),type = "FASTA",dbFile = dbConn,identifier = "")
        x <- dbGetQuery(dbConn,"select description from Seqs")$description
        Add2DB(myData = data.frame(identifier=x,stringsAsFactors = FALSE),dbFile = dbConn)
        AA <- SearchDB(dbConn,nameBy = "identifier")
        BrowseSeqs(AA,htmlFile = paste0(getwd(),"/",paste0("output/msa/",file_name,".html")),openURL = F)
        dbDisconnect(dbConn)
        system(paste0("sed -i 's/<\\/head>/<hl>",file_name,"<\\/hl><\\/head>/g' output/msa/",file_name,".html"))
        Sys.sleep(1)
      }

    }
     
      })
    system(paste0("cat ",paste0("output/msa/ORF",1:dim(blast_json$BlastOutput2)[1],"_",input$intervirus,".html",collapse = " ")," > www/",input$intervirus,".html"))
    system(paste0("pdftk ",paste0("www/ORF",1:dim(blast_json$BlastOutput2)[1],"_",input$intervirus,".pdf",collapse = " "),"cat output"," > www/",input$intervirus,".pdf"))

  })

  output$msa <- renderUI({
    return(includeHTML(paste0("output/msa/",input$interorf,".html")))
  })

  output$phylo <- renderUI({
    tags$iframe(style="height:600px; width:100%", src=paste0(input$interorf,".pdf"))
  })
  
  output$downloadPDF <- downloadHandler(
    filename = paste0(input$intervirus,".pdf"),
    content = function(file){
      file.copy(paste0("www/",input$intervirus,".pdf"),file)
    }
  )

  output$downloadHTML <- downloadHandler(
    filename = paste0(input$intervirus,".html"),
    content = function(file){
      file.copy(paste0("www/",input$intervirus,".html"),file)
    }
  )



})
