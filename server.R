library(shiny)
library(data.table)
library(jsonlite)
library(shinydashboard)
library(seqinr)
library(DECIPHER)
library(DT)
library(Gviz)
source("config.R")

server <-shinyServer(function(input, output) {

  
  
  output$preparation <- renderPrint({
    cat("extract files...\n")
    inFile <- input$file
    if (is.null(inFile)) return(NULL)
    files <- unzip(inFile$datapath)
    cat("extract genbank numbers...\n")
    Species <- gsub(".fasta","",gsub("./.*/","",files))
    write.table(Species,gb_id,quote = F,col.names = F,row.names = F )
    
    ## virus_db = "virus_db.fasta"
    cat("build virus genome sequence database...\n")
    system(paste("cat",paste0(files,collapse = " "),">",virus_db))
    
    ## virus_id = "virus.id"
    cat("extract virus genome description...\n")
    system(paste0("grep '>' ",virus_db," | sed 's/>//g' > ",virus_id))
    
    ## input file is virus_db.fasta,output file is virus_msa.fasta 
    cat("multiple alignment virus genome sequence...\n")
    system(paste0(mafft_path," ",virus_db," > ",virus_msa_db))
    
  })
  
  # plot raw phylo tree
  output$plot1 <- renderPlot({
    

    #########code from DECIPHER to plot phylo tree###########
    input$goButton
    dbConn <- dbConnect(SQLite(),":memory:")
    Seqs2DB(virus_msa_db,type = "FASTA",dbFile = dbConn,identifier = "")
    query_desc <- dbGetQuery(dbConn,"select description from Seqs")$description
    # Add2DB(myData = data.frame(identifier=query_desc,stringsAsFactors = FALSE),dbFile = dbConn)
    
    ## replace "Tomato yellow leaf curl" by "TYLC" to simplify identifiers
    query_desc<-gsub("gi\\|.*\\|.*\\|.*\\| ","",query_desc) #simplify tree identifier
    query_desc<-gsub(input$raw_name,input$new_name,query_desc)
    
    Add2DB(myData = data.frame(identifier=query_desc,stringsAsFactors = FALSE),dbConn)
    cons <- IdConsensus(dbConn,threshold=0.3,minInformation=0.1)
    d <- DistanceMatrix(cons,correction = "Jukes-Cantor")
    dend <- IdClusters(d,method="ML",type = "dendrogram",myXStringSet=cons)
    p <- par(mar=c(1,1,1,30),xpd=TRUE)
    plot(dend,yaxt="n",horiz = TRUE)
    par(p)
    
  })

# upload zip files
  output$orfdatabase <- DT::renderDataTable({
    
    #################make orf database##########################
    
    input$refresh
    #system(paste0("cat ",gb_id,"| parallel --verbose --eta ",orf_path," -id {} -ml ",input$minorflen," -s 0 -out ORF/{}.orfs"))
    #system(paste0("cat ",gb_id,"| parallel --verbose --eta ",orf_path," -id {} -ml ",input$minorflen," -s 0 -outfmt 2 -out ORF/{}.asn"))
    ## orf_db ="orf_db.fasta"
    system(paste0("cat ORF/*.orfs > ",orf_db))
    ## makeblastdb, input file is "orf_db.fasta"
    system(paste0(blastdb_path," -in ",orf_db," -parse_seqids -dbtype prot"))
   
    Species <- read.table(gb_id,stringsAsFactors = F)
    virus_length <- c()
    orf_num <- c()
    for (i in 1:dim(Species)[1]){
      virus <- read.fasta(files[i])
      orf <- read.fasta(paste0("ORF/",Species$V1[i],".orfs"))
      orf_num <- c(orf_num,length(orf))
      virus_length <- c(virus_length,length(virus[[1]]))
    }
    orfdb <- data.frame(Species,virus_length,orf_num,stringsAsFactors = F)
    colnames(orfdb)[1] <- "Species"
    write.table(orfdb,"database/virus_length.txt",sep = "\t",quote = F,row.names = F)
    datatable(data.frame(Species=paste0("<a href='#blastp'",
                                        "alt='",orfdb$Species,"'",                                                 
                                        "onclick=\"",
                                        "tabs = $('.tabbable .nav.nav-tabs li');",
                                        "tabs.each(function() {",
                                        "$(this).removeClass('active')",
                                        "});",
                                        "$(tabs[1]).addClass('active');",
                                        "tabsContents = $('.tabbable .tab-content .tab-pane');",
                                        "tabsContents.each(function() {",
                                        "$(this).removeClass('active')",
                                        "});",
                                        "$(tabsContents[1]).addClass('active');",
                                        "$('#blastp').trigger('change').trigger('shown');",
                                        "Shiny.onInputChange('species', getAttribute('alt'));",
                                        "\">",
                                        orfdb$Species,
                                        "</a>"),orfdb$virus_length,orfdb$orf_num),
              escape = FALSE)
  })
  
  
  output$blastp <- DT::renderDataTable({
    server = F
    if(is.null(input$species)){
      datatable(orfdb)
    }else{
      #outfmt = json 
      system(paste0(blastp_path," -query ORF/",input$species,".orfs -db ",orf_db," -outfmt 15 -evalue ",input$evalue," > blastp/",input$species,".blast"))
      blast_json <- fromJSON(paste0("blastp/",input$species,".blast"))
      hits_number <- c()
      for (i in 1:length(blast_json$BlastOutput2$report$program)) {
        hits_number<- c(hits_number,dim(blast_json$BlastOutput2$report$results$search$hits[[i]])[1])
      }

      orfViewer(input$species)
      orf_table <- read.table(paste0("ORF/",input$species,".table"),sep = "\t",header = T,stringsAsFactors = F)
      blast_result <- data.frame(cbind(orf_table),hits_number)

      datatable(
        cbind(' ' = '&oplus;', blast_result), escape = -0,
        options = list(
          columnDefs = list(
            list(visible = FALSE, targets = c(7)),
            list(orderable = FALSE, className = 'details-control', targets = 1)
          )
        ),
        callback = JS("
                      table.column(1).nodes().to$().css({cursor: 'pointer'});
                      var format = function(d) {
                      return '<div style=\"background-color:#eee; padding: .5em;\"> orf_sequence ' +
                              d[7]  + '</div>';
                      };
                      table.on('click', 'td.details-control', function() {
                      var td = $(this), row = table.row(td.closest('tr'));
                      if (row.child.isShown()) {
                      row.child.hide();
                      td.html('&oplus;');
                      } else {
                      row.child(format(row.data())).show();
                      td.html('&CircleMinus;');
                      }
                      });"
))
      
    }
  })
  
  #plot ORF Viewer
  output$plot2 <- renderPlot({
    orfdb <- read.table("database/virus_length.txt",sep = "\t",header = T)
    orf_model <- read.table(paste0("ORF/",input$species,".viewer"),sep = "\t",header = T,stringsAsFactors = F)
    axisTrack <- GenomeAxisTrack(range = IRanges(start=1,end = orfdb$virus_length[orfdb$Species==input$species],names="virus genome"))
    aTrack <- AnnotationTrack(start = orf_model$start,width = orf_model$width, chromosome = "chrX", strand = c(orf_model$strand), id = orf_model$symbol,genome = "hg19", name = "ORF Viewer",transcriptAnnotation="symbol")
    plotTracks(list(axisTrack,aTrack), from = 1, to = orfdb$virus_length[orfdb$Species==input$species],showId=F,fontcolor="black",add53=T,add35=T,cex = 1,labelPos="above",cex.id=1.5,col.id="black",main = input$species,col=NULL,featureAnnotation="id",fontcolor.feature = "black")
    
  })


  #show alignment of interest ORF
  output$msa <- renderUI({
    interorf <- input$blastp_rows_selected
    blast_json <- fromJSON(paste0("blastp/",input$species,".blast"))
      for (i in interorf) {
        json_abb <- blast_json$BlastOutput2$report$results$search$hits[[i]] ###abbreviate code
        ### pdf phylo tree for each orf in interest virus
        file_name <- paste0("ORF",i,"_",input$species)
        #file_name <- paste0("ORF",i,"_","AJ489258.1")
        orf_hit_desc<-c()
        
        ####write ORFxxxx.fasta,which comprising querys that hits ORFxxxxx,
        ####ORFxxxx.fasta is the basic of multiple sequence alignment
        system(paste0("echo '' > output/msa/",file_name,".fasta"))
        for (j in 1:dim(json_abb)[1]){
          id <- json_abb$description[[j]]$id
          
          #####add space or query showing in a wrong way
          n_gap_fill <- paste0(rep("-",json_abb$hsps[[j]]$query_from-1),collapse = "")
          query <- paste0(n_gap_fill,json_abb$hsps[[j]]$hseq)
          query_title <- paste0(">",id)
          
          #####echo processed query to ORFxxxxx.fasta
          system(paste0("echo '",query_title,"' >> ",paste0("output/msa/",file_name,".fasta")))
          system(paste0("echo '",query,"' >> ",paste0("output/msa/",file_name,".fasta")))
        
        ####show alignment
        dbConn <-dbConnect(SQLite(),":memory:")
        Seqs2DB(paste0("output/msa/",file_name,".fasta"),type = "FASTA",dbFile = dbConn,identifier = "")
        x <- dbGetQuery(dbConn,"select description from Seqs")$description
        Add2DB(myData = data.frame(identifier=x,stringsAsFactors = FALSE),dbFile = dbConn)
        AA <- SearchDB(dbConn,nameBy = "identifier")
        BrowseSeqs(AA,htmlFile = paste0(getwd(),"/",paste0("output/msa/",file_name,".html")),openURL = F)
        dbDisconnect(dbConn)
          } 
        system(paste0("sed -i 's/<\\/head>/<hl>",file_name,"<\\/hl><\\/head>/g' output/msa/",file_name,".html"))
        Sys.sleep(1)
      }
    
    
   
    ##merge alignment results
    system(paste0("cat ",paste0("output/msa/ORF",interorf,"_",input$species,".html",collapse = " ")," > www/",input$species,".html"))
    return(includeHTML(paste0("www/",input$species,".html")))
      
  })
    
  
  output$phylo <- renderUI({
    
    interorf <- input$blastp_rows_selected
    virus_id <- fread(virus_id)
    virus_id <- data.frame("gi"= virus_id$V2,"gb"=virus_id$V4,"title"=virus_id$V5)
    blast_json <- fromJSON(paste0("blastp/",input$species,".blast"))
    dbConn <-dbConnect(SQLite(),":memory:")
    Seqs2DB(virus_msa_db,type = "FASTA",dbFile = dbConn,identifier = "")
    query_desc<- dbGetQuery(dbConn,"select description from Seqs")$description
    query_desc<-gsub("gi\\|.*\\|.*\\|.*\\| ","",query_desc) #simplify tree identifier
    query_desc<-gsub(input$raw_name,input$new_name,query_desc)
    Add2DB(myData = data.frame(identifier=query_desc,stringsAsFactors = FALSE),dbConn)
    cons <- IdConsensus(dbConn,threshold=0.3,minInformation=0.1)
    d <- DistanceMatrix(cons,correction = "Jukes-Cantor")
    dend<-IdClusters(d,method="ML",type = "dendrogram",myXStringSet=cons)
    
    pdf(paste0("www/",input$species,".pdf"),width = 11,height = 8.5)
    for (i in interorf) {
      json_abb <- blast_json$BlastOutput2$report$results$search$hits[[i]] ###abbreviate code
      file_name <- paste0("ORF",i,"_",input$species)
      orf_hit_desc<-c()
      for (j in 1:dim(json_abb)[1]){
        id <- json_abb$description[[j]]$id
        orf_hit_desc_pre <- substring(strsplit(id,":")[[1]][1],gregexpr("_",id)[[1]][1]+1)#simplify query description
        orf_hit_trans <- as.character(virus_id$title[which(orf_hit_desc_pre == virus_id$gb)])
        orf_hit_trans <- gsub(input$raw_name,input$new_name,orf_hit_trans)
        orf_hit_desc <- c(orf_hit_desc,orf_hit_trans)
      }
      
      ####change color of phylo tree
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
      })
      ####plot tree with color change
      dL <- dendrapply(dend, colLab)
      p<- par(mar=c(1,1,1,30),xpd=TRUE)
      plot(dL,yaxt="n",horiz = TRUE,main =file_name)
      par(p)
    }
    dev.off()
    
    tags$iframe(style="height:600px; width:100%", src=paste0(input$species,".pdf"))
  })
  
  output$downloadHTML <- downloadHandler(
    filename = paste0(input$species,".html"),
    content = function(file){
      file.copy(paste0("www/",input$species,".html"),file)
    }
  )
  
  
  output$downloadPDF <- downloadHandler(
    filename = paste0(input$species,".pdf"),
    content = function(file){
      file.copy(paste0("www/",input$species,".pdf"),file)
    }
  )
  
})

