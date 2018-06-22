library(shiny)
library(data.table)
library(jsonlite)
library(shinydashboard)
library(seqinr)
library(DECIPHER)
source("config.R")

server <-shinyServer(function(input, output) {
    
    
    # upload zip files
    output$uploadfiles <- renderTable({
        
        inFile <- input$file
        if (is.null(inFile)) return(NULL)
        files <- unzip(inFile$datapath)
        gsub("./.*/","",files)
        
    })
    
    # the select box for choosing the interest virus, default is AJ489258.1
    output$intervirus <- renderUI({
        inFile <- input$file
        files <- unzip(inFile$datapath)
        selectInput("intervirus", "Interest Virus:",
        choices = gsub(".fasta","",gsub("./.*/","",files)),
        selected = gsub(".fasta","",gsub("./.*/","",files))[7]
        )
    })
    # the select box for choosing interest virus, default is ORF1
    output$interorf <- renderUI({
        blast_json <- fromJSON(paste0("blastp/",input$intervirus,".blast"))
        selectInput("interorf", "Interest ORF:",
        choices = paste0("ORF",1:length(blast_json$BlastOutput2$report$program),"_",input$intervirus),
        selected =paste0("ORF",1:length(blast_json$BlastOutput2$report$program),"_",input$intervirus)[1]
        )
    })

    
    # plot raw phylo tree
    output$plot1 <- renderPlot({
        
        inFile <- input$file
        if (is.null(inFile)) return(NULL)
        files <- unzip(inFile$datapath)
        
        ## merge all virus fasta files into one fasta file named "virus_db.fasta"
        system(paste("cat",paste0(files,collapse = " "),">",virus_db))
        
        ## extract sequence name from "virus_db.fasta" and write into file "virus.id",
        ## which would be used to generate identifier for phylo tree
        system(paste0("grep '>' ",virus_db," | sed 's/>//g' > ",virus_id))
        
        ## mafft is applied to do multiple sequence alignment,which is necessary for plotting phylo tree.
        ## input file is virus_db.fasta,output file is virus_msa.fasta
        system(paste0(mafft_path," ",virus_db," > ",virus_msa_db))
        
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
    
    
    # call ORF via NCBI ORFfinder
    output$orfdb <- renderTable({
        inFile <- input$file
        files <- unzip(inFile$datapath)
        virus_gb <- c()
        orf_number <- c()
        
        ##function "withProgress","incProgress" => show progress
        withProgress("Calculating Progress",value = 0,min=0,max = length(files),{
            for (i in 1:length(files)){
                input$refresh ###refresh button
                incProgress(1, detail = paste("Analyzing virus ", i)) ###show progress
                
                ###get genbank id for each virus
                gb <- gsub(".fasta","",gsub("./.*/","",files[i]))
                ###exec ORFfinder to get orfs for each virus
                orf_cmd <- paste0(orf_path," -id ",gb," -ml ",input$minorflen," -s 0 -out ",paste0("ORF/",gb,".orfs"))
                system(orf_cmd)
                
                ###show orf summary table,which comprising of virus genbank id and orf number of each virus
                
                if (file.info(paste0("ORF/",gb,".orfs"))$size != 0){
                    virus <- read.fasta(paste0("ORF/",gb,".orfs"))
                    virus_gb <- c(gb,virus_gb)
                    orf_number <- c(length(virus),orf_number)
                }
                
            }
        })
        
        ## merge ORF files of all virus into one big ORF file named "orf_db.fasta"
        system(paste0("cat ORF/*.orfs > ",orf_db))
        ## makeblastdb, input file is "orf_db.fasta"
        system(paste0(blastdb_path," -in ",orf_db," -parse_seqids -dbtype prot"))
        
        orfdb <- data.frame(virus_gb,orf_number)
        
    })
    
    
    #blastp
    output$blastp <- renderTable({
        
        #outfmt = json
        system(paste0(blastp_path," -query ORF/",input$intervirus,".orfs -db ",orf_db," -outfmt 15"," -evalue ",input$evalue," > blastp/",input$intervirus,".blast"))
        blast_json <- fromJSON(paste0("blastp/",input$intervirus,".blast"))
        hits_number <- c()
        for (i in 1:length(blast_json$BlastOutput2$report$program)) {
            hits_number<- c(hits_number,dim(blast_json$BlastOutput2$report$results$search$hits[[i]])[1])
        }
        
        #show blast summary table
        blast_result <- data.frame("query_title"=blast_json[["BlastOutput2"]][["report"]][["results"]][["search"]][["query_title"]],
        "query_length"=blast_json[["BlastOutput2"]][["report"]][["results"]][["search"]][["query_len"]],
        "hits_number"= hits_number
        )
        
    })
    
    # plot phylo tree and show multiple sequence alignment
    output$cal <- renderTable({
        ##code from DECIPHER:plot raw phylo tree, same as in ouput$plot1
        virus_id <- fread(virus_id)
        virus_id <- data.frame("gi"= virus_id$V2,"gb"=virus_id$V4,"title"=virus_id$V5)
        blast_json <- fromJSON(paste0("blastp/",input$intervirus,".blast"))
        dbConn <-dbConnect(SQLite(),":memory:")
        Seqs2DB(virus_msa_db,type = "FASTA",dbFile = dbConn,identifier = "")
        query_desc<- dbGetQuery(dbConn,"select description from Seqs")$description
        Add2DB(myData = data.frame(identifier=query_desc,stringsAsFactors = FALSE),dbFile = dbConn)
        query_desc<-gsub("gi\\|.*\\|.*\\|.*\\| ","",query_desc) #simplify tree identifier
        query_desc<-gsub(input$raw_name,input$new_name,query_desc)
        #query_desc <- gsub("Tomato","TYLVCCCCCCCC",query_desc)
        Add2DB(myData = data.frame(identifier=query_desc,stringsAsFactors = FALSE),dbConn)
        cons <- IdConsensus(dbConn,threshold=0.3,minInformation=0.1)
        d <- DistanceMatrix(cons,correction = "Jukes-Cantor")
        dend<-IdClusters(d,method="ML",type = "dendrogram",myXStringSet=cons)
        
        
        
        withProgress("Calculating Progress",value = 0,min=0,max = dim(blast_json$BlastOutput2)[1],{
            
            i_stat<-c()
            for (i in 1:dim(blast_json$BlastOutput2)[1]) {
                json_abb <- blast_json$BlastOutput2$report$results$search$hits[[i]] ###abbreviate code
                incProgress(1, detail = paste("Doing part", i))
                
                ### pdf phylo tree for each orf in interest virus
                
                if (length(json_abb)>1){
                    i_stat <- c(i,i_stat)
                    pdf(paste0("www/ORF",i,"_",input$intervirus,".pdf"),width = 11,height = 8.5)
                    file_name <- paste0("ORF",i,"_",input$intervirus)
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
                        msa_query_title_cmd <- paste0("echo '",query_title,"' >> ",paste0("output/msa/",file_name,".fasta"))
                        msa_query_cmd <- paste0("echo '",query,"' >> ",paste0("output/msa/",file_name,".fasta"))
                        system(msa_query_title_cmd)
                        system(msa_query_cmd)
                        orf_hit_desc_pre <- substring(strsplit(id,":")[[1]][1],gregexpr("_",id)[[1]][1]+1)#simplify query description
                        orf_hit_trans <- as.character(virus_id$title[which(orf_hit_desc_pre == virus_id$gb)])
                        orf_hit_trans <- gsub(input$raw_name,input$new_name,orf_hit_trans)
                        #orf_hit_trans <- gsub("Tomato","TYLVCCCCCCCC",orf_hit_trans)
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
                    }
                    )
                    
                    ####plot tree with color change
                    dL <- dendrapply(dend, colLab)
                    p<- par(mar=c(1,1,1,30),xpd=TRUE)
                    plot(dL,yaxt="n",horiz = TRUE,main =file_name)
                    par(p)
                    dev.off()
                    
                    ####show alignment
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
        
        ##merge alignment results
        system(paste0("cat ",paste0("output/msa/ORF",i_stat,"_",input$intervirus,".html",collapse = " ")," > www/",input$intervirus,".html"))
        
        ##merge phylo tree
        system(paste0("pdftk ",paste0("www/ORF",i_stat,"_",input$intervirus,".pdf",collapse = " ")," cat output"," www/",input$intervirus,".pdf"))
        
    })
    
    
    #show alignment of interest ORF
    output$msa <- renderUI({
        return(includeHTML(paste0("output/msa/",input$interorf,".html")))
    })
    
    #show phylo tree of interest ORF
    output$phylo <- renderUI({
        tags$iframe(style="height:600px; width:100%", src=paste0(input$interorf,".pdf"))
    })
    
    #download merged phylo tree result
    output$downloadPDF <- downloadHandler(
    filename = paste0(input$intervirus,".pdf"),
    content = function(file){
        file.copy(paste0("www/",input$intervirus,".pdf"),file)
    }
    )
    
    #download merged alignment result
    output$downloadHTML <- downloadHandler(
    filename = paste0(input$intervirus,".html"),
    content = function(file){
        file.copy(paste0("www/",input$intervirus,".html"),file)
    }
    )
    
    
    
})
