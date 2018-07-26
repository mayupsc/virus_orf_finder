library(shiny)
library(data.table)
library(jsonlite)
library(shinydashboard)
library(seqinr)
library(DECIPHER)
ui <- dashboardPage(


  dashboardHeader(title = "Virus ORF Alignment",disable = F),
  dashboardSidebar(
    fileInput("file", "input virus genome files(.zip)"),
  sidebarMenu(
    menuItem("Work Flow", tabName = "workflow", icon = icon("th"))
    ),
  
    textInput("minorflen","Minimal ORF length (nt) for ORF finder:","50"),
    textInput("evalue","Expectation value (E) threshold for BLAST:","1e-10"),
    actionButton(class = "btn-sm", "refresh",tagList(icon("refresh"), "Refresh")),
    
 
   sidebarMenu(
      menuItem("Make Database", tabName = "database", icon = icon("file",lib = "font-awesome"))
   ),
   
    
    div(style = "padding-left: 15px; padding-top: 70px;",
        h4("contact "),
        p(class="small"," "),
        p(class = "small","Xing Fu Ph.D. "),
        p(class = "small","Bioinfomatics Core Facility")

    )
   ),
  dashboardBody(
    
    tabItems(
      tabItem(tabName = "workflow",
              h3("Work Flow"),
              fluidRow(
              img(src="workflow.png",height="70%",width="80%"),
              verbatimTextOutput("preparation")
              )
      ),
      tabItem(tabName = "database",
              h3("ORF database"),
              fluidRow(
                tabBox(width = 12, height = NULL,

                  tabPanel("orfdatabase", 
                                h3("ORF database"),
                                fluidRow(
                                  box(width = 8,
                                      title = "Tree Identifier Simplifier", 
                                      solidHeader = TRUE,
                                      textInput("raw_name","full name of virus genus","Tomato yellow leaf curl"),
                                      textInput("new_name","short name of virus genus:","TYLC")
                                  )
                                ),
                                fluidRow(
                                  box(
                                    width = 5,
                                    title = "ORF database:",
                                    dataTableOutput("orfdatabase")
                                  ),
                                  box( width = 7, title = "Phylogenetic Tree",solidHeader = TRUE,
                                       collapsible = TRUE, plotOutput("plot1",height = "880px"))
                                )                                  
                                 ),
                       
                  tabPanel("blast",              
                        h3("BLAST"),
                        fluidRow(
                                  box( width = 12, title = "ORF Viewer",solidHeader = TRUE,
                                       collapsible = TRUE, plotOutput("plot2",height = "440px")
                                       ),
                                  box(width = 12,collapsible = T,
                                      dataTableOutput("blastp")
                                       )
                                )
                           ),
                  tabPanel("Hits Alignment & ORF Distribution",
                           fluidRow(
                             box(downloadButton("downloadHTML","DownloadHTML")),
                             box(width = 12, title = "ORF Hits Alignment", collapsible = T,solidHeader = TRUE,htmlOutput("msa"))
                           ),
      
                           fluidRow(
                                  box(downloadButton("downloadPDF","DownloadPDF")),
                                  box(width = 12,title = "Phylogenetic Tree",solidHeader = T,collapsible = T,uiOutput("phylo"))
                                    )
                           )
                                
                                
                       )
              )
                )


)
)
)
