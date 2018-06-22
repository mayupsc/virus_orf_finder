library(shiny)
library(data.table)
library(jsonlite)
library(shinydashboard)
library(seqinr)
library(DECIPHER)
ui <- dashboardPage(


  dashboardHeader(title = "Virus ORF Alignment"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("WorkFlow", tabName = "makedatabase", icon = icon("th")),
      menuItem("Uploading Files", tabName = "uploadfiles", icon = icon("file",lib = "font-awesome")),
      menuItem("ORFs Database",tabName = "inorfs",icon = icon("database")),
      menuItem("BLAST", tabName = "blast", icon = icon("cogs")),
      menuItem("Calculator",tabName = "cal",icon = icon("calculator")),
      menuItem("Hits Alignment and Distribution", tabName = "HAD", icon = icon("bar-chart-o"))

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
      tabItem(tabName = "makedatabase",
              h3("Work Flow"),
              fluidRow(
                
              img(src="workflow.png",height="70%",width="80%")
              )
      ),
      tabItem(tabName = "uploadfiles",
              h3("Virus database"),
              fluidRow(
               box(width = 4,
                       fileInput("file", "input virus genome files(.zip)"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                    img(src="giphy.gif"))
                ),
              box(width = 8,
                title = "Tree Identifier Simplifier", 
                solidHeader = TRUE,
                textInput("raw_name","full name of virus genus","Tomato yellow leaf curl"),
                textInput("new_name","short name of virus genus:","TYLC")
                  )
             
              ),
              fluidRow(
                box(
                  width = 4, status = "info",
                  title = "Virus Genome Files:",
                  tableOutput("uploadfiles")
                ),
                box( width = 8, title = "Phylogenetic Tree",solidHeader = TRUE,
                     collapsible = TRUE, plotOutput("plot1",height = "880px"))

              )
      ),
      tabItem(tabName = "blast",
              h3("BLAST"),
              fluidRow(
                box(title = "Select Interest Virus", solidHeader = TRUE,uiOutput("intervirus"))
                    ),
               fluidRow(
                 column(8,
                       textInput("evalue","Expectation value (E) threshold for BLAST:","1e-10")

                      ),
                 box(width = 8,collapsible = T,

                     tableOutput("blastp")
                 )

      )
    ),
    tabItem(tabName = "cal",
            h3("Calculate"),
            h4("Waiting..."),
            tableOutput("cal"),
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             img(src="wait.gif")
            )
            ),
    tabItem(tabName = "inorfs",
            h3("Make ORF Database"),
            fluidRow(
              box(
                width = 3,
                textInput("minorflen","Minimal ORF length (nt) for ORF finder:","50"),
                actionButton(class = "btn-sm", "refresh",
                              tagList(icon("refresh"), "Refresh"))
              ),
              box(
                width = 8,
                tableOutput("orfdb")
              )
              )

            ),

    tabItem(tabName = "HAD",
            fluidRow(
              box(width = 4,title = "Select ORF",solidHeader = T,uiOutput("interorf"))
            ),
            
            fluidRow(
              box(downloadButton("downloadHTML","DownloadHTML")),
              box(width = 12, title = "ORF Hits Alignment", collapsible = T,solidHeader = TRUE,
                  htmlOutput("msa")
                  
              )
              
            ),
            
            fluidRow(
              box(downloadButton("downloadPDF","DownloadPDF")),
              box(width = 12,title = "Phylogenetic Tree",solidHeader = T,collapsible = T,
                  uiOutput("phylo"))
              
              
            )
            
            )

)
)
)