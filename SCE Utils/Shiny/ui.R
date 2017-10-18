# Shiny app for SCE classification
# ui.R User interface
#
# Please put in the names of your users on the lines marked with  '#!#'.

library(shiny)
library(shinythemes)
library(data.table)
library(ggplot2)
library(shinyjs)

source("global.R")

# Define UI
shinyUI(
  fluidPage(

    # To enable plotting two plots on top of each other
    tags$head(tags$style(".fixPosition{position: absolute;}")),

    # Cool buttons
    theme = shinytheme("spacelab"),

    #Showing and hiding stuff
    useShinyjs(),

    # First page
    div(
      id = "signup",
      titlePanel("SCE Classification App"),
      tags$p("Who are you?"),
      tags$p(
        selectInput("name", "Name:",
                    c("Just testing" = "Just testing", # You can leave this in. The results won't be stored.
                      "Ashley" = "Ashley",   #!# Input the names of your users here.
                      "Karen" = "Karen",     #!#
                      "Sascha" = "Sascha",   #!#
                      "Venla" = "Venla")),   #!#
        actionButton("submitButton", "Submit")
      ),
      br(),
      tags$p(tags$b("Note:"), " Please choose your own name, as it is used to determine which chromosomes to show you. You can test the app by choosing ", tags$i("Just testing"), "as your username. Testers will be started on a random chromosome, and the SCEs they pick are not permanently stored."),
      tags$p("Calculations are performed first, so loading the first plot might take some time. Thank you for your patience.")
    ),

    # Main page, hidden at first
    hidden(
      div(
        useShinyjs(),
        id = "mainContent",

        # Application title
        titlePanel("SCE Classification"),

        mainPanel(
          tabsetPanel(
            type = "tabs",
            id = "Tabs",
            tabPanel(
              "Plot",
              # Superfluous bells and whistles
              #hidden(
              #  div(
              #    id="firstPage",
              #    h3(textOutput("greeting"))
              #  )
              #),
              #tags$p("Please click on all the breakpoints (blue dashed lines) that you think border
              #a sister chromatid exchange. Once you've selected all the breakpoints, click on the button below the plot."),
              #hr(),
              fluidRow(
                column(
                  width = 11,
                  offset = 0,
                  div(
                    plotOutput("plot", height = 575),
                    class = "fixPosition"
                  ),
                  div(
                    plotOutput("bkpPlot", height = 575, click = "plotClick"),
                    class = "fixPosition"

                  ),
                  style = paste0("width: 100% ; height: ", (575 + 10),"px")
                ),
                fluidRow(
                  column(11,
                    verbatimTextOutput("info")
                    #    tags$p(verbatimTextOutput("bkps_info"))
                  ),
                  column(1,
                    actionButton("yesButton", "Breakpoints chosen. Next!"),
                    style = "display: inline-block; vertical-align: top"
                  )
                )
              )
            ),
            tabPanel(
              title = "Leaderboard",
              value = "score_tab",
              verbatimTextOutput("scores")
            )
          )
        )
      )
    ),
    padding = 10
  )
)
