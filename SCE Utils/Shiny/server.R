# Shiny app for SCE classification
# server.R - Server side functionality
#
# Lines marked with '#!#' are where you should put in the details of your data.
# You need to input the path and name to your mosaiCatcher output in this file.

# Libraries
library(shiny)
library(shinythemes)
library(data.table)
library(ggplot2)
library(htmlwidgets)
library(shinyjs)

source("global.R")

# Define server logic
shinyServer(function(input, output, session) {

  ## Data to classify. This should be a mosaiCatcher output file. CHANGE THIS!
  f_in <- "zcat your_path/count.table.txt.gz"  #!#
  # Name or short description of the dataset to be shown under the plots.
  # Example: "Skin fibroblasts, 500kb bins"
  data_name <<- "" #!#

  ## Storing data. These don't need to be changed, unless you plan to use several different datasets.
  # File suffix to use for the result files the result filename will be in the form of
  #`<user_name>[result_suffix].txt` This suffix is also used to filter out result files for the leaderboard.
  result_suffix <<- "_results"    #!#
  # Where to store the result files. This doesn't need to be changed.
  result_folder <<- "./Results/"  #!#

  ## Variables
  # Other
  var_cellsToShow <<- 5
  var_plotWidth <<- 1200
  var_plotHeight <<- 575
  var_UIPadding <<- 10
  # How close the click has to be (in bp)
  var_pad <<- 2e6


  # Initialize reactive variables
  # Where we are
  r_counter <<- 1
  makeReactiveBinding("r_counter")
  # User
  r_name <-NULL
  makeReactiveBinding("r_name")
  # File to write results into
  r_resultFile <- ""
  makeReactiveBinding("r_resultFile")
  # How many cells were shown to the user in the last plot. This is only relevant when there are too few chromosomes to display.
  r_cellsShown <<- var_cellsToShow
  makeReactiveBinding("r_cellsShown")
  # Next cells to show
  r_nextData <- NULL
  makeReactiveBinding("r_nextData")
  # Flag for tester account
  r_tester <- FALSE
  makeReactiveBinding("r_tester")
  # Leaderboard
  r_results <- NULL
  makeReactiveBinding("r_results")
  #Flag for having gone through it all
  the_end <<- F
  makeReactiveBinding("the_end")


  # Read all counts
  #message("Reading in data...")
  d <- fread(f_in, showProgress = F)
  #message("Please wait, calculating everything needed...")
  d <- initialize_data(d)

  # All breakpoints ..
  potential_SCEs <- d[, .(end=max(end)), by=.(cell_name, cell, sample, chrom, cnsc)]
  # ..where there are at least two breakpoints per cell..
  potential_SCEs <- merge(potential_SCEs[(duplicated(potential_SCEs, by=c("chrom", "cell_name"), fromLast=F))],
                          potential_SCEs[(duplicated(potential_SCEs, by=c("chrom", "cell_name"), fromLast=T))], all = T)
  #.. and not looking at breakpoints that denote the end of the chromosome
  potential_SCEs <- potential_SCEs[potential_SCEs[, .I[end != max(end)], by=chrom]$V1]

  # Choose all chromosomes on all the cells (ignore bins)
  sample_cell_list <- unique(d[order(chrom, cell_name),.(chrom, cell_name, sample)])
  # Take only the ones that have at least one breakpoint
  sample_cell_list <- sample_cell_list[paste(cell_name, chrom, sep="_") %in% paste(potential_SCEs$cell_name, potential_SCEs$chrom, sep="_")]
  # 20 of each chromosome for all the donors #!!
  sample_cell_list <- sample_cell_list[, .SD[1:20], by=list(chrom, sample)] #!!
  # Remove cases where there weren't enough eligible chromosomes, and set counter (id)
  sample_cell_list <- sample_cell_list[!is.na(sample_cell_list$cell_name)][, id := 1:.N]

  # Only looking at breakpoints from the cells selected above
  potential_SCEs <- potential_SCEs[paste(cell_name, chrom, sep="_") %in% paste(sample_cell_list$cell_name, sample_cell_list$chrom, sep="_")]

  # This will store the breakpoints the user selects
  # This is here instead of with the other variables, because it needs 'potential_SCEs'
  r_selectedSCEs <- potential_SCEs[0,]
  makeReactiveBinding("r_selectedSCEs")

  # Workaround to stop drawing the plot twice for each click
  clickSaved <- reactiveValues(singleclick = NULL)
  observeEvent(eventExpr = input$plotClick, handlerExpr = { clickSaved$singleclick <- input$plotClick })

  ## Show login page, and after getting the name of the user, initialize counter, then show main page
  # User initialization
  observeEvent(input$submitButton, {
    r_name <<- input$name
    r_cellsShown <<- var_cellsToShow

    if (r_name != "Just testing") {
      # A result file for each user
      r_resultFile <<- paste0(result_folder, r_name, result_suffix, ".txt")

      # If there are already results on file, carry on from where we left off
      if (file.exists(r_resultFile)) {
        prev_results <<- fread(r_resultFile, header = T, showProgress = F)
        if (r_name %in% prev_results$user) {
          r_counter <<- counter <- max(prev_results[user == r_name]$counter)
        } else {
          r_counter <<- 1
          fwrite(list("cell_name", "cell", "sample", "chrom", "cnsc", "end", "user", "time", "counter"), r_resultFile, append = F, showProgress = F)
        }
        # New user
      } else {
        r_counter <<- 1
        fwrite(list("cell_name", "cell", "sample", "chrom", "cnsc", "end", "user", "time", "counter"), r_resultFile, append = F, showProgress = F)
      }
    } else {
      # User is a tester, create random starting point. Don't write results.
      r_tester <<- TRUE
      r_counter <<- sample.int(nrow(sample_cell_list), 1)
    }

    # Switch from the signup page to the main UI
    hide("signup")
    show("mainContent")
    # Show greeting to user on first page.
    #show("firstPage")
  })

  # Surplus bells and whistles
  # Input greeting to user on first page.
  #output$greeting <- renderText({
  #  input$submitButton
  #  isolate({
  #    greeting <- paste0("Hello ", input$name, ".")
  #  })
  #})

  # On 'Yes'-button click, write selected SCEs into memory, then reset counters
  observeEvent(input$yesButton, {

    # Reconcile actual cells shown
    if (r_cellsShown == var_cellsToShow) {
      r_counter <<- r_counter + var_cellsToShow
    } else {
      r_counter <<- r_counter + r_cellsShown
      r_cellsShown <<- var_cellsToShow
    }

    # Add user info
    r_selectedSCEs$user <- rep(r_name)
    r_selectedSCEs$date <- Sys.time()
    r_selectedSCEs$counter <- r_counter
    # Should add field for dataset! and a corresponding check at user initialization

    # Save results for everyone except testers
    if (!r_tester) {
      fwrite(r_selectedSCEs, r_resultFile, append = T, showProgress = F)
    }

    # The greeting
    #hide("firstPage")

    # Reset variables
    r_selectedSCEs <<- potential_SCEs[0,]
    clickSaved$singleclick <<- NULL
  })

  # Create plot of segment colours and binned W and C reads
  output$plot <- renderPlot({

    # Redraw on every click of the button, but not otherwise
    input$yesButton
    isolate({

      # Number of cells to show:
      if (sample_cell_list[chrom == sample_cell_list$chrom[r_counter], max(id)] == r_counter) {
        last_cell <- r_counter
        r_cellsShown <<- 1
      } else {
        max_id <- sample_cell_list[chrom == sample_cell_list$chrom[r_counter], max(id)]
        if (!is.null(max_id)) {
          last_cell <- min((r_counter + var_cellsToShow - 1), max_id)
          if((r_counter + var_cellsToShow - 1) > max_id) {
            r_cellsShown <<- (last_cell - r_counter)
          }
        } else {
          last_cell <- (r_counter + var_cellsToShow - 1)
          r_cellsShown <- var_cellsToShow
        }
      }
      #print(last_cell)
      # Select cell "r_counter" and a few more
      r_nextData <<- d[sample_cell_list[r_counter:last_cell, ], on = .(chrom, cell_name)]

      # Plotting function in global.R
      if (last_cell == nrow(sample_cell_list)) {
        the_end <<- T
      }
      main_plot(r_nextData)
    })
  }, height = var_plotHeight, width = var_plotWidth)

  # Create plot of breakpoints to overlay on top of the base plot. This speeds things up, as only the top plot needs to be
  # redrawn on breakpoint selection clicks
  output$bkpPlot <- renderPlot({

    # Draw if and only if plot area clicked or yes button pressed
    clickSaved$singleclick
    input$yesButton
    isolate({
      g <- clickSaved$singleclick
      if (!is.null(g)) {

        # Find out where the click was, and if hit any of the breakpoints
        close_lines <- potential_SCEs[0,]
        new_line <- potential_SCEs[0,]
        close_lines <- potential_SCEs[chrom==r_nextData$chrom[1] &
                                        end >= (g$x - var_pad) &
                                        end <= (g$x + var_pad) &
                                        sample == g$panelvar1 &
                                        cell == g$panelvar2]
        if (!is.null(close_lines) && nrow(close_lines)>0) {
          if ( nrow(close_lines) == 1 ) {
            new_line <- close_lines[1]
          } else {
            new_line <- close_lines[which.min(abs(end - g$x)),][1,]
          }
          r_selectedSCEs <<- rbind(r_selectedSCEs, new_line)
          # remove _all_ duplicates if any (toggle mode)
          r_selectedSCEs <<- r_selectedSCEs[!(duplicated(r_selectedSCEs) | duplicated(r_selectedSCEs, fromLast = TRUE)), ]
        }
      }

      # Plotting function in global.R
      bkp_plot(r_nextData, potential_SCEs, r_selectedSCEs)
    })
  }, height = var_plotHeight, width = var_plotWidth, bg="transparent")

  ## Surplus bells and whistles.
  # Info about selected breakpoints in a text field under the plot. Cursor hover info would probably be more useful.
  #  output$bkps_info <- renderText({
  #    clickSaved$singleclick
  #    input$yesButton
  #    isolate({
  #      xy_str <- function(g) {
  #        if(nrow(g)==0)
  #          return(" [no breakpoints selected] \n")
  #        g <- g[order(sample, cell, end)]
  #        info_string <- ""
  #        for(i in 1:nrow(g)) {
  #          info_string <- paste0(info_string, "[", i , "][Sample: ", g$sample[i], ", cell: ", g$cell[i], ", at ", round(g$end[i]/1e6), "MB] \n", collapse = "")
  #        }
  #        return(info_string)
  #      }
  #      paste0("You have selected the following breakpoints: \n\n", xy_str(r_selectedSCEs))
  #    })
  # })

  # Listen for click on the leaderboard tab. Only refresh scores on click.
  observeEvent(input$Tabs, {
    if (input$Tabs == "score_tab") {
      r_results <<- get_results(result_folder, result_suffix)
    }
  })

  # Find out high scores!
  output$scores <- renderText({
    if (!is.null(r_results)) {
      r_results <- r_results[, count_cells:=max(counter), by=user]
      r_results <- r_results[, count_bkps:=.N, by=user]
      r_results <- r_results[order(count_cells, decreasing = T)]
      cells = "\tUser\t\tChromosomes\tBreakpoints\n"
      rank = 1
      for (i in unique(r_results$user)) {
        cells = paste0(cells, rank, ".\t", i, "\t\t", (r_results[user==i]$count_cells[1]-1), "\t\t", r_results[user==i]$count_bkps[1], "\n")
        rank = rank+1
      }
      paste0(cells)
    }
  })

  output$info <- renderText({
    paste0("Data: ", data_name)
  })

  session$onSessionEnded(stopApp)
})
