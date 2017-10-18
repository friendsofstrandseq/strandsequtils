# Shiny App for Expert SCE Annotation
A web application for semi-manual SCE annotation. 
## Installation
Download the `server.R`, `ui.R` and `global.R` files and place them in the same folder with each other.
## Customization
Edit `server.R` with the details of your data on the lines marked with `'#!#'`. Edit `ui.r` with the names 
of your users on the lines marked with `'#!#'`. More information is provided in the aforementioned files. **NOTE** You need to 
have a mosaiCatcher output file, and you need to edit the location of said file into the `server.R`, for this app to work. 
## Running the App
**From a terminal:** execute the following, where `~/shinyapp` should be replaced with the path to your application:
```
R -e "shiny::runApp('~/shinyapp', launch.browser=TRUE)"
```
**From RStudio:** open `server.R` in your RStudio, and click on the _Run App_ text in the top right of the source code pane.
