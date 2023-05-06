library(shiny)
library(ggplot2)
library(colourpicker)
library(tidyverse)
library(DT)
library(shinythemes)
library(clusterProfiler)
library(ReactomePA)
library(pathview)
source("function.R")

options(shiny.maxRequestSize=100*1024^2)
# Define UI for the samples plot sub-tab
ui_sample_plot <- fluidPage(
  # layout
  sidebarLayout(
    # side bar
    sidebarPanel(
      # select a plot type
      selectInput(inputId = "sample_plot_type",
                  label = "Select a plot you want to see",
                  choices = c("scatter plot", "violin plot","line plot"),
                  selected = "scatter plot"),
      # select a column for x axis
      uiOutput(outputId = "plot_x"),
      # select a column for y axis
      uiOutput(outputId = "plot_y"),
      # select a column to group the samples
      uiOutput(outputId = "group_by"),
      # submit
      submitButton(width = "100%"),
    ),
    # main panel
    mainPanel(
      # plot
      plotOutput(outputId = "meta_sample_plot",height = "400px"),
    )
  )
)

# Define UI for the samples tab
ui_sample <- fluidPage(
  # layout
  sidebarLayout(
    # side bar
    sidebarPanel(
      # input box for file
      fileInput(inputId = "meta_file",
                label = "Load metadata file",
                accept = ".csv,.tsv",
                placeholder = ".csv or .tsv file"),
      # error message
      uiOutput(outputId = "error_message"),
      width = 3
      ),
    # main panel
    mainPanel(
      # sub tabs
      tabsetPanel(
        tabPanel("Table",
                 p("If a column is filled with the same value,
                 it will not be shown in the table. You can still view it in summary tab.",
                   br(),"Reanalyzed.by, SRA, and BioSample are not shown in the table."),
                 dataTableOutput(outputId = "meta_table")),
        tabPanel("Summary",
                 h3("Summary of common values for all samples"),
                 verbatimTextOutput(outputId = "summary")),
        tabPanel("Statistics Summary",
                 h3("Summary of each column in the trimmed metadata file"),
                 tableOutput(outputId = "meta_table_summary")),
        tabPanel("Plot", ui_sample_plot)
      ),
      width = 9
    )
  )
)

# Define UI for the PCA plot sub-tab
ui_pca_plot <- fluidPage(
  # layout
  sidebarLayout(
    # side bar
    sidebarPanel(
      # select a pca for x an y axis
      uiOutput(outputId = "pca_xy"),
      # submit
      submitButton(width = "100%"),
    width = 3),
    # main panel
    mainPanel(
      # plot
      plotOutput(outputId = "pca_plot",height = "400px"),
    )
  )
)

# Define UI for the counts tab
ui_counts <- fluidPage(
  # layout
  sidebarLayout(
    # side bar
    sidebarPanel(
      # input box for file
      fileInput(inputId = "count_file",
                label = "Load count matrix file",
                accept = ".csv,.tsv",
                placeholder = ".csv or .tsv file"),
      # error message
      uiOutput(outputId = "count_error_message"),
      # filter genes by variance
      sliderInput(inputId = "count_filter_var",
                  label = "keep genes with X percentile of variance",
                  min = 0,
                  max = 100,
                  value = 100,
                  step = 1,
                  ticks = TRUE,
                  width = "100%"),
      # filter genes by non-zero counts
      sliderInput(inputId = "count_filter_nonzero",
                  label = "Keep genes with at least X percetage of non-zero counts",
                  min = 0,
                  max = 100,
                  value = 50,
                  step = 1,
                  ticks = TRUE,
                  width = "100%"),
      # submit
      submitButton(width = "100%"),
      width = 3
    ),
    # main panel
    mainPanel(
      # sub tabs
      tabsetPanel(
        tabPanel("Filter Summary",
                 h3("Summary of data after filtering"),
                 tableOutput(outputId = "count_summary")),
        tabPanel("Filter Diagnostic Plot",
                 plotOutput(outputId = "count_filter_plot")),
        tabPanel("Heatmap",
                 plotOutput(outputId = "count_heatmap")),
        tabPanel("PCA Plot",
                 ui_pca_plot),
      ),
      width = 9
    )
  )
)

# Define UI for the DE tab
ui_de <- fluidPage(
  # layout
  sidebarLayout(
    # side bar
    sidebarPanel(
      # input box for file
      fileInput(inputId = "de_file",
                label = "Load DE analysis results file",
                accept = ".csv,.tsv",
                placeholder = ".csv or .tsv file"),
      # error message
      uiOutput(outputId = "de_error_message"),
      # select a column for x axis
      uiOutput(outputId = "de_input"),

      # submit
      submitButton(width = "100%"),
      width = 3
    ),
    # main panel
    mainPanel(
      # sub tabs
      tabsetPanel(
        tabPanel("Table",
                 # add a newline
                    br(),
                 dataTableOutput(outputId = "de_table")),
        tabPanel("Volcano Plot",
                 plotOutput(outputId = "de_volcano_plot")),
      ),
      width = 9
    )
  )
)

# Define UI for the gene enrichment plot sub-tab
ui_ge_plot <- fluidPage(
  # layout
  sidebarLayout(
    # side bar
    sidebarPanel(
      # select a plot type for enrichment
      radioButtons(inputId = "ge_plot_type",
                   label = "Select a plot type",
                   choices = c("Bar Plot","Enrichment Map"),
                   selected = "Bar Plot"),
      # submit
      submitButton(width = "100%"),
    width = 3),
    # main panel
    mainPanel(
      # plot
      plotOutput(outputId = "ge_plot",height = "600px"),
    )
  )
)

# Define UI for the gene enrichment tab
ui_ge <- fluidPage(
  # layout
  sidebarLayout(
    # side bar
    sidebarPanel(
      # select a gene list
      uiOutput(outputId = "enrichment_input"),
      # submit
      submitButton(width = "100%"),
    width = 3,),
    # main panel
    mainPanel(
      # sub tabs
      tabsetPanel(
        tabPanel("Table",
                 # add a newline
                 br(),
                 dataTableOutput(outputId = "enrichment_table")),
        tabPanel("KEGG Pathway Plot",
                 ui_ge_plot),
      ),
    )
    )
    )



# Define main UI for application
ui <- fluidPage(
  # theme
  theme = shinytheme("cerulean"),
  # title
  titlePanel("BF591 Final Project - Zedong"),
  # tabs layout
    tabsetPanel(
        # samples tab
        tabPanel("Samples",
                 h2("Please load metadata file of your samples"),
                    ui_sample
                    ),
        # counts tab
        tabPanel("Counts",
                 h2("Please load count matrix file of your samples"),
                    ui_counts
                  ),
        # differential expression tab
        tabPanel("Differential Expression",
                h2("Please load  differential expression analysis results file"),
                   ui_de
                ),

        # gene enrichment tab
        tabPanel("Gene Enrichment",
                 h3("Gene list for enrichment analysis is generated from DE results"),
                 p("If you have not loaded DE results, please go to DE tab and load it first."),
                 p("If you want to change the gene list, please go to DE tab and reload DE results."),
                ui_ge
                ),
        # about tab
        tabPanel("About",
                 h2("About this app"),
                 p("This app is designed to help you visualize and explore
                 differential expression results from RNA-seq data."))
                )

)




# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # reactive values
  rv <- reactiveValues(meta_df = NULL,
                       same_list = NULL,
                       count_tibble = NULL,
                       filtered_count_tibble = NULL,
                       de_tibble = NULL)

  # check and load metadata file
  observeEvent(input$meta_file, {
    file_type <- tools::file_ext(input$meta_file$name)
    if (file_type == "csv") {
      df <- read.csv(input$meta_file$datapath, header = TRUE, row.names = 1)
      df <- filter_same_cols(df)
      rv$meta_df <- df[[1]]
      rv$same_list <- df[[2]]
      # clear error message
      output$error_message <- NULL
    }
    else if (file_type == "tsv") {
      df <- read.delim(input$meta_file$datapath, header = TRUE, row.names = 1)
      df <- filter_same_cols(df)
      rv$meta_df <- df[[1]]
      rv$same_list <- df[[2]]
      # clear error message
      output$error_message <- NULL
    }
    else {
      # show error message
      output$error_message <- renderUI({
        tags$div(
          class = "alert alert-danger",
          "Please upload a .csv or .tsv file"
        )
      })
      rv$meta_df <- NULL
      rv$same_list <- NULL
    }
  })

  # check and load count matrix file
  observeEvent(input$count_file, {
    file_type <- tools::file_ext(input$count_file$name)
    if (file_type == "csv") {
      df <- read.csv(input$count_file$datapath, header = TRUE, row.names = 1)
      rv$count_tibble <- make_count_tibble(df)
      # clear error message
      output$count_error_message <- NULL
    }
    else if (file_type == "tsv") {
      df <- read.delim(input$count_file$datapath, header = TRUE, row.names = 1)
      rv$count_tibble <- make_count_tibble(df)
      # clear error message
      output$count_error_message <- NULL
    }
    else {
      # show error message
      output$count_error_message <- renderUI({
        tags$div(
          class = "alert alert-danger",
          "Please upload a .csv or .tsv file"
        )
      })
      rv$count_tibble <- NULL
    }
  })

  # check and load differential expression file
  observeEvent(input$de_file, {
        file_type <- tools::file_ext(input$de_file$name)
        if (file_type == "csv") {
        df <- read.csv(input$de_file$datapath, header = TRUE, row.names = 1)
        rv$de_tibble <- make_de_tibble(df)
        # clear error message
        output$de_error_message <- NULL
        }
        else if (file_type == "tsv") {
        df <- read.delim(input$de_file$datapath, header = TRUE, row.names = 1)
        rv$de_tibble <- make_de_tibble(df)
        # clear error message
        output$de_error_message <- NULL
        }
        else {
        # show error message
        output$de_error_message <- renderUI({
            tags$div(
            class = "alert alert-danger",
            "Please upload a .csv or .tsv file"
            )
        })
        rv$de_tibble <- NULL
        }
    })


# Server function for samples tab
  # render group_by input
  output$group_by <- renderUI({
    # only show group_by input when meta_df is not NULL
    if (!is.null(rv$meta_df)) {
      selectInput(inputId = "group_by_output",
                  label = "group the samples by (e.g. condition, time point, treatment)",
                  choices = c("All in one Group", colnames(rv$meta_df)),
                  selected = NULL)
    }
  })

  # render plot_x input
  output$plot_x <- renderUI({
    # only show plot_x input when meta_df is not NULL
    if (!is.null(rv$meta_df)) {
      selectInput(inputId = "plot_x_output",
                  label = "Select a column for x axis",
                  choices = colnames(rv$meta_df),
                  selected = NULL)
    }
  })

  # render plot_y input
  output$plot_y <- renderUI({
      # only show plot_y input when meta_df is not NULL
      if (!is.null(rv$meta_df)) {
      selectInput(inputId = "plot_y_output",
                  label = "Select a column for y axis",
                  choices = colnames(rv$meta_df),
                  selected = NULL)
      }
  })

  # render meta_table
  output$meta_table <- renderDataTable({
    if (!is.null(rv$meta_df)) {
      datatable(
      rv$meta_df,
      options = list(
        fixedHeader = TRUE,
        scrollX = TRUE
      ),
      rownames = FALSE
      )}
  })

  # render summary
  output$summary <- renderPrint({
      if (!is.null(rv$meta_df)) {
        meta_common_summary(rv$same_list)
      }
    })

  # render meta_table_summary
  output$meta_table_summary <- renderTable({
    if (!is.null(rv$meta_df)) {
      meta_table_summary(rv$meta_df)
    }
    })

  # render sample plot
  output$meta_sample_plot <- renderPlot({
    req(input$plot_x_output, input$plot_y_output, input$group_by_output, input$sample_plot_type, rv$meta_df)
    if (input$group_by_output== "All in one Group") {
      p <- plot_sample(rv$meta_df,
                  input$plot_x_output,
                  input$plot_y_output,
                  input$sample_plot_type )
      print("1")
    }
    else {
      p <- plot_sample(rv$meta_df,
                  input$plot_x_output,
                  input$plot_y_output,
                  input$sample_plot_type,
                  input$group_by_output)
      print("2")
    }
    return(p)
  },height = 400, width = 600)

# Server function for counts tab
  # render filter summary
  output$count_summary <- renderTable({
    req(rv$count_tibble,input$count_filter_nonzero, input$count_filter_var)
    count_tibble <- anno_filter(rv$count_tibble,input$count_filter_var,input$count_filter_nonzero)
    n_sample <- ncol(rv$count_tibble) - 5
    n_gene <- nrow(rv$count_tibble)
    n_pass <- sum(count_tibble$keep)
    n_fail <- n_gene - n_pass
    ret <- data.frame("Number of samples" = n_sample,
               "Number of genes" = n_gene,
               "Number of genes passing current filter" = n_pass,
               "Percent of genes passing current filter" = paste0(round(n_pass/n_gene*100, 2), "%"),
               "Number of genes not passing current filter" = n_fail,
               "Percent of genes not passing current filter" = paste0(round(n_fail/n_gene*100, 2), "%"))
    df_long <- ret %>% mutate_all(as.character)%>%
               pivot_longer(cols = everything(),
                            names_to = "Feature",
                            values_to = "Value")
    rv$filtered_count_tibble <- count_tibble
    return(df_long)
    })

  # render filter plot
  output$count_filter_plot <- renderPlot({
    req(rv$count_tibble,input$count_filter_nonzero, input$count_filter_var)
    rv$filtered_count_tibble <- anno_filter(rv$count_tibble,input$count_filter_var,input$count_filter_nonzero)
    p <- plot_filtered_counts(rv$filtered_count_tibble)
    return(p)
    },height = 450, width = 1000)

  # render count heatmap
  output$count_heatmap <- renderPlot({
  req(rv$count_tibble,input$count_filter_nonzero, input$count_filter_var)
  rv$filtered_count_tibble <- anno_filter(rv$count_tibble,input$count_filter_var,input$count_filter_nonzero)
  p <- plot_count_heatmap(rv$filtered_count_tibble)
  return(p)
  },height = 450, width = 900)

  # render pca option
  output$pca_xy <- renderUI({
    # only show pca_x and y input when count_tibble and fiters are not NULL
    req(rv$filtered_count_tibble, input$count_filter_nonzero, input$count_filter_var)
    rv$filtered_count_tibble <- anno_filter(rv$count_tibble,input$count_filter_var,input$count_filter_nonzero)

    df <- calc_pca(rv$filtered_count_tibble)
    rv$pca_df <- df
    list(
    selectInput(inputId = "pca_x_output",
                label = "Select a column for x axis",
                choices = colnames(df$importance),
                selected = NULL),
    selectInput(inputId = "pca_y_output",
                label = "Select a column for y axis",
                choices = colnames(df$importance),
                selected = NULL))
    })

  # render pca plot
  output$pca_plot <- renderPlot({
    req(rv$pca_df, input$pca_x_output, input$pca_y_output)
    p <- plot_pca(rv$pca_df, input$pca_x_output, input$pca_y_output)
    return(p)
    })

  # render de option
  output$de_input <- renderUI({
    if (!is.null(rv$de_tibble)) {
    option_list <- colnames(rv$de_tibble)
    ret<-list(
      selectInput(inputId = "de_x_output",
                  label = "Select a column for x axis",
                  choices = option_list,
                  selected = NULL),
      selectInput(inputId = "de_y_output",
                  label = "Select a column for y axis",
                  choices = option_list,
                  selected = NULL),
      selectInput(inputId = "de_p_output",
                  label = "Select the column for p-adjusted values",
                  choices = option_list,
                  selected = NULL),
      # base point color
      colourInput(inputId = "de_base",
                  label = "Base point color",
                  value = "#22577A",
                  closeOnClick = TRUE),
      # Highlight color
      colourInput(inputId = "de_highlight",
                  label = "Highlight point color",
                  value = "#FFCF56",
                  closeOnClick = TRUE),
      # p value
      sliderInput(inputId = "p_slider",
                  label = "Select the magnitude of the p adjusted coloring:",
                  min = -300,
                  max = 0,
                  step = 1,
                  value = -100)
    )
    }
    else {
      ret <- NULL
    }
    return(ret)
    })

  # render de table
  output$de_table <- renderDataTable({
      req(rv$de_tibble)
      if (!is.null(rv$de_tibble)&!is.null(input$de_p_output)&!is.null(input$p_slider)) {
        print(input$de_p_output)
        print(input$p_slider)
        ret <- dplyr::filter(rv$de_tibble,!!sym(input$de_p_output) <= 10^input$p_slider)
        ret <- format_number(ret)

      }
      else {
        ret <- format_number(rv$de_tibble)
      }
      return(ret)
    })

  # render de volcano plot
  output$de_volcano_plot <- renderPlot({
    req(rv$de_tibble, input$de_x_output, input$de_y_output, input$de_p_output, input$de_base, input$de_highlight, input$p_slider)
    p <- volcano_plot(rv$de_tibble,
                      input$de_x_output,
                      input$de_y_output,
                      input$de_p_output,
                      input$p_slider,
                      input$de_base,
                      input$de_highlight
                      )
    return(p)
    },height = 600, width = 800)

  # render enrichment option
  output$enrichment_input <- renderUI({
    if(!is.null(rv$de_tibble)) {
      option_list <- colnames(rv$de_tibble)
      ret<-list(
        selectInput(inputId = "enrichment_rank",
                    label = "Select a column to rank the genes from greatest to least",
                    choices = option_list,
                    selected = NULL),
        selectInput(inputId = "enrichment_genelist",
                    label = "Select the column for Ensemble IDs",
                    choices = option_list,
                    selected = NULL),
        selectInput(inputId = "enrichment_p_output",
                    label = "Select the column for p-adjusted values",
                    choices = option_list,
                    selected = NULL),

        # p value
        sliderInput(inputId = "enrichment_p_slider",
                    label = "Select the magnitude of the p-adjusted threshold:",
                    min = -300,
                    max = 0,
                    step = 1,
                    value = -5)
      )
    }
    else {
      ret <- NULL
    }
    return(ret)
  })

  # render enrichment table
    output$enrichment_table <- renderDataTable({
        req(rv$de_tibble)
        if (!is.null(rv$de_tibble)&
          !is.null(input$enrichment_p_output)&
          !is.null(input$enrichment_p_slider)&
          !is.null(input$enrichment_rank)&
          !is.null(input$enrichment_genelist)) {
        ret <- dplyr::filter(rv$de_tibble,!!sym(input$enrichment_p_output) <= 10^input$enrichment_p_slider)
        rv$ge_tibble <- make_gene_list(ret, input$enrichment_genelist, input$enrichment_rank)
        ret <- format_number(ret)
        }
        else {
        ret <- format_number(rv$de_tibble)
        }
        return(ret)
    })

  # render enrichment plot
  output$ge_plot <- renderPlot({
    req(rv$ge_tibble, input$ge_plot_type)
    if (input$ge_plot_type == "Bar Plot") {
      p <- kegg_bar_plot(rv$ge_tibble)
      return(p)
    }
    else if (input$ge_plot_type == "Enrichment Map") {
      p <- kegg_map_plot(rv$ge_tibble)
      return(p)
    }
  })
}




# Run the application
shinyApp(ui = ui, server = server)