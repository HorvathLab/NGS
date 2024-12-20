#' Generate Report
#'
#' This function generates interactive 3D dimensionality reduction plots, histograms, Cell Types Plot, and CopyKat Plot, with the option to save plots individually.
#'
#' @importFrom
#'
#' @param plot_object
#' @param snv_plot_object
#' @param hide_ind_plots 
#' @details
#' The function generates various visualizations for SNV data including:
#' - **3D Plots**: Visualizes metrics such as SNV.N, SNVCount, and TotalVAF.
#' - **Histograms**: Distribution of metrics including TotalVAF, MeanSNVsVAF, and N_VARreadCounts.
#' - **Cell Types Plot**: Shows cell types (e.g., custom classifications) with optional slingshot trajectories.
#' - **CopyKat Plot**: Depicts Copy Number Variations (CNVs) using CopyKat analysis.
#'
#' The output_dir parameter specifies where the plots will be saved if save_each_plot is enabled.
#'
#' @examples
#' # Example usage:
#'
#'
#' @export
generate_report <- function(plot_object, snv_plot_object = NULL, hide_ind_plots = TRUE){

    if(missing(plot_object)){
        stop("Use output object of plot_snv_data() function as first parameter")
    }

        #combine plots.
        #initalizing
    all_plots <- list()
    individual_SNV_html <- NULL
    

    for (plot_name in names(plot_object)) {
    all_plots <- append(all_plots, list(plot_object[[plot_name]]))
    names(all_plots)[length(all_plots)] <- plot_name
    }
   #if(!is.null(snv_plot_object)){
    #for (plot_name in names(snv_plot_object)) {
      #all_plots <- append(all_plots, list(snv_plot_object[[plot_name]]))
      #names(all_plots)[length(all_plots)] <- plot_name
      #}
    #}

    histograms <- plot_object$histograms 
    #plot descriptions
    plot_descriptions <- list(
        "N_sceSNVs" = "Absolute number of sceSNVs per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
        "N_VARreads" = "Sum of the absolute number of reads bearing the variant nucleotide across all sceSNVs in the submitted list per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
        "N_REFreads" = "Sum of the absolute number of reads bearing the reference nucleotide across all sceSNVs in the submitted list per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
        "Total_VAF_RNA" = "Total VAF_RNA is calculated by dividing the sum of N_VAR counts across the sceSNVs (loci) by the total reads covering the sceSNV loci (N_VAR + N_REF), and is intended to provide an assessment of SNV expression magnitude across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
        "Mean_VAF_RNA" = "Mean VAF_RNA is calculated as the mean of the VAF_RNA across all individual sceSNVs and is intended to provide an assessment of the variability of the SNV expression levels across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
        "Median_VAF_RNA" = "Median VAF_RNA is calculated as the median of the VAF_RNA across all individual sceSNVs and is intended to provide an assessment of the central tendency of the SNV expression levels across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
        "Cell types (scType)" = "Classification of individual cells into distinct types using the scType tool. scType employs predefined marker gene sets and a scoring algorithm to analyze single-cell RNA sequencing (scRNA-seq) data, assigning each cell to a specific type based on its gene expression profile. The expected tissue type for the analysis is expected to be submitted by the user.",
        "CNVs (CopyKat)" = "Copy number alterations (CNAs) detected using the CopyKat tool designed to infer genomic CNVs across individual cells using scRNA-seq data. By comparing gene expression profiles, CopyKat identifies regions of the genome that have been amplified or deleted, which can be indicative of genetic abnormalities such as those found in cancer cells.",
        "VAF_RNA" = "Expressed Variant Allele Fraction. For each individual sceSNV VAF_RNA is calculated as the ratio of the number of variant reads (N_VAR) divided by the total number of reads (N_VAR + N_REF) covering the sceSNV locus (VAF_RNA = N_VAR / (N_VAR + N_REF).",
        "N_VAR" = "Absolute number of reads bearing the variant nucleotide at the sceSNV locus.",
        "N_REF" = "Absolute number of reads bearing the reference nucleotide at the sceSNV locus."
        )
  if(!is.null(snv_plot_object)){
    snv_options <- snv_plot_object$snv_options
    plots_json <- snv_plot_object$plots_json
    individual_SNV_html <- paste0('
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <title>Individual sceSNV plot selection</title>
            <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <style>
            body {text-align: center; margin-bottom: 50px;}
            .center-container {display: inline-block; text-align: left;}
            .plot-container {margin-top: 20px; display: grid; grid-template-columns: repeat(2, 1fr); gap: 10px; padding: 0 20px;}
            .plot-item {width: 800px; height: 500px; margin-bottom: 50px;}
            .plot-title {
                cursor: pointer;
                color: blue;
                font-weight: bold;
                text-decoration: none;
                display: inline-block;
            }
            </style>
        </head>
        <body>
            <div class="center-container">
            <h2>Individual sceSNV plot selection</h2>
            <label for="snv">Select SNV:</label>
            <select id="snv">
                ', paste0('<option value="', snv_options, '">', snv_options, '</option>', collapse = ''), '
            </select>
            <br><br>
            <button id="updatePlot">Update Plot</button>
            <p id="selectedSNV"></p>
            </div>
            <br><br>
            <div class="plot-container">
            <div id="plotContainer_VAF" class="plot-item">
                <div class="plot-title" data-plot-id="VAF_RNA">VAF</div>
            </div>
            <div id="plotContainer_N_VAR" class="plot-item">
                <div class="plot-title" data-plot-id="N_VAR">N_VAR</div>
            </div>
            <div id="plotContainer_N_REF" class="plot-item">
                <div class="plot-title" data-plot-id="N_REF">N_REF</div>
            </div>
            </div>
            <script>
            const plots = {
                ', paste0(sapply(plots_json, function(p) {
                paste0('"', p$id, '": ', p$json)
                }), collapse = ", "), '
            };

            function updatePlot() {
                var snv = $("#snv").val();
                var plotIds = ["VAF", "N_VAR", "N_REF"];
                plotIds.forEach(function(plotType) {
                var plotId = "plot_" + plotType + "_" + snv.replace(/:/g, "_");
                var plotData = plots[plotId];
                if (plotData) {
                    Plotly.newPlot("plotContainer_" + plotType, plotData.data, plotData.layout, {
                    responsive: true,
                    autosize: true,
                    width: document.getElementById("plotContainer_" + plotType).offsetWidth,
                    height: document.getElementById("plotContainer_" + plotType).offsetHeight
                    });
                }
                });
                $("#selectedSNV").text("Selected SNV: " + snv);
            }

            $(document).ready(function(){
                $("#updatePlot").click(updatePlot);
                $("#snv").val($("#snv option:first").val());
                updatePlot();
            });
            </script>
        </body>
        </html>'
      )
  }
    
###
    link_html <- '<br><button id="toggleButton" onclick="toggleHistograms()" style="font-size: 20px; padding: 10px 20px; margin: 40px 0;">View Histograms</button><br>'

    toggle_histogram_script <- '
        <script>
            // Ensure histograms are hidden on page load and set button text correctly
            document.addEventListener("DOMContentLoaded", function() {
                var histograms = document.getElementById("histograms");
                var toggleButton = document.getElementById("toggleButton");

                // Initially hide histograms and set button text
                if (histograms.style.display === "" || histograms.style.display === "none") {
                    histograms.style.display = "none";
                    toggleButton.textContent = "View Histograms";
                } else {
                    histograms.style.display = "grid";
                    toggleButton.textContent = "Hide Histograms";
                }
            });

            function toggleHistograms() {
                var histograms = document.getElementById("histograms");
                var toggleButton = document.getElementById("toggleButton");

                // Toggle visibility and button text
                if (histograms.style.display === "none") {
                    histograms.style.display = "grid";
                    toggleButton.textContent = "Hide Histograms";
                } else {
                    histograms.style.display = "none";
                    toggleButton.textContent = "View Histograms";
                }
            }
        </script>
    '
    # JavaScript code for plot title pop-ups
    title_popup_script <- '
        <script>
        document.addEventListener("DOMContentLoaded", function() {
            document.querySelectorAll(".plot-title").forEach(function(element) {
            element.addEventListener("click", function() {
                var plotId = element.getAttribute("data-plot-id");
                var descriptions = {
                "N_sceSNVs": "Absolute number of sceSNVs per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
                "N_VARreads": "Sum of the absolute number of reads bearing the variant nucleotide across all sceSNVs in the submitted list per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
                "N_REFreads": "Sum of the absolute number of reads bearing the reference nucleotide across all sceSNVs in the submitted list per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
                "Total_VAF_RNA": "Total VAF_RNA is calculated by dividing the sum of N_VAR counts across the sceSNVs (loci) by the total reads covering the sceSNV loci (N_VAR + N_REF), and is intended to provide an assessment of SNV expression magnitude across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
                "Mean_VAF_RNA": "Mean VAF_RNA is calculated as the mean of the VAF_RNA across all individual sceSNVs and is intended to provide an assessment of the variability of the SNV expression levels across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
                "Median_VAF_RNA": "Median VAF_RNA is calculated as the median of the VAF_RNA across all individual sceSNVs and is intended to provide an assessment of the central tendency of the SNV expression levels across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
                "Cell types (scType)": "Classification of individual cells into distinct types using the scType tool. scType employs predefined marker gene sets and a scoring algorithm to analyze single-cell RNA sequencing (scRNA-seq) data, assigning each cell to a specific type based on its gene expression profile. The expected tissue type for the analysis is expected to be submitted by the user.",
                "CNVs (CopyKat)": "Copy number alterations (CNAs) detected using the CopyKat tool designed to infer genomic CNVs across individual cells using scRNA-seq data. By comparing gene expression profiles, CopyKat identifies regions of the genome that have been amplified or deleted, which can be indicative of genetic abnormalities such as those found in cancer cells.",
                "VAF_RNA": "Expressed Variant Allele Fraction. For each individual sceSNV VAF_RNA is calculated as the ratio of the number of variant reads (N_VAR) divided by the total number of reads (N_VAR + N_REF) covering the sceSNV locus (VAF_RNA = N_VAR / (N_VAR + N_REF).",
                "N_VAR": "Absolute number of reads bearing the variant nucleotide at the sceSNV locus.",  
                "N_REF": "Absolute number of reads bearing the reference nucleotide at the sceSNV locus."
                };

                var existingTooltip = document.querySelector(`.tooltip-box[data-plot-id="${plotId}"]`);
                if (existingTooltip) {
                document.body.removeChild(existingTooltip);
                return;
                }

                var tooltip = document.createElement("div");
                tooltip.className = "tooltip-box";
                tooltip.setAttribute("data-plot-id", plotId);
                tooltip.innerHTML = descriptions[plotId];

                document.body.appendChild(tooltip);

                var rect = element.getBoundingClientRect();
                var tooltipWidth = tooltip.offsetWidth;
                var tooltipHeight = tooltip.offsetHeight;
                var pageWidth = window.innerWidth;
                var pageHeight = window.innerHeight;

                // Adjust tooltip positioning to ensure it remains within page boundaries
                if (rect.right + tooltipWidth > pageWidth) {
                tooltip.style.left = (rect.left - tooltipWidth - 10) + "px";
                } else {
                tooltip.style.left = (rect.right + 10) + "px";
                }

                if (rect.bottom + tooltipHeight > pageHeight) {
                tooltip.style.top = (rect.top - tooltipHeight - 10) + "px";
                } else {
                tooltip.style.top = (rect.bottom + window.scrollY) + "px";
                }

                window.addEventListener("scroll", function() {
                var updatedRect = element.getBoundingClientRect();
                tooltip.style.left = (updatedRect.left + window.scrollX) + "px";
                tooltip.style.top = (updatedRect.bottom + window.scrollY) + "px";
                }, { passive: true });
            });
            });
        });
        </script>
    '

    # CSS for tooltip and plot titles
    tooltip_css <- '
        <style>
        .tooltip-box {
            position: absolute;
            background-color: white;
            border: 1px solid black;
            padding: 10px;
            z-index: 1000;
            max-width: 300px;
            box-shadow: 0px 0px 10px rgba(0,0,0,0.5);
            font-size: 16px; /* Adjust text size */
        }
        .plot-title {
            cursor: pointer;
            color: blue;
            font-weight: bold;
            text-decoration: none;
            display: inline-block;
        }
        </style>
        '   
    plot_divs <- lapply(seq_along(all_plots), function(i) {
    if(names(all_plots)[i] != "histograms"){
    plot <- all_plots[[i]]
    plot_id <- names(plot_descriptions)[i]
    div(class = "grid-item plot3d", 
        div(class = "plot-title", `data-plot-id` = plot_id, plot_id),  # plot_id for title
        as_widget(plot)
        )}})

    # CSS layout style
grid_html <- tags$html(
  tags$head(
    tags$style(HTML("
      .grid-container {
        display: grid;
        grid-template-columns: repeat(2, 1fr);  /* 1st arg. specifies no of columns */
        gap: 20px 5px;   /* adjust the gap: 1st value for row spacing, 2nd for column spacing */
        justify-content: center; 
      }
      .grid-item {
        padding: 0; 
      }
      .grid-item.plot3d {
        width: 100%;   /* width for snv plots */
        height: 100%;  /* height for snv plots */
        text-align: center; /* center align titles */
      }
      .grid-item.histogram {
        width: 100%;   /* width for histograms */
        height: 100%;  /* height for histograms */
      }
      #histograms {
        display: none;
        grid-template-columns: repeat(2, 1fr); /* 2-column grid for histograms */
        gap: 10px;
        justify-content: center;
      }
      .button-container {
        width: 100%;         
        display: flex;    /* flexbox is for centering */
        justify-content: center;  
        margin: 20px 0;     
      }
      #toggleButton {
        display: inline-block; /* inline-block for better control */
        font-size: 20px;
        padding: 10px 20px;
      }
      body {
        font-size: 20px; 
      }
    ")),
    HTML(tooltip_css)
  ),
  tags$body(
    div(class = "grid-container", plot_divs),
    div(class = "button-container",
      HTML(link_html)
    ),
    div(id = "histograms", 
      lapply(histograms, function(plot) {
        div(class = "grid-item histogram", as_widget(plot))
      })
    ),
    if (!is.null(individual_SNV_html) && !hide_ind_plots) {
      HTML(individual_SNV_html)
    },
    HTML(toggle_histogram_script),
    HTML(title_popup_script)
  )
    )

    # saving the final html file with 3D plots, histograms, and individual sceSNV plots
    dir.name <- getwd()
    output_file <- paste0(dir.name, '/Exploratory_Combined_Plots.html')
    save_html(grid_html, file = output_file)





}
