# AMOC


## Introduction
This is a selection of the most relevant files from my thesis project, which assessed the trend in observations of the AMOC (Atlantic Meridional Overturning Circulation)
strength using climate simulations to represent natural variability. Note, this collection would not produce the expected output by themselves; instead I'm keeping some 
scripts in a private repository. The available scripts are more for illustrative purposes.

## Data Access
The data used is all freely-available. Observations of the AMOC are found at rapid.ac.uk. Simulated data is available at esgf-node.llnl.gov.


## Data Extraction & Cleaning script(s)
- function_make_AMOC_mat_file.m     -- extracts AMOC strength at 26N from downloaded NetCDF files, cleans, and saves to .mat format

## Data Transformation script(s)
- function_x12_filter.m             -- extracts the interannual component of the timeseries, discarding the annual cycle and subannual variation
- function_apply_henderson.m        -- helps to identify the interannual component
- function_apply_seasonal_filter.m  -- identifies the annual cycle

## Data Analysis script(s)
- function_trend_id_pdf.m            -- breaks multi-century timeseries into segments of a given length, finds the linear trend in those segments, and 
                                              fits a PDF (probability distribution function)
- function_sig_vs_lgth.m             -- finds the significance of the observed trend relative to the trend distribution of segments in simulations
- script_make_PDFs_for_tableau.m     -- focuses on coarsening the results (i.e. segment lengths in multiples of 5) to make Tableau interactivity more responsive
- analysis_spectra_of_timeseries.m   -- finds and displays the spectra of simulated and observed timeseries

## Data Visualisation script(s)
- thesis_figure_AMOC_timeseries.m                    -- 4 panel figure showing array location, frequency components of observations etc.
- thesis_figure_PDF_scenario_comparison.m            -- 2 panel figure showing variability of simulations, and difference between pre-industrial and present-day
- thesis_figure_PDFs_for_changing_trend_durations.m  -- 2 panel figure showing how segment length impacts how long a trend must persist to become significant
- thesis_figure_duration_for_significance.m          -- 2 panel figure showing the distribution of long- and short-term variability in models, and how PDFs narrow

