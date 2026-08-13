library(here)
library(deSolve)
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)
library(readr)
library(plotly)
library(zoo)
library(lubridate)
library(BayesianTools)
library(posterior)
library(HDInterval)
library(viridis)
library(socialmixr)
library(readxl)
library(parallel)
library(tibble)
library(patchwork)
library(Rcpp)
library(coda)

# helper for pathogen colours
pathogen_cols <- c(
  # Scottish and RCGP spellings
  "Respiratory syncytial virus"= "#FABA39FF", "RSV"= "#FABA39FF",
  "Influenza" = "#36AAF9FF",
  "Influenza - Type A (any subtype)" = "#1AE4B6FF", "Influenza A" = "#1AE4B6FF",
  "Influenza - Type B" = "#72FE5EFF", "Influenza B" = "#72FE5EFF",
  "Rhinovirus" = "#F66B19FF", "Rhino"  = "#F66B19FF", "Entero- or Rhinovirus"= "#F66B19FF",
  "Entero"= "#7A0403FF", 
  "Adenovirus" = "#30123BFF", "Adeno" = "#30123BFF",
  "Human metapneumovirus" = "#4662D7FF", "hMPV" = "#4662D7FF", "HMPV" = "#4662D7FF",
  "Parainfluenza virus" = "#C7EF34FF", "Parainfluenza (Any Type)" = "#C7EF34FF",
  "Seasonal coronavirus (Non-SARS-CoV-2)" = "#CB2A04FF", "OtherCoronavirus" = "#CB2A04FF", "Seasonal coronavirus" = "#CB2A04FF", "Other Coronavirus" = "#CB2A04FF"
)
