library(tidyverse)

gppdi_flags <- read_csv("GPPDI_617/data/GPPDI_ClassC_Flags_5164.csv")
gppdi_data <- read_csv("GPPDI_617/data/GPPDI_ClassC_NPP_5164.csv")

readLines("GPPDI_617/data/GPPDI_ClassC_NPP_5164.csv", 10)

gppdi_data
