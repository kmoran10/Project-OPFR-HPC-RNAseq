
# specific gene investigations 



library(tidyverse)
library(lme4)
library(lmerTest)
library(WGCNA)


TREAT <- readRDS("results/TREATMENT_limma_results.RDS")
males <- readRDS("results/Elevel_M_limma_results.RDS")
low <- readRDS("results/Elevel_low_limma_results.RDS")
high <- readRDS("results/Elevel_high_limma_results.RDS")
oil <- readRDS("results/Oil_limma_results.RDS")
opfr <- readRDS("results/OPFR_limma_results.RDS")


# gaba related genes
# gaba a R
TREAT %>% filter(symbol == "Gabra1") # no
TREAT %>% filter(symbol == "Gabra2") # yes

oil$M_v_low %>% filter(symbol == "Gabra2") #no
oil$M_v_high %>% filter(symbol == "Gabra2") #no
oil$low_v_high %>% filter(symbol == "Gabra2") #no
opfr$M_v_low %>% filter(symbol == "Gabra2") #no
opfr$M_v_high %>% filter(symbol == "Gabra2") #yes 
opfr$low_v_high %>% filter(symbol == "Gabra2") #close --- so, Gabra2 is esp higher in higher e-level

TREAT %>% filter(symbol == "Gabra3") # no 
TREAT %>% filter(symbol == "Gabra4") # no
TREAT %>% filter(symbol == "Gabra5") # no
TREAT %>% filter(symbol == "Gabrb1") # no 
TREAT %>% filter(symbol == "Gabrb2") # no 
TREAT %>% filter(symbol == "Gabrb3") # no
TREAT %>% filter(symbol == "Gabrd") #  no
TREAT %>% filter(symbol == "Gabrg1") # no
TREAT %>% filter(symbol == "Gabrg2") # no
TREAT %>% filter(symbol == "Gabrg3") # no
TREAT %>% filter(symbol == "Gabrr2") # no
# gaba b R
TREAT %>% filter(symbol == "Gabbr1") # no
TREAT %>% filter(symbol == "Gabbr2") # no
# gaba other - only 
TREAT %>% filter(symbol == "Gabarap") # no
TREAT %>% filter(symbol == "Gad1") # no
TREAT %>% filter(symbol == "Gad2") # no
TREAT %>% filter(symbol == "Gpr37") # no
TREAT %>% filter(symbol == "Slc32a1") # no

oil$M_v_low %>% filter(symbol == "Slc32a1") # close -- but not in high. so not really.
oil$M_v_high %>% filter(symbol == "Slc32a1") #no
oil$low_v_high %>% filter(symbol == "Slc32a1") #no
opfr$M_v_low %>% filter(symbol == "Slc32a1") #no
opfr$M_v_high %>% filter(symbol == "Slc32a1") #no
opfr$low_v_high %>% filter(symbol == "Slc32a1") #no


# oxytocin?
TREAT %>% filter(symbol == "Oxt") # no
TREAT %>% filter(symbol == "Oxtr") # no
TREAT %>% filter(symbol == "Cd38") # no

oil$M_v_low %>% filter(symbol == "Oxtr") #no
oil$M_v_high %>% filter(symbol == "Oxtr") #no
oil$low_v_high %>% filter(symbol == "Oxtr") #no
opfr$M_v_low %>% filter(symbol == "Oxtr") #no
opfr$M_v_high %>% filter(symbol == "Oxtr") #no
opfr$low_v_high %>% filter(symbol == "Oxtr") #no

oil$M_v_low %>% filter(symbol == "Cd38") #no
oil$M_v_high %>% filter(symbol == "Cd38") #no
oil$low_v_high %>% filter(symbol == "Cd38") #no
opfr$M_v_low %>% filter(symbol == "Cd38") #no
opfr$M_v_high %>% filter(symbol == "Cd38") #no
opfr$low_v_high %>% filter(symbol == "Cd38") #no

# estrogen R
TREAT %>% filter(symbol == "Esrra") # no
TREAT %>% filter(symbol == "Esrrb") # no

oil$M_v_low %>% filter(symbol == "Esrra") #no
oil$M_v_high %>% filter(symbol == "Esrra") #no
oil$low_v_high %>% filter(symbol == "Esrra") #no
## kinda weird that in oil animals, there is no diff in Esrra between e-level??


oil$M_v_low %>% filter(symbol == "Esrrb") #no
oil$M_v_high %>% filter(symbol == "Esrrb") #no
oil$low_v_high %>% filter(symbol == "Esrrb") #no
## kinda weird that in oil animals, there is no diff in Esrr between e-level??


opfr$M_v_low %>% filter(symbol == "Esrra") #no
opfr$M_v_high %>% filter(symbol == "Esrra") #no
opfr$low_v_high %>% filter(symbol == "Esrra") #no


#other 
TREAT %>% filter(symbol == "Hsd11b2") # no
TREAT %>% filter(symbol == "Sod2") # no
