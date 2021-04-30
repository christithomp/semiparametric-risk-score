#####################################################################
# This program combines the images to both approaches into one image.
#####################################################################
rm(list = ls())

library(ggplot2)
library(gridExtra)

# Load estimates and derivative plots
load("./final_code/results/est_deriv_plotlist_PL.rda")
load("./final_code/results/est_deriv_plotlist_jags.rda")
glist <- c(glist_PL, glist_jags)
fdlist <- c(fdlist_PL, fdlist_jags)

# Load simulation plots
load("./final_code/results/simplot_PL_list.rda")
load("./final_code/results/simplot_jags_list.rda")
simplot <- c(simplot_PL, simplot_jags)

# Plot estimate results together
#pdf(file = "./final_code/figures/Estimates_plot.pdf", width = 10, height = 9)
tiff("./final_code/figures/CuiFigure1.tiff", units="in", width=10, height=9, res=800)
do.call(grid.arrange, c(glist, nrow = 2))
dev.off()

# Plot first derivative results together
#pdf(file = "./final_code/figures/First_deriv_plot.pdf", width = 10, height = 9)
tiff("./final_code/figures/CuiFigure2.tiff", units="in", width=10, height=9, res=800)
do.call(grid.arrange, c(fdlist, nrow = 2))
dev.off()

# Plot simulation results together
#pdf(file = "./final_code/figures/simulation_plot.pdf", width = 10, height = 10)
tiff("./final_code/figures/CuiFigure3.tiff", units="in", width=10, height=10, res=800)
do.call(grid.arrange, c(simplot, nrow = 2))
dev.off()