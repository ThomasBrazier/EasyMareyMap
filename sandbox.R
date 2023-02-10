source("R/mareyMap.R")
m = data.frame(set = "Arabidopsis",
               map = c(1, 1, 2, 3),
               mkr = c("a", "b", "c", "d"),
               phys = c(2, 5, 6, 9),
               gen = c(10, 20, 40, 56),
               vld = c(T, T, T, F))
m

df = mareyMap(m)
df

m = read.table("data/Arabidopsis_thaliana.txt", header = T)

df = mareyMap(m)
df

source("R/outlierSelection.R")
source("R/removeChromosome.R")

# df = outlierSelection(df, "1")
# df$mareyMap$vld
#
# df = removeChromosome(df, "1")

source("R/plot.R")
plot(df)
p = plot(df)
p



df = mareyMap(m)
source("R/summary.R")
summary(df)

source("R/recombinationMap.R")
source("R/calibrateSmoothing.R")
source("R/cvResampling.R")
source("R/fitLoess.R")
source("R/fitSpline.R")

res = recombinationMap(df, chromosome = "1", method = "loess")
res = recombinationMap(res, chromosome = "2", method = "loess")
res = recombinationMap(res, chromosome = "3", method = "spline")
summary(res)

source("R/recombinationPlot.R")
recombinationPlot(res)
