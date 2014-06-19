# Defining the survival.plot function for plotting survival curves
survival.plot <- function(x, title, p, position = "topright", levels, ...){
        plot(x, cex = 2, main = title, cex.main = 1.75,
             xlab = "Follow-Up (Months)", ylab = "Survival Function", cex.lab = 1.5,
             col =c(1,2,4,3), mark = c(2,0,5,1), lty = c(2,1,3,4))
        text(paste("P (log-rank test) = ", p), x = 17, y = 0.01, cex = 1.25)
        legend(x = position, legend = levels, pch = c(2,0,5,1), lty = c(2,1,3,4),
               col = c(1,2,4,3), bty = "n", cex = 1.25)
}
# x = must be a survival object
# title = corresponds to the plot's main title
# p = corresponds to the P value from comparing the curves
# position = corresponds to the legend position
# levels = corresponds to the curves denomination for the legend
# ... = additional arguments to be passed