# Data Analisis for the CYCLIN A1 IN pT1 UC project
# Opening the Dataset
ca1 <- read.csv("ca1.csv")
# Attaching the Dataset
attach(ca1)
# Loading the Libraries that will be used
library(gmodels)  # For frequencies and contingency tables (CrossTable)
library(survival) # For Survival Analyses
# Describing the Data
CrossTable(pt.biopsy, digits = 0, format = c("SPSS")) # pT1 at Biopsy
CrossTable(cis.biopsy, digits = 0, format = c("SPSS")) # CIS at Biopsy
CrossTable(treatment, digits = 0, format = c("SPSS")) # Initial Treatment
CrossTable(recurrence.outcome, digits = 0, format = c("SPSS")) # Recurrence ratio
CrossTable(dx.final, digits = 0, format = c("SPSS")) # Distribution of TNM stages
describe(recurrence.number) # Describing the number of recurrence episodes
IQR(recurrence.number, na.rm = TRUE) # Estimates the IQR of recurrence episodes
CrossTable(recurrence.groups, digits = 0, format = c("SPSS")) # Recurrence episodes by groups
# Cyclin A1 Expression and Outcome
describe(ca1.extension) # Describing the % of cyclin A1 expression
IQR(ca1.extension, na.rm = TRUE) # Estimates the IQR of % of cyclin A1 expression
CrossTable(ca1.positive.10, digits = 0, format = c("SPSS")) # % Cyclin A1 Positivity
# Table 1
CrossTable(pt.biopsy, ca1.positive.10, digits = 0, fisher = TRUE, format = c("SPSS")) # Cyclin A1 by pT at Biopsy
CrossTable(cis.biopsy, ca1.positive.10, digits = 0, fisher = TRUE, format = c("SPSS")) # Cyclin A1 by CIS at Biopsy
CrossTable(recurrence, ca1.positive.10, digits = 0, fisher = TRUE, format = c("SPSS")) # Cyclin A1 by Recurrence
CrossTable(progression, ca1.positive.10, digits = 0, fisher = TRUE, format = c("SPSS")) # Cyclin A1 by Progression
1-(1-.05)^(1/4) # Adjusting P value using Sidak's correction
# Defining Survival Objects
surv.recurrence <- with(ca1, Surv(time.event, recurrence)) # For tumor recurrence
surv.progression <- with(ca1, Surv(time.event, progression)) # For tumor progression
# Estimating Survival Functions for Tumor Recurrence
pt1.recurrence <- survfit(surv.recurrence ~ pt.biopsy, data = ca1) # pT1 and recurrence
cis.recurrence <- survfit(surv.recurrence ~ cis.biopsy, data = ca1) # CIS and recurrence
treatment.recurrence <- survfit(surv.recurrence ~ treatment, data = ca1) # Initial treatment and recurrence
ca1.recurrence <- survfit(surv.recurrence ~ ca1.positive.10, data = ca1) # Cyclin A1 and recurrence
# Figure 2A
survdiff(surv.recurrence ~ pt.biopsy) # Log-rank test
tiff("fig2a.tiff", 600, 450, pointsize = 16, compression = c("zip"))
  plot(pt1.recurrence, main = "pT at Biopsy", cex.main = 1.75, col =c(1,2), mark = c(2,0), lty = c(2,1))
  text("P (log-rank test) = 0.02", x = 17, y = 0.01)
  legend("topright", c("pT1a", "pT1b"), pch = c(2,0), lty = c(2,1), col = c(1,2), bty = "n")
dev.off()
# Figure 2B
survdiff(surv.recurrence ~ cis.biopsy) # Log-rank test
tiff("fig2b.tiff", 600, 450, pointsize = 16, compression = c("zip"))
  plot(cis.recurrence, main = "CIS at Biopsy", cex.main = 1.75, col =c(1,2), mark = c(2,0), lty = c(2,1))
  text("P (log-rank test) = 0.32", x = 17, y = 0.01)
  legend("topright", c("Absent", "Present"), pch = c(2,0), lty = c(2,1), col = c(1,2), bty = "n")
dev.off()
# Figure 2C
survdiff(surv.recurrence ~ treatment) # Log-rank test
tiff("fig2c.tiff", 600, 450, pointsize = 16, compression = c("zip"))
  plot(treatment.recurrence, main = "Initial Treatment", cex.main = 1.75, col =c(1,2,4), mark = c(2,0,5), lty = c(2,1,3))
  text("P (log-rank test) = 0.47", x = 17, y = 0.01)
  legend("topright", c("TURB+BCG", "TURB+MC", "TURB alone"), pch = c(2,0,5), lty = c(2,1,3), col = c(1,2,4), bty = "n")
dev.off()
# Figure 2D
survdiff(surv.recurrence ~ ca1.positive.10) # Log-rank test
tiff("fig2d.tiff", 600, 450, pointsize = 16, compression = c("zip"))
  plot(ca1.recurrence, main = "Cyclin A1", cex.main = 1.75, col =c(1,2), mark = c(2,0), lty = c(2,1))
  text("P (log-rank test) = 0.68", x = 17, y = 0.01)
  legend("topright", c("Negative", "Positive"), pch = c(2,0), lty = c(2,1), col = c(1,2), bty = "n")
dev.off()
# Estimating Survival Functions for Tumor Progression
pt1.progression <- survfit(surv.progression ~ pt.biopsy, data = ca1) # pT1 and progression
cis.progression <- survfit(surv.progression ~ cis.biopsy, data = ca1) # CIS and progression
treatment.progression <- survfit(surv.progression ~ treatment, data = ca1) # Initial treatment and progression
ca1.progression <- survfit(surv.progression ~ ca1.positive.10, data = ca1) # Cyclin A1 and progression
# Figure 2E
survdiff(surv.progression ~ pt.biopsy) # Log-rank test
tiff("fig2e.tiff", 600, 450, pointsize = 16, compression = c("zip"))
  plot(pt1.progression, main = "pT at Biopsy", cex.main = 1.75, col =c(1,2), mark = c(2,0), lty = c(2,1))
  text("P (log-rank test) = 1.11e-06", x = 19, y = 0.01)
  legend("topright", c("pT1a", "pT1b"), pch = c(2,0), lty = c(2,1), col = c(1,2), bty = "n")
dev.off()
# Figure 2F
survdiff(surv.progression ~ cis.biopsy) # Log-rank test
tiff("fig2f.tiff", 600, 450, pointsize = 16, compression = c("zip"))
  plot(cis.progression, main = "CIS at Biopsy", cex.main = 1.75, col =c(1,2), mark = c(2,0), lty = c(2,1))
  text("P (log-rank test) = 0.03", x = 17, y = 0.01)
  legend("bottomright", c("Absent", "Present"), pch = c(2,0), lty = c(2,1), col = c(1,2), bty = "n")
dev.off()
# Figure 2G
survdiff(surv.progression ~ treatment) # Log-rank test
tiff("fig2g.tiff", 600, 450, pointsize = 16, compression = c("zip"))
  plot(treatment.progression, main = "Initial Treatment", cex.main = 1.75, col =c(1,2,4), mark = c(2,0,5), lty = c(2,1,3))
  text("P (log-rank test) = 0.47", x = 17, y = 0.01)
  legend("topright", c("TURB+BCG", "TURB+MC", "TURB alone"), pch = c(2,0,5), lty = c(2,1,3), col = c(1,2,4), bty = "n")
dev.off()
# Figure 2H
survdiff(surv.progression ~ ca1.positive.10) # Log-rank test
tiff("fig2h.tiff", 600, 450, pointsize = 16, compression = c("zip"))
  plot(ca1.progression, main = "Cyclin A1", cex.main = 1.75, col =c(1,2), mark = c(2,0), lty = c(2,1))
  text("P (log-rank test) = 0.004", x = 17, y = 0.01)
  legend("topright", c("Negative", "Positive"), pch = c(2,0), lty = c(2,1), col = c(1,2), bty = "n")
dev.off()
# Estimating Hazard Ratios for Tumor Recurrence
cox.pt1.recurrence <- coxph(surv.recurrence ~ pt.biopsy, data = ca1) # pT1 and recurrence
cox.cis.recurrence <- coxph(surv.recurrence ~ cis.biopsy, data = ca1) # CIS and recurrence
cox.treatment.recurrence <- coxph(surv.recurrence ~ treatment, data = ca1) # Initial treatment and recurrence
cox.ca1.recurrence <- coxph(surv.recurrence ~ ca1.positive.10, data = ca1) # Cyclin A1 and recurrence
cox.ca1pt1.recurrence <- coxph(surv.recurrence ~ ca1.positive.10 + pt.biopsy, data = ca1) # Cyclin A1 & pT1 at Biopsy and recurrence
# Estimating Hazard Ratios for Tumor Progression
cox.pt1.progression <- coxph(surv.progression ~ pt.biopsy, data = ca1) # pT1 and progression
cox.cis.progression <- coxph(surv.progression ~ cis.biopsy, data = ca1) # CIS and progression
cox.treatment.progression <- coxph(surv.progression ~ treatment, data = ca1) # Initial treatment and progression
cox.ca1.progression <- coxph(surv.progression ~ ca1.positive.10, data = ca1) # Cyclin A1 and progression
cox.ca1pt1.progression <- coxph(surv.progression ~ ca1.positive.10 + pt.biopsy, data = ca1) # Cyclin A1 & pT1 at Biopsy and progression
#Table 2
# Tumor Recurrence
summary(cox.pt1.recurrence) # pT1 and recurrence
summary(cox.cis.recurrence) # CIS and recurrence
summary(cox.treatment.recurrence) # Treatment and recurrence
summary(cox.ca1.recurrence) # Cyclin A1 and recurrence
summary(cox.ca1pt1.recurrence) # Cyclin A1 & pT1 and recurrence
# Tumor Progression
summary(cox.pt1.progression) # pT1 and progression
summary(cox.cis.progression) # CIS and progression
summary(cox.treatment.progression) # Treatment and progression
summary(cox.ca1.progression) # Cyclin A1 and progression
summary(cox.ca1pt1.progression) # Cyclin A1 & pT1 and progression
# Evaluating the interaction between pT1 and Cyclin A1
cox.ca1bypt1.progression <- coxph(surv.progression ~ ca1.positive.10*pt.biopsy, data = ca1) 
summary(cox.ca1bypt1.progression)
# Comparison of Model 1 (only pT1) vs. Model 2 (pT1 + Cyclin A1)
# vs. Model 3 (only Cyclin A1) for tumor progression
anova(cox.pt1.progression, cox.ca1pt1.progression) # Model 1 vs. Model 2
anova(cox.ca1.progression, cox.ca1pt1.progression) # Model 3 vs. Model 2