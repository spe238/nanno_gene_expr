# Figure 1a + b

ls(all.names = TRUE)
rm(list = ls(all.names = TRUE))

fl_growth <- read.csv2("fl-growth.csv")
etr <- read.csv2("ETR.csv")

library(reshape)
library(plotrix)
library(sciplot)

fl_melt <- melt(fl_growth, id=c("Time"))
fl_bind <- cbind(fl_melt, colsplit(fl_melt$variable, split="[_]", names=c("Sample","Replicate")))

etr_melt <- melt(etr, id=c("I"))
cols <- c("#9c5514", "#7f2a55", "#569e15", "#138f93")

etr_table <- data.frame(c("0.158 a", "0.143 a", "0.171 a", "0.177 a", "0.000", "yes"),
                        c("0.353 b", "0.358 b", "0.415 a", "0.386 ab", "0.000", "yes"),
                        c("51.563 b", "60.950 ab", "72.093 a", "53.850 b", "0.000", "yes"),
                        c("433.658 a", "477.586 a", "459.726 a", "400.771 a", "0.001", "yes"))
rownames(etr_table) <- c("CL", "FL5", "FL50", "FL500", "Pr > F", "significant")
colnames(etr_table) <- c(expression(alpla), "Fv'/Fm'", expression(ETR[max]), expression(E[k]))

par(mfrow=c(1,2), xpd = NA)
lineplot.CI(x.factor = fl_bind$Time, response = fl_bind$value, group = fl_bind$Sample, x.cont = TRUE, legend = TRUE, ncol = 1, pch = c(19,17,15,18),
            x.leg = 0, y.leg = 2.0, lty = 1, col = cols, type = "b", err.col = "black",
            ci.fun = function(x) {c(mean(x)-sd(x),mean(x)+sd(x))}, xlim = c(0,10), ylim = c(0,2), cex = 2, cex.leg = 1.3, cex.axis = 1.5,
            cex.lab = 1.5, err.width = 0.04, xaxp = c(0,10,5), lwd = 2, las = 1, bty = "n",
            xlab = "Time (d)", ylab = "Biomass (g/L)", fixed = TRUE, xaxt = "n")
axis(side = 1, cex.axis = 1.5)
legend(legend = "a", x = -3.5, y = 2.4, bty = "n", cex = 2.5)
legend(legend = "b", x = 10.0, y = 2.4, bty = "n", cex = 2.5)
text(x = 7, y = 1.55, expression(atop(µ*" = "*0.146~d^-1, R^2*" = 0.985")), pos = 4)
text(x = 7, y = 1.00, expression(atop(µ*" = "*0.126~d^-1, R^2*" = 0.995")), pos = 4)
text(x = 7, y = 0.75, expression(atop(µ*" = "*0.057~d^-1, R^2*" = 0.960")), pos = 4)
text(x = 7, y = 0.20, expression(atop(µ*" = "*0.025~d^-1, R^2*" = 0.870")), pos = 4)
lineplot.CI(x.factor = etr_melt$I, response = etr_melt$value, group = etr_melt$variable, x.cont = TRUE, legend = TRUE, ncol = 1, pch = c(19,17,15,18),
            x.leg = 50, y.leg = 70, lty = 1, col = cols, type = "b", err.col = "black",
            ci.fun = function(x) {c(mean(x)-sd(x),mean(x)+sd(x))}, xlim = c(44,1572), ylim = c(0,70), cex = 2, cex.leg = 1.3, cex.axis = 1.5,
            cex.lab = 1.5, err.width = 0.04, lwd = 2, las = 1, bty = "n",
            xlab = expression(Light~intensity~(mol~m^-2~s^-1)), ylab = "ETR", fixed = TRUE, xaxt = "n")
axis(side = 1, cex.axis = 1.5)
addtable2plot(250, 5, etr_table, bty = "o", hlines=FALSE, display.rownames = TRUE, display.colnames = FALSE, xpad = 0.2, ypad = 1.5)
text(x = 600, y = 30.5, pos = 4, expression(alpha*"          Fv'/Fm'       "*ETR[max]~~~~~~~~~~~~~~E[k]))
par(mfrow=c(1,1))

# Figure 1c

ls(all.names = TRUE)
rm(list = ls(all.names = TRUE))

npq <- read.csv2("NPQ.csv", sep = ",")
colnames(npq) <- c("Treatment", "Replicate", 0, 44, 132, 261, 421, 610, 823, 1188, 1572)

library(reshape)
library(sciplot)

npq_melt <- melt(npq, id=c("Treatment", "Replicate"))
npq_melt$value <- as.numeric(npq_melt$value)

cols <- c("#9c5514", "#7f2a55", "#569e15", "#138f93")

lineplot.CI(x.factor = npq_melt$variable, response = npq_melt$value, group = npq_melt$Treatment, x.cont = TRUE, legend = TRUE, ncol = 1, pch = c(19,17,15,18),
            x.leg = 0, y.leg = 1.0, lty = 1, col = cols, type = "b", err.col = "black",
            ci.fun = function(x) {c(mean(x)-sd(x),mean(x)+sd(x))}, xlim = c(0,2000), ylim = c(0,1.0), cex = 2, cex.leg = 1.3, cex.axis = 1.5,
            cex.lab = 1.5, err.width = 0.04, xaxp = c(0,1500,5), lwd = 2, las = 1, bty = "n",
            xlab = expression(PAR~(µmol~photons~m^-2~s^-1)), ylab = "NPQ (-)", fixed = TRUE, xaxt = "n")
axis(side = 1, cex.axis = 1.5)
legend(legend = "c", x = -300.5, y = 1.25, bty = "n", cex = 2.5, xpd = NA)

# Figure 3

ls(all.names = TRUE)
rm(list = ls(all.names = TRUE))

library(reshape)

gene_data <- read.csv2(file = "all_genes_data_feb20.csv", sep = ",")

colnames(gene_data) <- c("Gene", "CL", "FL5", "FL50", "FL500")

for (i in 1:(nrow(gene_data)/3)) {gene_data$CL_mean[i*3] <- mean(gene_data$CL[(3*i-2):(3*i)])}
for (i in 1:(nrow(gene_data)/3)) {gene_data$CL_mean[i*3-1] <- mean(gene_data$CL[(3*i-2):(3*i)])}
for (i in 1:(nrow(gene_data)/3)) {gene_data$CL_mean[i*3-2] <- mean(gene_data$CL[(3*i-2):(3*i)])}

gene_data$CL_norm <- gene_data$CL/gene_data$CL_mean
gene_data$FL5_norm <- gene_data$FL5/gene_data$CL_mean
gene_data$FL50_norm <- gene_data$FL50/gene_data$CL_mean
gene_data$FL500_norm <- gene_data$FL500/gene_data$CL_mean

gene_norm <- subset(gene_data, select = c(Gene, CL_norm, FL5_norm, FL50_norm, FL500_norm))
colnames(gene_norm) <- c("Gene", "CL", "FL5", "FL50", "FL500")

sub_photo <- droplevels(subset(gene_norm, Gene == "VCP1" |
                                 Gene == "VDE1" |
                                 Gene == "VDE2" |
                                 Gene == "ZE" |
                                 Gene == "LCYB" |
                                 Gene == "ALAD" |
                                 Gene == "POR1"))

sub_fatty_syn <- droplevels(subset(gene_norm, Gene == "BC" |
                                     Gene == "MCT" |
                                     Gene == "KAS" |
                                     Gene == "HAD" |
                                     Gene == "KAR1" |
                                     Gene == "KAR2" |
                                     Gene == "ENR1" |
                                     Gene == "ENR2"))

sub_fatty_traf <- droplevels(subset(gene_norm, Gene == "LCS1" |
                                      Gene == "LCS2" |
                                      Gene == "O6FADes-1" |
                                      Gene == "D6S"))

sub_glycero <- droplevels(subset(gene_norm, Gene == "GPAT" |
                                   Gene == "LPAT1" |
                                   Gene == "LPAT2" |
                                   Gene == "PAP1" |
                                   Gene == "PAP2" |
                                   Gene == "PAP3"))

sub_DGAT <- droplevels(subset(gene_norm, Gene == "DGAT1" |
                                Gene == "DGAT2" |
                                Gene == "DGAT3" |
                                Gene == "DGAT4" |
                                Gene == "DGAT5"))

sub_GPD <- droplevels(subset(gene_norm, Gene == "GPD1" |
                               Gene == "GPD2" |
                               Gene == "GPD3"))

sub_starch <- droplevels(subset(gene_norm, Gene == "AGPase1" |
                                  Gene == "AGPase2"))

sub_nitro <- droplevels(subset(gene_norm, Gene == "NAR1" |
                                 Gene == "NR"))

melt_1 <- melt(sub_photo, id=c("Gene"))
melt_2 <- melt(sub_fatty_syn, id=c("Gene"))
melt_3 <- melt(sub_fatty_traf, id=c("Gene"))
melt_4 <- melt(sub_glycero, id=c("Gene"))
melt_5 <- melt(sub_DGAT, id=c("Gene"))
melt_6 <- melt(sub_GPD, id=c("Gene"))
melt_7 <- melt(sub_starch, id=c("Gene"))
melt_8 <- melt(sub_nitro, id=c("Gene"))

all_sub <- list(melt_1, melt_2, melt_3, melt_4, melt_5, melt_6, melt_7, melt_8)

main <- c("Photosynthesis and photoprotection",
          "Fatty acid biosynthesis",
          "Fatty acid trafficking and desaturation",
          "Glycerolipid synthesis",
          "DGAT genes",
          "GPD genes",
          "Starch biosynthesis",
          "Nitrogen transport and assimilation")

sub_fig <- c("a", "b", "c", "d", "e", "f", "g", "h")
marg <- c(4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1)
space <- c(0.6, 0.6, 2, 1, 1, 2, 3, 3)
boxwex <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6)
cols <- c("darkgoldenrod4", "#7f2a55", "#569e15", "#138f93")

par(mfrow=c(4,2), xpd = NA)
for (i in 1:8) {
  (myplot <- boxplot(value ~ variable*Gene, data=all_sub[[i]], 
                     ylab="Fold change gene expression", xlab = "",
                     main=main[i], ylim = c(0.4, 3.0), xlim = c(1, 30),
                     col=cols, xaxt="n", frame.plot = FALSE, whisklty = "solid", medlwd = 1, boxwex = boxwex[i]))
  for(j in seq(4.5, nrow(all_sub[[i]])/3, 4)){ 
    abline(v=j,lty=1, col="grey", xpd = FALSE)
  }
  my_names <- sapply(strsplit(myplot$names, '\\.'), function(x) x[[2]])
  my_names <- my_names[seq(1, length(my_names), 4)]
  axis(1, at = seq(2.5, nrow(all_sub[[i]])/3, 4), labels = my_names, tick=FALSE , cex=0.3)
  text(-3.0, 3.5, labels = sub_fig[i], cex = 2)
}
legend(x = -15, y = +0.05, fill = cols, legend = unique(all_sub[[i]]$variable), horiz = TRUE, bty = "n", cex = 1.3)
par(mfrow=c(1,1))

# Figure S2

ls(all.names = TRUE)
rm(list = ls(all.names = TRUE))

library(reshape)
library(sciplot)

gene_data <- read.csv2(file = "all_genes_data_feb20.csv", sep = ",")

colnames(gene_data) <- c("Gene", "CL", "FL5", "FL50", "FL500")

sub_photo <- droplevels(subset(gene_data, Gene == "VCP1" |
                                 Gene == "VDE1" |
                                 Gene == "VDE2" |
                                 Gene == "ZE" |
                                 Gene == "LCYB" |
                                 Gene == "ALAD" |
                                 Gene == "POR1"))

sub_fatty_syn <- droplevels(subset(gene_data, Gene == "BC" |
                                     Gene == "MCT" |
                                     Gene == "KAS" |
                                     Gene == "HAD" |
                                     Gene == "KAR1" |
                                     Gene == "KAR2" |
                                     Gene == "ENR1" |
                                     Gene == "ENR2"))

sub_fatty_traf <- droplevels(subset(gene_data, Gene == "LCS1" |
                                      Gene == "LCS2" |
                                      Gene == "O6FADes-1" |
                                      Gene == "D6S"))

sub_glycero <- droplevels(subset(gene_data, Gene == "GPAT" |
                                   Gene == "LPAT1" |
                                   Gene == "LPAT2" |
                                   Gene == "PAP1" |
                                   Gene == "PAP2" |
                                   Gene == "PAP3"))

sub_DGAT <- droplevels(subset(gene_data, Gene == "DGAT1" |
                                Gene == "DGAT2" |
                                Gene == "DGAT3" |
                                Gene == "DGAT4" |
                                Gene == "DGAT5"))

sub_GPD <- droplevels(subset(gene_data, Gene == "GPD1" |
                               Gene == "GPD2" |
                               Gene == "GPD3"))

sub_starch <- droplevels(subset(gene_data, Gene == "AGPase1" |
                                  Gene == "AGPase2"))

sub_nitro <- droplevels(subset(gene_data, Gene == "NAR1" |
                                 Gene == "NR"))

melt_1 <- melt(sub_photo, id=c("Gene"))
melt_2 <- melt(sub_fatty_syn, id=c("Gene"))
melt_3 <- melt(sub_fatty_traf, id=c("Gene"))
melt_4 <- melt(sub_glycero, id=c("Gene"))
melt_5 <- melt(sub_DGAT, id=c("Gene"))
melt_6 <- melt(sub_GPD, id=c("Gene"))
melt_7 <- melt(sub_starch, id=c("Gene"))
melt_8 <- melt(sub_nitro, id=c("Gene"))

all_sub <- list(melt_1, melt_2, melt_3, melt_4, melt_5, melt_6, melt_7, melt_8)

main <- c("Photosynthesis and photoprotection",
          "Fatty acid biosynthesis",
          "Fatty acid trafficking and desaturation",
          "Glycerolipid synthesis",
          "DGAT genes",
          "GPD genes",
          "Starch biosynthesis",
          "Nitrogen transport and assimilation")

sub_fig <- c("a", "b", "c", "d", "e", "f", "g", "h")
marg <- c(4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1)
space <- c(0.6, 0.6, 2, 1, 1, 2, 3, 3)
width <- c(1, 0.8, 1, 1, 1, 1, 1, 1)
cols <- c("darkgoldenrod4", "#7f2a55", "#569e15", "#138f93")

par(mfrow=c(4,2), xpd = NA)
for (i in 1:8) {
  par(mar=c(5.1, marg[i], 4.1, 2.1))
  bargraph.CI(x.factor = Gene,
              response = value,
              group = variable, data = all_sub[[i]],
              col = cols,
              xlab = "",
              ylab = "Relative gene expression (-)",
              main = "",
              ylim = c(0.0, 1.6),
              xlim = c(1,31),
              ci.fun = function(x) {c(if (mean(x) < 0) {mean(x)-sd(x)} else {mean(x)}, if (mean(x) > 0) {mean(x)+sd(x)} else {mean(x)})},
              las = 1,
              lc = FALSE,
              uc = TRUE,
              cex.names = 1.3,
              cex.axis = 1.3,
              cex.lab = 1.3,
              err.width = 0.04,
              width = width[i],
              space = c(0, space[i]))
  text(-3.0, 1.8, labels = sub_fig[i], cex = 2)
  text(9, 1.8, labels = main[i], cex = 1.5)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}
legend(x = -15, y = -0.4, fill = cols, legend = unique(all_sub[[i]]$variable), horiz = TRUE, bty = "n", cex = 1.3)
par(mfrow=c(1,1))

# Figure S4

ls(all.names = TRUE)
rm(list = ls(all.names = TRUE))

library(reshape)
library(sciplot)

gene_data <- read.csv2(file = "all_genes_data_feb20.csv", sep = ",")

colnames(gene_data) <- c("Gene", "CL", "FL5", "FL50", "FL500")

s4 <- droplevels(subset(gene_data, Gene == "GAPDH" |
                          Gene == "BT" |
                          Gene == "AT1" |
                          Gene == "ACT" |
                          Gene == "TEF1" |
                          Gene == "UFP" |
                          Gene == "CYP"))

melt_9 <- melt(s4, id=c("Gene"))

cols <- c("darkgoldenrod4", "#7f2a55", "#569e15", "#138f93")

bargraph.CI(x.factor = Gene,
            response = value,
            group = variable, data = melt_9,
            col = cols,
            xlab = "",
            ylab = "Relative gene expression (-)",
            main = "",
            ylim = c(0.0, 1.6),
            xlim = c(1,31),
            fun = function(x) {mean(x)},
            ci.fun = function(x) {c(if (mean(x) < 0) {mean(x)-sd(x)} else {mean(x)}, if (mean(x) > 0) {mean(x)+sd(x)} else {mean(x)})},
            las = 1,
            lc = TRUE,
            uc = TRUE,
            cex.names = 1.3,
            cex.axis = 1.3,
            cex.lab = 1.3,
            err.width = 0.04,
            width = 1,
            space = c(0,0.6))
legend(x = 15, y = 1.5, fill = cols, legend = unique(melt_9$variable), horiz = TRUE, bty = "n", cex = 1.3)

# Tukey Test

ls(all.names = TRUE)
rm(list = ls(all.names = TRUE))

library(multcompView)
library(reshape)

gene_data <- read.csv2(file = "all_genes_data_feb20.csv", sep = ",")
colnames(gene_data) <- c("Gene", "CL", "FL5", "FL50", "FL500")

data <- melt(gene_data, id=c("Gene"))

model=lm(data$value ~ data$variable)
ANOVA=aov(model)

TUKEY <- TukeyHSD(x=ANOVA, 'data$variable', conf.level=0.95)

plot(TUKEY , las=1 , col="brown")

generate_label_df <- function(TUKEY, variable){
  
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  Tukey.labels$variable=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$variable) , ]
  return(Tukey.labels)
}

LABELS <- generate_label_df(TUKEY , "data$variable")

my_colors <- c("darkgoldenrod4", "#7f2a55", "#569e15", "#138f93")

a <- boxplot(data$value ~ data$variable, ylim=c(min(data$value), 1.1*max(data$value)), col=my_colors, ylab="", main="Tukey test", xlab = "Treatment", outline = FALSE,
             boxwex=0.2)
over <- 0.1*max(a$stats[nrow(a$stats),] )
text(c(1:nlevels(data$variable)), a$stats[nrow(a$stats),]+over, LABELS[,1],)
