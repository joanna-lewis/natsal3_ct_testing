rm(list=ls())

library(rstan)
library(foreign)
library(survey)
library(MASS)

setwd('/Users/Joanna/OneDrive - Imperial College London/backup/ct_trends/Natsal_symptoms_testing/modelling')

natsal3 <- read.dta("/Users/Joanna/OneDrive - Imperial College London/backup/Natsal-3/UKDA-7799-stata11/stata11/eul_natsal_2010_for_archive.dta")

natsal3$tested <- NA
natsal3$tested[ (natsal3$whnchlam == "Less than 1 year ago") | (natsal3$chlmtest == "Yes") ] <- TRUE
natsal3$tested[ natsal3$chlmtest == "No" ] <- FALSE
natsal3$tested[natsal3$totlife %in% c(0,9999)] <- NA

natsal3$why <- factor(NA, levels = levels(natsal3$chtstwy1))
natsal3$why[!(natsal3$chtstwy1 %in% c("Not applicable", "Not answered"))] <- natsal3$chtstwy1[!(natsal3$chtstwy1 %in% c("Not applicable", "Not answered"))]
natsal3$why[!(natsal3$chltstwy %in% c("Not applicable", "Not answered"))] <- natsal3$chltstwy[!(natsal3$chltstwy %in% c("Not applicable", "Not answered"))]
natsal3$why <- droplevels(natsal3$why)

natsal3$why2 <- NA
natsal3$why2[natsal3$why == "I had symptoms"] <- "symptoms"
natsal3$why2[natsal3$why %in% c("I was notified because a partner was diagnosed with Chlamydia")] <- "partner"
natsal3$why2[natsal3$why == "Check up after previous positive test"] <- "re-test"
natsal3$why2[natsal3$why %in% c("I wanted a general sexual health check-up", "I had no symptoms but I was worried about the risk of Chlamydia", "I was offered a routine test")] <- "screen"
natsal3$why2[natsal3$why %in% c("Other", "My partner had symptoms")] <- "other"
natsal3$why2 <- factor(natsal3$why2, levels = names(table(natsal3$why2))[order(table(natsal3$why2), decreasing=TRUE)], ordered = TRUE)

sub <- subset(natsal3, natsal3$totnewy3 %in% c("0 new het.&/or sam. partners in last year","1 new het.&/or sam. partner in last year", "2+ new het.&/or sam. partners in last year"))
sub$totnewy3 <- factor(as.character(sub$totnewy3))

init0 <- list(p_symp = 0.5, lambda_slow = 0.74, foi = rep(0.001, times=3), scr = rep(0.001, times=3), trt = 0.001)

# data set with stratification by number of partners
dt_str <- list(
  
  N = nrow(sub[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested),]),
  N_strata = length(unique(sub$totnewy3)),
  str = (as.numeric(sub$totnewy3[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)])),
  wt = sub$total_wt[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)],
  tested = sub$tested[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)],
  symp = sub$why2[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)] == "symptoms",
  partner = sub$why2[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)] == "partner",
  diag = sub$whnchlam[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)] == "Less than 1 year ago"
  
)

dt_str$symp[is.na(dt_str$symp)] <- -99
dt_str$partner[is.na(dt_str$partner)] <- -99

# data set without stratification by number of partners
dt_ustr <- dt_str
dt_ustr$str <- 1 + 0*(as.numeric(sub$totnewy3[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)]))

fit_str <-stan(
  file = '3-comp-model-stratified-PN-women.stan',
  data = dt_str,
  chains = 1,
  iter = 10000,
  warmup = 1000,
  init = list(init0),
  seed = 12345
)
op_str <- extract(fit_str)

fit_ustr <-stan(
  file = '3-comp-model-stratified-PN-women.stan',
  data = dt_ustr,
  chains = 1,
  iter = 10000,
  warmup = 1000,
  init = list(init0),
  seed = 12345
)
op_ustr <- extract(fit_ustr)


quartz(height = 9, width= 9 )
par(mfrow=c(3,3), mar=c(5,4,2,2))

h_s <- hist(100*op_str$p_symp, plot=FALSE)
h_u <- hist(100*op_ustr$p_symp, plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = 'Proportion of infections symptomatic (%)', ylab='Density', main='', xlim=c(0,100), ylim=c(0,0.12), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(seq(0,100,1), dbeta(seq(0,1,0.01), 27, 90)/100)
legend("topright", "A", bty="n", cex=3)

h_s <- hist(op_str$trt, plot=FALSE)
h_u <- hist(op_ustr$trt, plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = expression("Symptomatic treatment rate (" ~ year^{-1} ~ ")"), ylab = 'Density', main='', xlim=c(0,35), ylim=c(0,0.11), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(seq(0,100,0.1), dgamma(seq(0,100,0.1),14,1))
legend("topright", "B", bty="n", cex=3)

h_s <- hist(op_str$lambda_slow, plot=FALSE)
h_u <- hist(op_ustr$lambda_slow, plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = expression('Natural clearance rate (' ~ year^{-1} ~ ')'), ylab='Denisty', main='', xlim=c(0,1.5), ylim=c(0,6), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(seq(0,10,0.01), dnorm(seq(0,10,0.01), 0.74, (0.89-0.61)/3.919928))
legend("topright", "C", bty="n", cex=3)

h_s <- hist(op_str$pn[,1], plot=FALSE)
h_u <- hist(op_ustr$pn[,1], plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = expression('Partner notification rate (uninfected women) (' ~ year^{-1} ~ ')'), ylab='Denisty', main='', xlim=c(0,0.025), ylim=c(0,200), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(seq(0,1,0.001), dexp(seq(0,1,0.001),0.001))
legend("topright", "D", bty="n", cex=3)

h_s <- hist(op_str$pn[,2], plot=FALSE)
h_u <- hist(op_ustr$pn[,2], plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = expression('Partner notification rate (infected women) (' ~ year^{-1} ~ ')'), ylab='Denisty', main='', xlim=c(0,1), ylim=c(0,10), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(seq(0,1,0.001), dexp(seq(0,1,0.001),0.001))
legend("topright", "E", bty="n", cex=3)

plot.new()

plot(0,0,pch='',xlab = expression("Force of infection (" ~ year^{-1} ~ ")"), ylab='Density', main='', xlim=c(0,0.5), ylim=c(0,50), bty='n')
lines(seq(0,1,0.01), dexp(seq(0,1,0.01),0.001))
legend("topright", "F", bty="n", cex=3)
fh1 <- hist(op_str$foi[,1], plot=FALSE)
fh2 <- hist(op_str$foi[,2], plot=FALSE)
fh3 <- hist(op_str$foi[,3], plot=FALSE)
fh <- hist(op_ustr$foi[,1], plot=FALSE)
lines(rep(fh$breaks, each=2), c(0, rep(fh$density, each=2), 0), col='grey')
lines(rep(fh1$breaks, each=2), c(0, rep(fh1$density, each=2), 0), col='darkgreen')
lines(rep(fh2$breaks, each=2), c(0, rep(fh2$density, each=2), 0), col='blue')
lines(rep(fh3$breaks, each=2), c(0, rep(fh3$density, each=2), 0), col='red')

plot(0,0,pch='', xlab = expression("Screening rate (" ~ year^{-1} ~ ")"), ylab='Density', main='', xlim=c(0,1.5), ylim=c(0,20), bty='n')
lines(seq(0,1.5,0.01), dexp(seq(0,1.5,0.01),0.001))
legend("topright", "G", bty="n", cex=3)
sh1 <- hist(op_str$scr[,1], plot=FALSE)
sh2 <- hist(op_str$scr[,2], plot=FALSE)
sh3 <- hist(op_str$scr[,3], plot=FALSE)
sh <- hist(op_ustr$scr[,1], plot=FALSE)
lines(rep(sh$breaks, each=2), c(0, rep(sh$density, each=2), 0), col='grey')
lines(rep(sh1$breaks, each=2), c(0, rep(sh1$density, each=2), 0), col='darkgreen')
lines(rep(sh2$breaks, each=2), c(0, rep(sh2$density, each=2), 0), col='blue')
lines(rep(sh3$breaks, each=2), c(0, rep(sh3$density, each=2), 0), col='red')

plot(0,0,pch='', xlab='Chlamydia prevalence (%)', ylab = 'Density', main='', xlim=c(0,15), ylim=c(0,0.7), bty='n')
legend("topright", "H", bty="n", cex=3)
ph1 <- hist(100*op_str$prev[,1], plot=FALSE)
ph2 <- hist(100*op_str$prev[,2], plot=FALSE)
ph3 <- hist(100*op_str$prev[,3], plot=FALSE)
ph <- hist(100*op_ustr$prev[,1], plot=FALSE)
lines(rep(ph$breaks, each=2), c(0, rep(ph$density, each=2), 0), col='grey')
lines(rep(ph1$breaks, each=2), c(0, rep(ph1$density, each=2), 0), col='darkgreen')
lines(rep(ph2$breaks, each=2), c(0, rep(ph2$density, each=2), 0), col='blue')
lines(rep(ph3$breaks, each=2), c(0, rep(ph3$density, each=2), 0), col='red')
arrows(2.2,0.7,4.3,0.7,angle=90,code=3,length=0.02, col='grey')
arrows(1.3,0.65,3.7,0.65,angle=90,code=3,length=0.02, col='darkgreen')
arrows(1.2,0.6,6.3,0.6,angle=90,code=3,length=0.02, col='blue')
arrows(3.5,0.55,9.8,0.55,angle=90,code=3,length=0.02, col='red')
points(3.1,0.7,pch=16, col='grey') # all women aged 16-24
points(2.2,0.65,pch=16, col='darkgreen') # women aged 16-24 reporting 0 new partners
points(2.8,0.6, pch=16, col='blue') # women aged 16-24 reporting 1 new partner
points(5.9,0.55, pch=16, col='red') # women aged 16-24 reporting 2+ new partners

# posterior summaries

t(apply(op0$foi, 2, quantile, p=c(0.5, 0.025, 0.975)))
t(apply(op0$scr, 2, quantile, p=c(0.5, 0.025, 0.975)))
t(apply(op0$prev, 2, quantile, p=c(0.5, 0.025, 0.975)))

# plot to show positivity

par(mfrow=c(1,1))

mli <- which(op0$lp__ == max(op0$lp__))[1]

trt <- op0$trt[mli] # treatment seeking rate in symptomatic people
sc <- op0$lambda_slow[mli] # rate of natural recovery (self-clear)
p_symp <- op0$p_symp[mli] # proportion of incident infections that develop symptoms

foi <- matrix(rep(seq(0,0.4,0.0005), times = 751), nrow=751, byrow=TRUE) # 801 values
scr <- matrix(rep(seq(0,1.5,0.002), each = 801), nrow=751, byrow=TRUE) # 751 values

alpha_AU <- scr + sc
alpha_UA <- foi * (1 - p_symp)
alpha_SU <- scr + sc + trt
alpha_US <- foi * p_symp

S <- alpha_AU*alpha_US/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))
A <- alpha_SU*alpha_UA/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))

tr <- scr + trt*S
dr <- scr*(S + A) + trt*S

positivity <- dr/tr

image(foi[1,], scr[,1], t(positivity), 
      col=gray(seq(1,0,-0.001)), 
      useRaster=TRUE, 
      xlab='Force of infection', ylab='Screening rate', 
      main = "Test positivity in women",
      xlim = c(0, 0.4), ylim=c(0, 1.5), zlim=c(0,1)
)

cs <- contourLines(foi[1,], scr[,1], t(positivity), 
                   levels=c(0.015, 0.050, 0.025, 0.091, 0.073, 0.156))

contour(foi[1,], scr[,1], t(positivity), 
        add=TRUE, 
        levels=c(0.015, 0.050, 0.025, 0.091, 0.073, 0.156), 
        lwd=3,
        col=rep(c('darkgreen','blue','red'),each=2),
        drawlabels=FALSE)

for(i in 1:6)
  text(cs[[i]][["x"]][1 + 40*ceiling(i/2)], 
       cs[[i]][["y"]][1 + 40*ceiling(i/2)], 
       cs[[i]][['level']], 
       pos=1, offset=1, col=rep(c('darkgreen','blue','red'),each=2)[i])

z1 <- ks::kde(matrix(c(op0$foi[,1], op0$scr[,1]), ncol=2))
z2 <- ks::kde(matrix(c(op0$foi[,2], op0$scr[,2]), ncol=2))
z3 <- ks::kde(matrix(c(op0$foi[,3], op0$scr[,3]), ncol=2))

contour(z1$eval.points[[1]], z1$eval.points[[2]], z1$estimate, 
        levels = z1$cont['5%'], 
        add=TRUE, drawlabels=FALSE, lty=2, lwd=3, col='darkgreen')
#text(0.05, 0.3, '0 new \npartners', adj=c(0,1))
contour(z2$eval.points[[1]], z2$eval.points[[2]], z2$estimate, 
        levels = z2$cont['5%'], 
        add=TRUE, drawlabels=FALSE, lty=2, lwd=3, col='blue')
#text(0.07, 0.45, '1 new \npartner', adj=c(0,1))
contour(z3$eval.points[[1]], z3$eval.points[[2]], z3$estimate, 
        levels = z3$cont['5%'], 
        add=TRUE, drawlabels=FALSE, lty=2, lwd=3, col='red')
#text(0.13, 0.65, '2+ new \npartners', adj=c(0,1))

legend('bottomright', inset = c(0.1,0.3),
       col = c('darkgreen', 'blue', 'red'),
       lwd = 3,
       legend = c('0 new partners', '1 new partner', '2+ new partners'),
       bty = 'n')

# box for comparison to men's plot
polygon(c(0, 0.2, 0.2, 0), c(0, 0, 0.8, 0.8))

