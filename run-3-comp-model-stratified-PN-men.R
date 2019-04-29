rm(list=ls())

library(rstan)
library(foreign)
library(survey)
library(MASS)

####################
# data
####################

setwd("~/OneDrive - Imperial College London/backup/papers/reasons_and_venues/natsal3_ct_testing")

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

####################
# run Stan model
####################

##########
# with stratification
##########

dt_str <- list(
  
  N = nrow(sub[sub$agrp == "16-24" & sub$rsex == "Male" & !is.na(sub$tested),]),
  N_strata = length(unique(sub$totnewy3)),
  str = (as.numeric(sub$totnewy3[sub$agrp == "16-24" & sub$rsex == "Male" & !is.na(sub$tested)])),
  wt = sub$total_wt[sub$agrp == "16-24" & sub$rsex == "Male" & !is.na(sub$tested)],
  tested = sub$tested[sub$agrp == "16-24" & sub$rsex == "Male" & !is.na(sub$tested)],
  symp = sub$why2[sub$agrp == "16-24" & sub$rsex == "Male" & !is.na(sub$tested)] == "symptoms",
  partner = sub$why2[sub$agrp == "16-24" & sub$rsex == "Male" & !is.na(sub$tested)] == "partner",
  diag = sub$whnchlam[sub$agrp == "16-24" & sub$rsex == "Male" & !is.na(sub$tested)] == "Less than 1 year ago"
  
)

dt_str$symp[is.na(dt_str$symp)] <- -99
dt_str$partner[is.na(dt_str$partner)] <- -99

##########
# without stratification
##########

dt_ustr <- dt_str
dt_ustr$str <- 1 + 0*(as.numeric(sub$totnewy3[sub$agrp == "16-24" & sub$rsex == "Male" & !is.na(sub$tested)]))

fit_str <-stan(
  file = '3-comp-model-stratified-PN-men.stan',
  data = dt_str,
  chains = 1,
  iter = 10000,
  warmup = 1000,
  init = list(init0),
  seed = 12345
)
op_str <- extract(fit_str)

fit_ustr <-stan(
  file = '3-comp-model-stratified-PN-men.stan',
  data = dt_ustr,
  chains = 1,
  iter = 10000,
  warmup = 1000,
  init = list(init0),
  seed = 12345
)
op_ustr <- extract(fit_ustr)

####################
# plot posteriors
####################

quartz(height = 9, width= 9 )
par(mfrow=c(3,3), mar=c(5,4,2,2))

h_s <- hist(100*op_str$p_symp, plot=FALSE)
h_u <- hist(100*op_ustr$p_symp, plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = 'Proportion of infections symptomatic (%)', ylab='Density', main='', xlim=c(0,100), ylim=c(0,0.12), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(100*seq(0,1,0.01), dbeta(seq(0,1,0.01), 11, 5)/100)
legend("topright", "A", bty="n", cex=2)

h_s <- hist(op_str$trt, plot=FALSE)
h_u <- hist(op_ustr$trt, plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = expression("Symptomatic treatment rate (" ~ year^{-1} ~ ")"), ylab = 'Density', main='', xlim=c(0,35), ylim=c(0,0.11), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(seq(0,100,0.1), dgamma(seq(0,100,0.1),14,1))
legend("topright", "B", bty="n", cex=2)

h_s <- hist(op_str$lambda_slow, plot=FALSE)
h_u <- hist(op_ustr$lambda_slow, plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = expression('Natural clearance rate (' ~ year^{-1} ~ ')'), ylab='Denisty', main='', xlim=c(0,1.5), ylim=c(0,6), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(seq(0,10,0.01), dlnorm(seq(0,10,0.01),log(0.42), 0.4))
legend("topright", "C", bty="n", cex=2)

h_s <- hist(op_str$pn[,1], plot=FALSE)
h_u <- hist(op_ustr$pn[,1], plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = expression('Partner notification rate (uninfected men) (' ~ year^{-1} ~ ')'), ylab='Denisty', main='', xlim=c(0,0.025), ylim=c(0,200), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(seq(0,1,0.001), dexp(seq(0,1,0.001),0.001))
legend("topright", "D", bty="n", cex=2)

h_s <- hist(op_str$pn[,2], plot=FALSE)
h_u <- hist(op_ustr$pn[,2], plot=FALSE)
plot(rep(h_u$breaks, each=2), c(0, rep(h_u$density, each=2), 0), type='l', col='grey', xlab = expression('Partner notification rate (infected men) (' ~ year^{-1} ~ ')'), ylab='Denisty', main='', xlim=c(0,1), ylim=c(0,10), bty='n')
lines(rep(h_s$breaks, each=2), c(0, rep(h_s$density, each=2), 0))
lines(seq(0,1,0.001), dexp(seq(0,1,0.001),0.001))
legend("topright", "E", bty="n", cex=2)

plot.new()

plot(0,0,pch='',xlab = expression("Force of infection (" ~ year^{-1} ~ ")"), ylab='Density', main='', xlim=c(0,0.5), ylim=c(0,50), bty='n')
lines(seq(0,1,0.01), dexp(seq(0,1,0.01),0.001))
legend("topright", "F", bty="n", cex=2)
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
arrows(1.5,0.7,3.4,0.7,angle=90,code=3,length=0.02, col='grey')
arrows(0.9,0.65,3.8,0.65,angle=90,code=3,length=0.02, col='darkgreen')
arrows(0.2,0.6,2.5,0.6,angle=90,code=3,length=0.02, col='blue')
arrows(2.9,0.55,8.8,0.55,angle=90,code=3,length=0.02, col='red')
points(2.3, 0.7, pch=16, col='grey') # all men aged 16-24 
points(1.8,0.65,pch=16, col='darkgreen') # men aged 16-24 reporting 0 new partners
points(0.8, 0.6, pch=16, col='blue') # men aged 16-24 reporting 1 new partner
points(5.1, 0.55, pch=16, col='red') # men aged 16-24 reporting 2+ new partners

##########
# posterior summaries
##########

t(apply(op0$foi, 2, quantile, p=c(0.5, 0.025, 0.975)))
t(apply(op0$scr, 2, quantile, p=c(0.5, 0.025, 0.975)))
t(apply(op0$prev, 2, quantile, p=c(0.5, 0.025, 0.975)))
