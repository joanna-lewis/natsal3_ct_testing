rm(list=ls())

library(rstan)
library(foreign)
library(survey)

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
natsal3$why2[natsal3$why %in% c("My partner had symptoms", "I was notified because a partner was diagnosed with Chlamydia")] <- "partner"
natsal3$why2[natsal3$why == "Check up after previous positive test"] <- "re-test"
natsal3$why2[natsal3$why %in% c("I wanted a general sexual health check-up", "I had no symptoms but I was worried about the risk of Chlamydia", "I was offered a routine test")] <- "screen"
natsal3$why2[natsal3$why == "Other"] <- "other"
natsal3$why2 <- factor(natsal3$why2, levels = names(table(natsal3$why2))[order(table(natsal3$why2), decreasing=TRUE)], ordered = TRUE)

sub <- subset(natsal3, natsal3$totnewy3 %in% c("0 new het.&/or sam. partners in last year","1 new het.&/or sam. partner in last year", "2+ new het.&/or sam. partners in last year"))
sub$totnewy3 <- factor(as.character(sub$totnewy3))

init0 <- list(p_symp = 0.5, lambda_slow = 0.74, foi = rep(0.001, times=3), scr = rep(0.001, times=3), trt = 0.001)

dt0 <- list(
  
  N = nrow(sub[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested),]),
  N_strata = length(unique(sub$totnewy3)),
  str = as.numeric(sub$totnewy3[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)]),
  wt = sub$total_wt[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)],
  tested = sub$tested[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)],
  symp = sub$why2[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)] == "symptoms",
#  symp = sub$why2[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)] %in% c("symptoms","partner","re-test","other"),
  diag = sub$whnchlam[sub$agrp == "16-24" & sub$rsex == "Female" & !is.na(sub$tested)] == "Less than 1 year ago"
  
)

dt0$symp[is.na(dt0$symp)] <- -99

fit0 <-stan(
  file = '3-comp-model-stratified-women.stan',
  data = dt0,
  chains = 1,
  iter = 10000,
  warmup = 5000,
  init = list(init0),
  seed = 12345
)

op0 <- extract(fit0)

quartz(height = 5.83, width= 8.27)
par(mfrow=c(2,3), mar=c(5,4,2,2))

h <- hist(100*op0$p_symp, plot=FALSE)
plot(rep(h$breaks, each=2), c(0, rep(h$density, each=2), 0), type='l', xlab = 'Proportion of infections symptomatic (%)', ylab='Density', main='', xlim=c(0,100), ylim=c(0,0.12), bty='n')
lines(100*seq(0,1,0.01), dbeta(seq(0,1,0.01),27, 90)/100)
legend("topright", "A", bty="n", cex=3, inset = c(0,-0.1))

h <- hist(op0$trt, plot=FALSE)
plot(rep(h$breaks, each=2), c(0, rep(h$density, each=2), 0), type='l',  xlab = expression("Symptomatic treatment rate (" ~ year^{-1} ~ ")"), ylab = 'Density', main='', xlim=c(0,35), ylim=c(0,0.11), bty='n')
lines(seq(0,100,0.1), dgamma(seq(0,100,0.1),14,1))
legend("topright", "B", bty="n", cex=3, inset = c(0,-0.1))

h <- hist(op0$lambda_slow, plot=FALSE)
#par(mgp=c(4,1,0))
plot(rep(h$breaks, each=2), c(0, rep(h$density, each=2), 0), type='l', xlab = expression('Natural clearance rate (' ~ year^{-1} ~ ')'), ylab='Denisty', main='', xlim=c(0,1.5), ylim=c(0,6), bty='n')
lines(seq(0,10,0.01), dnorm(seq(0,10,0.01),0.74, (0.89-0.61)/3.919928))
legend("topright", "C", bty="n", cex=3, inset = c(0,-0.1))
#par(mgp=c(3,1,0))
#mtext("Density", 2, cex=par('cex'), line=3)

plot(0,0,pch='',xlab = expression("Force of infection (" ~ year^{-1} ~ ")"), ylab='Density', main='', xlim=c(0,0.5), ylim=c(0,50), bty='n')
lines(seq(0,1,0.01), dexp(seq(0,1,0.01),0.001))
legend("topright", "D", bty="n", cex=3, inset = c(0,-0.1))
# stratified
fh1 <- hist(op0$foi[,1], plot=FALSE)
fh2 <- hist(op0$foi[,2], plot=FALSE)
fh3 <- hist(op0$foi[,3], plot=FALSE)
lines(rep(fh1$breaks, each=2), c(0, rep(fh1$density, each=2), 0), col='darkgreen')
lines(rep(fh2$breaks, each=2), c(0, rep(fh2$density, each=2), 0), col='blue')
lines(rep(fh3$breaks, each=2), c(0, rep(fh3$density, each=2), 0), col='red')
# unstratified
fh <- hist(op0$foi[,1], plot=FALSE)
lines(rep(fh$breaks, each=2), c(0, rep(fh$density, each=2), 0))

plot(0,0,pch='', xlab = expression("Screening rate (" ~ year^{-1} ~ ")"), ylab='Density', main='', xlim=c(0,1.5), ylim=c(0,20), bty='n')
lines(seq(0,1.5,0.01), dexp(seq(0,1.5,0.01),0.001))
legend("topright", "E", bty="n", cex=3, inset = c(0,-0.1))
# stratified
sh1 <- hist(op0$scr[,1], plot=FALSE)
sh2 <- hist(op0$scr[,2], plot=FALSE)
sh3 <- hist(op0$scr[,3], plot=FALSE)
lines(rep(sh1$breaks, each=2), c(0, rep(sh1$density, each=2), 0), col='darkgreen')
lines(rep(sh2$breaks, each=2), c(0, rep(sh2$density, each=2), 0), col='blue')
lines(rep(sh3$breaks, each=2), c(0, rep(sh3$density, each=2), 0), col='red')
# unstratified
sh <- hist(op0$scr[,1], plot=FALSE)
lines(rep(sh$breaks, each=2), c(0, rep(sh$density, each=2), 0))

plot(0,0,pch='', xlab='Chlamydia prevalence (%)', ylab = 'Density', main='', xlim=c(0,15), ylim=c(0,0.6), bty='n')
legend("topright", "F", bty="n", cex=3, inset = c(0,-0.1))
# stratified
ph1 <- hist(100*op0$prev[,1], plot=FALSE)
ph2 <- hist(100*op0$prev[,2], plot=FALSE)
ph3 <- hist(100*op0$prev[,3], plot=FALSE)
lines(rep(ph1$breaks, each=2), c(0, rep(ph1$density, each=2), 0), col='darkgreen')
lines(rep(ph2$breaks, each=2), c(0, rep(ph2$density, each=2), 0), col='blue')
lines(rep(ph3$breaks, each=2), c(0, rep(ph3$density, each=2), 0), col='red')
arrows(1.3,0.6,3.7,0.6,angle=90,code=3,length=0.02, col='darkgreen')
arrows(1.2,0.5,6.3,0.5,angle=90,code=3,length=0.02, col='blue')
arrows(3.5,0.4,9.8,0.4,angle=90,code=3,length=0.02, col='red')
points(2.2,0.6,pch=16, col='darkgreen') # women aged 16-24 reporting 0 new partners
points(2.8, 0.5, pch=16, col='blue') # women aged 16-24 reporting 1 new partner
points(5.9, 0.4, pch=16, col='red') # women aged 16-24 reporting 2+ new partners
# unstratified
ph <- hist(100*op0$prev[,1], plot=FALSE)
lines(rep(ph$breaks, each=2), c(0, rep(ph$density, each=2), 0))
arrows(2.2,0.55,4.3,0.55,angle=90,code=3,length=0.02)
points(3.1, 0.55, pch=16) # all women aged 16-24 

# posterior summaries

t(apply(op0$foi, 2, quantile, p=c(0.5, 0.025, 0.975)))
t(apply(op0$scr, 2, quantile, p=c(0.5, 0.025, 0.975)))
t(apply(op0$prev, 2, quantile, p=c(0.5, 0.025, 0.975)))

# plot to show positivity

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 4.1))

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
      xlab=expression("Force of infection (" ~ year^{-1} ~ ")"), 
      ylab=expression("Screening rate (" ~ year^{-1} ~ ")"), 
      main = "Test positivity in women",
      xlim = c(0, 0.4), ylim=c(0, 1.5), zlim=c(0,1)
      )

cs <- contourLines(foi[1,], scr[,1], t(positivity), 
                   levels=c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.5))
  
contour(foi[1,], scr[,1], t(positivity), 
        add=TRUE, 
        levels=c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.5), 
        lwd=1,
        drawlabels=FALSE)

for(i in 1:5)
  text(cs[[i]][["x"]][which(cs[[i]][["y"]] < 1.5)[1]], 
       cs[[i]][["y"]][which(cs[[i]][["y"]] < 1.5)[1]], 
       paste(100*cs[[i]][['level']], '%', sep=""), 
       adj=c(0,2)
       )
       
for(i in 6:7)
  text(cs[[i]][["x"]][which(cs[[i]][["x"]] < 0.4)[1]], 
       cs[[i]][["y"]][which(cs[[i]][["x"]] < 0.4)[1]], 
       paste(100*cs[[i]][['level']], '%', sep=""), 
       adj=c(-2,0), pos=2
       )

z1 <- ks::kde(matrix(c(op0$foi[,1], op0$scr[,1]), ncol=2))
z2 <- ks::kde(matrix(c(op0$foi[,2], op0$scr[,2]), ncol=2))
z3 <- ks::kde(matrix(c(op0$foi[,3], op0$scr[,3]), ncol=2))

contour(z1$eval.points[[1]], z1$eval.points[[2]], z1$estimate, 
        levels = z1$cont['5%'], 
        add=TRUE, drawlabels=FALSE, lwd=3, col='darkgreen')
#text(0.05, 0.3, '0 new \npartners', adj=c(0,1))
contour(z2$eval.points[[1]], z2$eval.points[[2]], z2$estimate, 
        levels = z2$cont['5%'], 
        add=TRUE, drawlabels=FALSE, lwd=3, col='blue')
#text(0.07, 0.45, '1 new \npartner', adj=c(0,1))
contour(z3$eval.points[[1]], z3$eval.points[[2]], z3$estimate, 
        levels = z3$cont['5%'], 
        add=TRUE, drawlabels=FALSE, lwd=3, col='red')
#text(0.13, 0.65, '2+ new \npartners', adj=c(0,1))

legend('bottomright', inset = c(0.05,0.15),
       col = c('darkgreen', 'blue', 'red'),
       lwd = 3,
       legend = c('0 new partners', '1 new partner', '≥2 new partners'),
       bty = 'n')

mtext('B', side=3, line=0.5, at=0.43, cex=5)

# box for comparison to men's plot
polygon(c(0, 0.2, 0.2, 0), c(0, 0, 0.8, 0.8))

