# Analysis of Natsal data to compare positivity among non-symptomatic testers with population prevalence.

rm(list=ls())

library(foreign)
library(survey)

natsal3 <- read.dta("/Users/Joanna/OneDrive - Imperial College London/backup/Natsal-3/UKDA-7799-stata11/stata11/eul_natsal_2010_for_archive.dta")

################################
# analysis on people who were tested for chlamydia
################################

# to be included in the subset tested for chlamydia in the last year, participants had to be:
# sexually active (hetlife > 0 or samlife > 0)
# aged < 45 (otherwise, were not asked whether they had tested)
# EITHER diagnosed within the last year, OR reported testing within the last year

tested <- natsal3[ (natsal3$whnchlam == "Less than 1 year ago") | (natsal3$chlmtest == "Yes"), ]
tested <- tested[ tested$chlmtest != "No", ]
tested <- tested[ !(tested$totlife %in% c(0,9999)), ]

##### restrict to 16-24-year-olds
#tested <- tested[tested$agrp == "16-24",]

##### reason recorded for chlamydia testing
tested$why <- factor(NA, levels = levels(tested$chtstwy1))
tested$why[!(tested$chtstwy1 %in% c("Not applicable", "Not answered"))] <- tested$chtstwy1[!(tested$chtstwy1 %in% c("Not applicable", "Not answered"))]
tested$why[!(tested$chltstwy %in% c("Not applicable", "Not answered"))] <- tested$chltstwy[!(tested$chltstwy %in% c("Not applicable", "Not answered"))]
tested$why <- droplevels(tested$why)

##### amalgamate two categories into "partner", and three into "screen"
tested$why2 <- NA
tested$why2[tested$why == "I had symptoms"] <- "symptoms"
tested$why2[tested$why %in% c("My partner had symptoms", "I was notified because a partner was diagnosed with Chlamydia")] <- "partner"
tested$why2[tested$why == "Check up after previous positive test"] <- "re-test"
tested$why2[tested$why %in% c("I wanted a general sexual health check-up", "I had no symptoms but I was worried about the risk of Chlamydia", "I was offered a routine test")] <- "screen"
tested$why2[tested$why == "Other"] <- "other"
tested$why2 <- factor(tested$why2, levels = names(table(tested$why2))[order(table(tested$why2), decreasing=TRUE)], ordered = TRUE)

##### distinguish between patient-initiated and provider-initiated screens
tested$why3 <- as.character(tested$why2)
tested$why3[tested$why %in% c("I wanted a general sexual health check-up", "I had no symptoms but I was worried about the risk of Chlamydia")] <- "screen - patient init"
tested$why3[tested$why == "I was offered a routine test"] <- "screen - provider init"
tested$why3 <- factor(tested$why3, levels = names(table(tested$why3))[order(table(tested$why3), decreasing=TRUE)], ordered = TRUE)

##### yet another categorisation!
tested$why4 <- NA
tested$why4[tested$why == "I had symptoms"] <- "symptoms"
tested$why4[tested$why %in% c("My partner had symptoms", "I wanted a general sexual health check-up","I had no symptoms but I was worried about the risk of Chlamydia", "I was offered a routine test", "Other", "Check up after previous positive test")] <- "screen"
tested$why4[tested$why %in% c("I was notified because a partner was diagnosed with Chlamydia")] <- "other"
tested$why4 <- factor(tested$why4, levels = names(table(tested$why4))[order(table(tested$why4), decreasing=TRUE)], ordered = TRUE)

##### group deprivation quintiles
tested$imd_grouped <- factor(tested$adj_imd_quintile, levels=1:5, labels=c("1-2","1-2","3","4-5","4-5"))

# set up survey design
tdesign <- svydesign( 
		id = ~psu_scrm , 
		strata = ~stratagrp3 ,
		data = tested ,		
		weight = ~total_wt 
	)

##################################################
# tabulate numbers of people reporting each reason
##################################################

table(tested$why, tested$agrp, tested$rsex)

##################################################
# positivity, by new partners in last year
##################################################

dt <- data.frame(
  X = 1:(3*4*2) + floor(0:(3*4*2 - 1) / 3), # three symptom groups; three risk groups, plus combined; two sexes
  Sex = rep(c("Male","Female"), each = 3*4),
  Age.group = rep("16-24", 3*4*2),
  New.partners.last.year = rep(c("Combined", "0", "1", "2+"), each = 3, times=2),
  Subset = c("All tested", "Symptoms", "Screens"),
  Estimate = NA,
  X2.5 = NA,
  X97.5 = NA
)

ind <- 1
for(i in c("Male","Female")){
	
	print(i)
	
  print("Combined")
  op <- NA
  try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                       design = subset(tdesign, rsex == i & agrp == "16-24")
  )))
  dt$Estimate[ind] <- 100*op[1]
  dt$X2.5[ind] <- 100*attr(op, "ci")[1]
  dt$X97.5[ind] <- 100*attr(op, "ci")[2]
  ind <- ind + 1
  
  op <- NA
  try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                       design = subset(tdesign, rsex == i & agrp == "16-24" & ((why == "I had symptoms" & whnchlam == "Less than 1 year ago") | why == "I was notified because a partner was diagnosed with Chlamydia"))
  )))
  dt$Estimate[ind] <- 100*op[1]
  dt$X2.5[ind] <- 100*attr(op, "ci")[1]
  dt$X97.5[ind] <- 100*attr(op, "ci")[2]
  ind <- ind + 1
  
  op <- NA
  try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                       design = subset(tdesign, rsex == i & agrp == "16-24" & !((why == "I had symptoms" & whnchlam == "Less than 1 year ago") | why == "I was notified because a partner was diagnosed with Chlamydia"))
  )))
  dt$Estimate[ind] <- 100*op[1]
  dt$X2.5[ind] <- 100*attr(op, "ci")[1]
  dt$X97.5[ind] <- 100*attr(op, "ci")[2]
  ind <- ind + 1
  
  for(k in c("0 new het.&/or sam. partners in last year", "1 new het.&/or sam. partner in last year", "2+ new het.&/or sam. partners in last year") ){
		
		print(k)
	  
	  op <- NA
		try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
			design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k)
			)))
	  dt$Estimate[ind] <- 100*op[1]
	  dt$X2.5[ind] <- 100*attr(op, "ci")[1]
	  dt$X97.5[ind] <- 100*attr(op, "ci")[2]
	  ind <- ind + 1
	  
	  op <- NA
		try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
			design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k & ((why == "I had symptoms" & whnchlam == "Less than 1 year ago") | why == "I was notified because a partner was diagnosed with Chlamydia"))
			)))
	  dt$Estimate[ind] <- 100*op[1]
	  dt$X2.5[ind] <- 100*attr(op, "ci")[1]
	  dt$X97.5[ind] <- 100*attr(op, "ci")[2]
	  ind <- ind + 1
	  
	  op <- NA
		try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
			design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k & !((why == "I had symptoms" & whnchlam == "Less than 1 year ago") | why == "I was notified because a partner was diagnosed with Chlamydia"))
			)))
	  dt$Estimate[ind] <- 100*op[1]
	  dt$X2.5[ind] <- 100*attr(op, "ci")[1]
	  dt$X97.5[ind] <- 100*attr(op, "ci")[2]
	  ind <- ind + 1
	  
	}
	
}



#dt <- read.csv("/Users/Joanna/OneDrive - Imperial College London/backup/ct_trends/Natsal_symptoms_testing/positivity_by_nonocong_1.csv")

quartz(width=18, height=6)
plot(dt$X, dt$Estimate, col = as.numeric(dt$Subset) + 1, 
	#pch = as.numeric(dt$Subset)-1,
	xaxt = 'n', ylim=c(0,25), xlab = "", ylab = "Percentage")

arrows(dt$X,
	y0 = dt$X2.5, y1 = dt$X97.5,
	col = as.numeric(dt$Subset) + 1,
	code = 3, angle = 90, length = 0.05
	)
	
mtext("0", side=1, line=1, at = 2, adj = 0.5)
mtext("1", side=1, line=1, at = 6, adj = 0.5)
mtext("2+", side=1, line=1, at = 10, adj = 0.5)

mtext("0", side=1, line=1, at = 14, adj = 0.5)
mtext("1", side=1, line=1, at = 18, adj = 0.5)
mtext("2+", side=1, line=1, at = 22, adj = 0.5)

mtext("Men", side=1, line=3, at = 6, adj = 0.5)
mtext("Women", side=1, line=3, at = 18, adj = 0.5)

axis(side = 1, line = 2.5, at = c(1,11), labels = FALSE, tcl=0)
axis(side = 1, line = 2.5, at = c(13,23), labels = FALSE, tcl=0)

legend('topleft', bty='n', inset = c(0.01,0.05),
	col = c(2:4, 'darkgrey'),
	pch = c(rep(1,3),15), 
	lty = 1,
#	fill = c(NA,NA,NA,rgb(0,0,0,0.1)),
#	x.intersp=c(2,2,2,0.5),
	border = NA,
  cex = 0.8, y.intersp = 1.3,
	legend = c("All tested", "Reported NOT symptoms as reason for test", "Reported screening as reason for test", "Natsal-3 prevalence estimates")
	)
	
# # lines for population prevalence
# lines(c(1, 14), 100*c(0.031, 0.031))
# lines(c(1, 14), 100*c(0.022, 0.022), lty=2)
# lines(c(1, 14), 100*c(0.043, 0.043), lty=2)
# 
# lines(c(16, 29), 100*c(0.023, 0.023))
# lines(c(16, 29), 100*c(0.015, 0.015), lty=2)
# lines(c(16, 29), 100*c(0.034, 0.034), lty=2)

# shading for population prevalence by group
lines(c(0.5, 3.5),c(1.8, 1.8))
polygon(c(0.5,3.5,3.5,0.5), c(0.9,0.9,3.8,3.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(4.5, 7.5), c(0.8, 0.8))
polygon(c(4.5,7.5,7.5,4.5), c(0.2,0.2,2.5,2.5), col=rgb(0,0,0,0.1), border=NA)

lines(c(8.5, 11.5),c(5.1, 5.1))
polygon(c(8.5,11.5,11.5,8.5), c(2.9,2.9,8.8,8.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(12.5, 15.5),c(2.2, 2.2))
polygon(c(12.5,15.5,15.5,12.5), c(1.3,1.3,3.7,3.7), col=rgb(0,0,0,0.1), border=NA)

lines(c(16.5, 19.5), c(2.8, 2.8))
polygon(c(16.5,19.5,19.5,16.5), c(1.2,1.2,6.3,6.3), col=rgb(0,0,0,0.1), border=NA)

lines(c(20.5, 23.5),c(5.9, 5.9))
polygon(c(20.5,23.5,23.5,20.5), c(3.5,3.5,9.8,9.8), col=rgb(0,0,0,0.1), border=NA)

#############################################
# All tests and Not symptoms only
#############################################

plot(1:(2*4*2) + floor(0:(2*4*2 - 1) / 2), 
     dt$Estimate[dt$Subset %in% c("All tested", "Screens")], 
     pch = c(0,5),
     xaxt = 'n', ylim=c(0,25), xlab = "", ylab = "Percentage positive")

arrows(1:(2*4*2) + floor(0:(2*4*2 - 1) / 2),
       y0 = dt$X2.5[dt$Subset %in% c("All tested", "Screens")], 
       y1 = dt$X97.5[dt$Subset %in% c("All tested", "Screens")],
       code = 3, angle = 90, length = 0.05
)

mtext("All", side=1, line=1, at = 1.5, adj = 0.5)
mtext("0", side=1, line=1, at = 4.5, adj = 0.5)
mtext("1", side=1, line=1, at = 7.5, adj = 0.5)
mtext("2+", side=1, line=1, at = 10.5, adj = 0.5)

mtext("All", side=1, line=1, at = 13.5, adj = 0.5)
mtext("0", side=1, line=1, at = 16.5, adj = 0.5)
mtext("1", side=1, line=1, at = 19.5, adj = 0.5)
mtext("2+", side=1, line=1, at = 22.5, adj = 0.5)

mtext("Men", side=1, line=3, at = 6, adj = 0.5)
mtext("Women", side=1, line=3, at = 18, adj = 0.5)

axis(side = 1, line = 2.5, at = c(1,11), labels = FALSE, tcl=0)
axis(side = 1, line = 2.5, at = c(13,23), labels = FALSE, tcl=0)

legend('topleft', bty='n', inset = c(0.01,0.05),
       col = c(1, 1, 'darkgrey'),
       pch = c(0,5,15), 
       lty = 1,
       #	fill = c(NA,NA,NA,rgb(0,0,0,0.1)),
       #	x.intersp=c(2,2,2,0.5),
       border = NA,
       cex = 0.8, y.intersp = 1.3,
       legend = c("All tested", "Reported screening as reason for test", "Natsal-3 prevalence estimates")
)

# shading for population prevalence by group
lines(c(0.5, 2.5),c(2.3, 2.3)) # all men
polygon(c(0.5,2.5,2.5,0.5), c(1.5,1.5,3.4,3.4), col=rgb(0,0,0,0.1), border=NA)

lines(c(3.5, 5.5),c(1.8, 1.8))
polygon(c(3.5,5.5,5.5,3.5), c(0.9,0.9,3.8,3.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(6.5, 8.5), c(0.8, 0.8))
polygon(c(6.5,8.5,8.5,6.5), c(0.2,0.2,2.5,2.5), col=rgb(0,0,0,0.1), border=NA)

lines(c(9.5, 11.5),c(5.1, 5.1))
polygon(c(9.5,11.5,11.5,9.5), c(2.9,2.9,8.8,8.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(12.5, 14.5),c(3.1, 3.1))
polygon(c(12.5,14.5,14.5,12.5), c(2.2,2.2,4.3,4.3), col=rgb(0,0,0,0.1), border=NA)

lines(c(15.5, 17.5),c(2.2, 2.2))
polygon(c(15.5,17.5,17.5,15.5), c(1.3,1.3,3.7,3.7), col=rgb(0,0,0,0.1), border=NA)

lines(c(18.5, 20.5), c(2.8, 2.8))
polygon(c(18.5,20.5,20.5,18.5), c(1.2,1.2,6.3,6.3), col=rgb(0,0,0,0.1), border=NA)

lines(c(21.5, 23.5),c(5.9, 5.9))
polygon(c(21.5,23.5,23.5,21.5), c(3.5,3.5,9.8,9.8), col=rgb(0,0,0,0.1), border=NA)


#############################################
# All Tests and Screens (why4 definition) only
#############################################

##################################################
# positivity, by new partners in last year
##################################################

dt_why <- data.frame(
  X = 1:(3*3*2) + floor(0:(3*3*2 - 1) / 3), # three symptom groups; three risk groups; two sexes
  Sex = rep(c("Male","Female"), each = 3*3),
  Age.group = rep("16-24", 3*3*2),
  New.partners.last.year = rep(c("0", "1", "2+"), each = 3, times=2),
  Subset = c("All tested", "Screens", "Not screens"),
  Estimate = NA,
  X2.5 = NA,
  X97.5 = NA
)

ind <- 1
for(i in c("Male","Female")){
  
  print(i)
  
  for(k in c("0 new het.&/or sam. partners in last year", "1 new het.&/or sam. partner in last year", "2+ new het.&/or sam. partners in last year") ){
    
    print(k)
    
    op <- NA
    try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                         design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k)
    )))
    dt_why$Estimate[ind] <- 100*op[1]
    dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
    dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
    ind <- ind + 1
    
    op <- NA
    try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                         design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k & !(why == "I had symptoms" & whnchlam == "Less than 1 year ago"))
    )))
    dt_why$Estimate[ind] <- 100*op[1]
    dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
    dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
    ind <- ind + 1
    
    op <- NA
    try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                         design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k & why == "I had symptoms" & whnchlam == "Less than 1 year ago")
    )))
    dt_why$Estimate[ind] <- 100*op[1]
    dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
    dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
    ind <- ind + 1
    
  }
  
}



plot(1:(2*3*2) + floor(0:(2*3*2 - 1) / 2), 
     dt_why$Estimate[dt_why$Subset %in% c("All tested", "Screens")], 
     pch = c(0,5),
     xaxt = 'n', ylim=c(0,25), xlab = "", ylab = "Percentage positive",
     main = "Screens -> NOT symptoms if diagnosed")

arrows(1:(2*3*2) + floor(0:(2*3*2 - 1) / 2),
       y0 = dt_why$X2.5[dt_why$Subset %in% c("All tested", "Screens")], 
       y1 = dt_why$X97.5[dt_why$Subset %in% c("All tested", "Screens")],
       code = 3, angle = 90, length = 0.05
)

mtext("0", side=1, line=1, at = 1.5, adj = 0.5)
mtext("1", side=1, line=1, at = 4.5, adj = 0.5)
mtext("2+", side=1, line=1, at = 7.5, adj = 0.5)

mtext("0", side=1, line=1, at = 10.5, adj = 0.5)
mtext("1", side=1, line=1, at = 13.5, adj = 0.5)
mtext("2+", side=1, line=1, at = 16.5, adj = 0.5)

mtext("Men", side=1, line=3, at = 4.5, adj = 0.5)
mtext("Women", side=1, line=3, at = 13.5, adj = 0.5)

axis(side = 1, line = 2.5, at = c(1,8), labels = FALSE, tcl=0)
axis(side = 1, line = 2.5, at = c(10,17), labels = FALSE, tcl=0)

legend('topleft', bty='n', inset = c(0.01,0.05),
       col = c(1, 1, 'darkgrey'),
       pch = c(0,5,15), 
       lty = 1,
       #	fill = c(NA,NA,NA,rgb(0,0,0,0.1)),
       #	x.intersp=c(2,2,2,0.5),
       border = NA,
       cex = 0.8, y.intersp = 1.3,
       legend = c("All tested", "Reported screening as reason for test", "Natsal-3 prevalence estimates")
)

# shading for population prevalence by group
lines(c(0.5, 2.5),c(1.8, 1.8))
polygon(c(0.5,2.5,2.5,0.5), c(0.9,0.9,3.8,3.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(3.5, 5.5), c(0.8, 0.8))
polygon(c(3.5,5.5,5.5,3.5), c(0.2,0.2,2.5,2.5), col=rgb(0,0,0,0.1), border=NA)

lines(c(6.5, 8.5),c(5.1, 5.1))
polygon(c(6.5,8.5,8.5,6.5), c(2.9,2.9,8.8,8.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(9.5, 11.5),c(2.2, 2.2))
polygon(c(9.5,11.5,11.5,9.5), c(1.3,1.3,3.7,3.7), col=rgb(0,0,0,0.1), border=NA)

lines(c(12.5, 14.5), c(2.8, 2.8))
polygon(c(12.5,14.5,14.5,12.5), c(1.2,1.2,6.3,6.3), col=rgb(0,0,0,0.1), border=NA)

lines(c(15.5, 17.5),c(5.9, 5.9))
polygon(c(15.5,17.5,17.5,15.5), c(3.5,3.5,9.8,9.8), col=rgb(0,0,0,0.1), border=NA)


##################################################
# three groups: all, exclude symptoms, screens only
##################################################

dt_why <- data.frame(
  X = 1:(4*3*2) + floor(0:(4*3*2 - 1) / 3), # three symptom groups; three risk groups; two sexes
  Sex = rep(c("Male","Female"), each = 4*3),
  Age.group = rep("16-24", 4*3*2),
  New.partners.last.year = rep(c("0", "1", "2+"), each = 4, times=2),
  Subset = c("All tested", "Not symptoms", "Screens"),
  Estimate = NA,
  X2.5 = NA,
  X97.5 = NA
)

ind <- 1
for(i in c("Male","Female")){
  
  print(i)
  
  print("Combined")
  op <- NA
  try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                       design = subset(tdesign, rsex == i & agrp == "16-24")
  )))
  dt_why$Estimate[ind] <- 100*op[1]
  dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
  dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
  ind <- ind + 1
  
  op <- NA
  try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                       design = subset(tdesign, rsex == i & agrp == "16-24" & !(why == "I had symptoms" & whnchlam == "Less than 1 year ago") )
  )))
  dt_why$Estimate[ind] <- 100*op[1]
  dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
  dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
  ind <- ind + 1
  
  op <- NA
  try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                       design = subset(tdesign, rsex == i & agrp == "16-24" & why %in% c("I wanted a general sexual health check-up", "I had no symptoms but I was worried about the risk of Chlamydia", "I was offered a routine test"))
  )))
  dt_why$Estimate[ind] <- 100*op[1]
  dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
  dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
  ind <- ind + 1
  
  
  for(k in c("0 new het.&/or sam. partners in last year", "1 new het.&/or sam. partner in last year", "2+ new het.&/or sam. partners in last year") ){
    
    print(k)
    
    op <- NA
    try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                         design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k)
    )))
    dt_why$Estimate[ind] <- 100*op[1]
    dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
    dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
    ind <- ind + 1
    
    op <- NA
    try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                         design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k & !(why == "I had symptoms" & whnchlam == "Less than 1 year ago") )
    )))
    dt_why$Estimate[ind] <- 100*op[1]
    dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
    dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
    ind <- ind + 1
    
    op <- NA
    try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                         design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k & why %in% c("I wanted a general sexual health check-up", "I had no symptoms but I was worried about the risk of Chlamydia", "I was offered a routine test")))))
    dt_why$Estimate[ind] <- 100*op[1]
    dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
    dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
    ind <- ind + 1
    
  }
  
}


# quartz(width=12, height=6)
plot(dt_why$X, 
     dt_why$Estimate, 
     pch = c(1,0,5),
     xaxt = 'n', ylim=c(0,25), xlab = "", ylab = "Percentage positive"
     )

arrows(dt_why$X,
       y0 = dt_why$X2.5, 
       y1 = dt_why$X97.5,
       code = 3, angle = 90, length = 0.05
)

mtext("All", side=1, line=1, at = 2, adj = 0.5)
mtext("0", side=1, line=1, at = 6, adj = 0.5)
mtext("1", side=1, line=1, at = 10, adj = 0.5)
mtext("≥2", side=1, line=1, at = 14, adj = 0.5)

mtext("All", side=1, line=1, at = 18, adj = 0.5)
mtext("0", side=1, line=1, at = 22, adj = 0.5)
mtext("1", side=1, line=1, at = 26, adj = 0.5)
mtext("≥2", side=1, line=1, at = 30, adj = 0.5)

mtext("Men", side=1, line=3, at = 8, adj = 0.5)
mtext("Women", side=1, line=3, at = 24, adj = 0.5)

axis(side = 1, line = 2.5, at = c(1,15), labels = FALSE, tcl=0)
axis(side = 1, line = 2.5, at = c(17,31), labels = FALSE, tcl=0)

# legend('topleft', bty='n', inset = c(0.01,0.05),
       # col = c(1, 1, 1, 'darkgrey'),
       # pch = c(1,0,5,15), 
       # lty = 1,
       # #	fill = c(NA,NA,NA,rgb(0,0,0,0.1)),
       # #	x.intersp=c(2,2,2,0.5),
       # border = NA,
       # cex = 0.8, y.intersp = 1.3,
       # legend = c("All tested", "Exclude symptomatic diagnoses", "Screens only", "Natsal-3 prevalence estimates")
# )

# legend
points(c(1,1,1), c(24,22,20), pch=c(1,0,5))
lines(c(0.5,1.5),c(18,18))
polygon(c(0.5, 1.5, 1.5, 0.5), c(17, 17, 19, 19), col=rgb(0,0,0,0.1), border=NA)
text(c(1.5,1.5,1.5,1.5), c(24,22,20,18), pos=4, labels = c("All tested", "Exclude symptomatic diagnoses", "Screens only", "Natsal-3 prevalence estimates"))

# shading for population prevalence by group
lines(c(0.5, 3.5),c(2.3, 2.3))
polygon(c(0.5,3.5,3.5,0.5), c(1.5,1.5,3.4,3.4), col=rgb(0,0,0,0.1), border=NA)

lines(c(4.5, 7.5),c(1.8, 1.8))
polygon(c(4.5,7.5,7.5,4.5), c(0.9,0.9,3.8,3.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(8.5, 11.5), c(0.8, 0.8))
polygon(c(8.5,11.5,11.5,8.5), c(0.2,0.2,2.5,2.5), col=rgb(0,0,0,0.1), border=NA)

lines(c(12.5, 15.5),c(5.1, 5.1))
polygon(c(12.5,15.5,15.5,12.5), c(2.9,2.9,8.8,8.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(16.5, 19.5),c(3.1, 3.1))
polygon(c(16.5,19.5,19.5,16.5), c(2.2,2.2,4.3,4.3), col=rgb(0,0,0,0.1), border=NA)

lines(c(20.5, 23.5),c(2.2, 2.2))
polygon(c(20.5,23.5,23.5,20.5), c(1.3,1.3,3.7,3.7), col=rgb(0,0,0,0.1), border=NA)

lines(c(24.5, 27.5), c(2.8, 2.8))
polygon(c(24.5,27.5,27.5,24.5), c(1.2,1.2,6.3,6.3), col=rgb(0,0,0,0.1), border=NA)

lines(c(28.5, 31.5),c(5.9, 5.9))
polygon(c(28.5,31.5,31.5,28.5), c(3.5,3.5,9.8,9.8), col=rgb(0,0,0,0.1), border=NA)

##################################################
# three groups: all, exclude symptoms and PN, screens only
##################################################

dt_why <- data.frame(
  X = 1:(4*3*2) + floor(0:(4*3*2 - 1) / 3), # three symptom groups; three risk groups; two sexes
  Sex = rep(c("Male","Female"), each = 4*3),
  Age.group = rep("16-24", 4*3*2),
  New.partners.last.year = rep(c("0", "1", "2+"), each = 4, times=2),
  Subset = c("All tested", "Not symptoms", "Screens"),
  Estimate = NA,
  X2.5 = NA,
  X97.5 = NA
)

ind <- 1
for(i in c("Male","Female")){
  
  print(i)
  
  print("Combined")
  op <- NA
  try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                       design = subset(tdesign, rsex == i & agrp == "16-24")
  )))
  dt_why$Estimate[ind] <- 100*op[1]
  dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
  dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
  ind <- ind + 1
  
  op <- NA
  try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                       design = subset(tdesign, rsex == i & agrp == "16-24" & !(why == "I had symptoms" & whnchlam == "Less than 1 year ago") & (why != "I was notified because a partner was diagnosed with Chlamydia") )
  )))
  dt_why$Estimate[ind] <- 100*op[1]
  dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
  dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
  ind <- ind + 1
  
  op <- NA
  try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                       design = subset(tdesign, rsex == i & agrp == "16-24" & why %in% c("I wanted a general sexual health check-up", "I had no symptoms but I was worried about the risk of Chlamydia", "I was offered a routine test"))
  )))
  dt_why$Estimate[ind] <- 100*op[1]
  dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
  dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
  ind <- ind + 1
  
  
  for(k in c("0 new het.&/or sam. partners in last year", "1 new het.&/or sam. partner in last year", "2+ new het.&/or sam. partners in last year") ){
    
    print(k)
    
    op <- NA
    try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                         design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k)
    )))
    dt_why$Estimate[ind] <- 100*op[1]
    dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
    dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
    ind <- ind + 1
    
    op <- NA
    try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                         design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k & !(why == "I had symptoms" & whnchlam == "Less than 1 year ago") & (why != "I was notified because a partner was diagnosed with Chlamydia") )
    )))
    dt_why$Estimate[ind] <- 100*op[1]
    dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
    dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
    ind <- ind + 1
    
    op <- NA
    try(op <- (svyciprop(~(whnchlam == "Less than 1 year ago"),
                         design = subset(tdesign, rsex == i & agrp == "16-24" & totnewy3 == k & why %in% c("I wanted a general sexual health check-up", "I had no symptoms but I was worried about the risk of Chlamydia", "I was offered a routine test")))))
    dt_why$Estimate[ind] <- 100*op[1]
    dt_why$X2.5[ind] <- 100*attr(op, "ci")[1]
    dt_why$X97.5[ind] <- 100*attr(op, "ci")[2]
    ind <- ind + 1
    
  }
  
}


quartz(width=12, height=6)
plot(dt_why$X, 
     dt_why$Estimate, 
     pch = c(1,0,5),
     xaxt = 'n', ylim=c(0,25), xlab = "", ylab = "Percentage positive"
     )

arrows(dt_why$X,
       y0 = dt_why$X2.5, 
       y1 = dt_why$X97.5,
       code = 3, angle = 90, length = 0.05
)

mtext("All", side=1, line=1, at = 2, adj = 0.5)
mtext("0", side=1, line=1, at = 6, adj = 0.5)
mtext("1", side=1, line=1, at = 10, adj = 0.5)
mtext("≥2", side=1, line=1, at = 14, adj = 0.5)

mtext("All", side=1, line=1, at = 18, adj = 0.5)
mtext("0", side=1, line=1, at = 22, adj = 0.5)
mtext("1", side=1, line=1, at = 26, adj = 0.5)
mtext("≥2", side=1, line=1, at = 30, adj = 0.5)

mtext("Men", side=1, line=3, at = 8, adj = 0.5)
mtext("Women", side=1, line=3, at = 24, adj = 0.5)

axis(side = 1, line = 2.5, at = c(1,15), labels = FALSE, tcl=0)
axis(side = 1, line = 2.5, at = c(17,31), labels = FALSE, tcl=0)

# legend('topleft', bty='n', inset = c(0.01,0.05),
       # col = c(1, 1, 1, 'darkgrey'),
       # pch = c(1,0,5,15), 
       # lty = 1,
       # #	fill = c(NA,NA,NA,rgb(0,0,0,0.1)),
       # #	x.intersp=c(2,2,2,0.5),
       # border = NA,
       # cex = 0.8, y.intersp = 1.3,
       # legend = c("All tested", "Exclude symptomatic diagnoses", "Screens only", "Natsal-3 prevalence estimates")
# )

# legend
points(c(1,1,1), c(24,22,20), pch=c(1,0,5))
lines(c(0.5,1.5),c(18,18))
polygon(c(0.5, 1.5, 1.5, 0.5), c(17, 17, 19, 19), col=rgb(0,0,0,0.1), border=NA)
text(c(1.5,1.5,1.5,1.5), c(24,22,20,18), pos=4, labels = c("All tested", "Exclude symptomatic diagnoses and partner notifications", "Screens only", "Natsal-3 prevalence estimates"))

# shading for population prevalence by group
lines(c(0.5, 3.5),c(2.3, 2.3))
polygon(c(0.5,3.5,3.5,0.5), c(1.5,1.5,3.4,3.4), col=rgb(0,0,0,0.1), border=NA)

lines(c(4.5, 7.5),c(1.8, 1.8))
polygon(c(4.5,7.5,7.5,4.5), c(0.9,0.9,3.8,3.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(8.5, 11.5), c(0.8, 0.8))
polygon(c(8.5,11.5,11.5,8.5), c(0.2,0.2,2.5,2.5), col=rgb(0,0,0,0.1), border=NA)

lines(c(12.5, 15.5),c(5.1, 5.1))
polygon(c(12.5,15.5,15.5,12.5), c(2.9,2.9,8.8,8.8), col=rgb(0,0,0,0.1), border=NA)

lines(c(16.5, 19.5),c(3.1, 3.1))
polygon(c(16.5,19.5,19.5,16.5), c(2.2,2.2,4.3,4.3), col=rgb(0,0,0,0.1), border=NA)

lines(c(20.5, 23.5),c(2.2, 2.2))
polygon(c(20.5,23.5,23.5,20.5), c(1.3,1.3,3.7,3.7), col=rgb(0,0,0,0.1), border=NA)

lines(c(24.5, 27.5), c(2.8, 2.8))
polygon(c(24.5,27.5,27.5,24.5), c(1.2,1.2,6.3,6.3), col=rgb(0,0,0,0.1), border=NA)

lines(c(28.5, 31.5),c(5.9, 5.9))
polygon(c(28.5,31.5,31.5,28.5), c(3.5,3.5,9.8,9.8), col=rgb(0,0,0,0.1), border=NA)

