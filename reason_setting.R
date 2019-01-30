# Analysis of Natsal data to investigate correlations between reason for test and test setting

rm(list=ls())

library(foreign)
library(survey)
library(vcd)

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

##### symptoms/screen/other
tested$why4 <- as.character(tested$why2)
tested$why4[tested$why2 %in% c("partner", "re-test", "other")] <- "other"
tested$why4 <- factor(tested$why4, levels = c("other", "symptoms", "screen"), ordered = TRUE)

##### where were tests carried out?
tested$where <- factor(NA, levels = unique(c(levels(tested$chtstwh1), levels(tested$chltstwh))))
tested$where[!(tested$chtstwh1 %in% c("Not applicable", "Not answered"))] <- tested$chtstwh1[!(tested$chtstwh1 %in% c("Not applicable", "Not answered"))]
tested$where[!(tested$chltstwh %in% c("Not applicable", "Not answered"))] <- tested$chltstwh[!(tested$chltstwh %in% c("Not applicable", "Not answered"))]
tested$where[tested$where == "NHS Family planning clinic / contraceptive clinic / reproductive health clinic"] <- "NHS FP clinic / contraceptive clinic / reproductive health clinic"
tested$where <- droplevels(tested$where)

##### combine categories
tested$where2 <- NA
tested$where2[tested$where == "Sexual health clinic (GUM clinic)"] <- "sexual health"
tested$where2[tested$where == "General practice (GP) surgery"] <- "GP"
tested$where2[tested$where == "NHS FP clinic / contraceptive clinic / reproductive health clinic"] <- "family planning"
tested$where2[tested$where == "School / college / university"] <- "education"
tested$where2[tested$where %in% c("Ante-natal Clinic / midwife", "Private non-NHS clinics or doctor", "Youth advisory clinic (e.g. Brook Clinic)", "Termination of pregnancy (abortion) clinic", "Hospital accident and emergency (A&E) department", "Pharmacy / chemist", "Internet", "Other non-health care place, e.g. youth club, festival, bar", "Somewhere else")] <- "elsewhere"

tested$where2 <- factor(tested$where2, levels = names(table(tested$where2))[order(table(tested$where2), decreasing=TRUE)], ordered = TRUE)

##### only three categories
tested$where3 <- tested$where2
tested$where3[tested$where2 %in% c("family planning", "education", "elsewhere")] <- "elsewhere"
tested$where3 <- factor(tested$where3, levels = c("sexual health", "GP", "elsewhere"), ordered = TRUE)


##### group deprivation quintiles
tested$imd_grouped <- factor(tested$adj_imd_quintile, levels=1:5, labels=c("1-2","1-2","3","4-5","4-5"))

# set up survey design
tdesign <- svydesign( 
		id = ~psu_scrm , 
		strata = ~stratagrp3 ,
		data = tested ,		
		weight = ~total_wt 
	)

##################
# mosaic plots
##################

cols <- c(rgb(135,93,167, maxColorValue = 256), rgb(168,91,74, maxColorValue = 256), rgb(143,183,113, maxColorValue = 256))

par(mfrow=c(1,2))

venue_reason_table_m <- svytable(~ where3 + why4, design=subset(tdesign, rsex == "Male"), Ntotal = 1 )
mosaicplot(venue_reason_table_m, color=cols, xlab="Test setting", ylab = "Reason for test", main= "Men, 16-24 years")

venue_reason_table_f <- svytable(~ where3 + why4, design=subset(tdesign, rsex == "Female"), Ntotal = 1  )
mosaicplot(venue_reason_table_f, color=cols, xlab="Test setting", ylab = "Reason for test", main= "Women, 16-24 years")

##################
# some proportions
##################
svytable(~ where3, design=subset(tdesign, rsex == "Male"), Ntotal = 1 )
svytable(~ where3, design=subset(tdesign, rsex == "Female"), Ntotal = 1 )

svytable(~ why4, design=subset(tdesign, rsex == "Male"), Ntotal = 1 )
svytable(~ why4, design=subset(tdesign, rsex == "Female"), Ntotal = 1 )

svytable(~ where3, design=subset(tdesign, rsex == "Male" & why4 == "symptoms"), Ntotal = 1 )
svytable(~ where3, design=subset(tdesign, rsex == "Female" & why4 == "symptoms"), Ntotal = 1 )


##################
# statistical tests for association
##################
svychisq(~ where3 + why4, design=subset(tdesign, rsex == "Male") )
svychisq(~ where3 + why4, design=subset(tdesign, rsex == "Female") )

