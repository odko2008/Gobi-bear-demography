#formatting Gobi bear data for captures from 2009, 2013, and 2017 
#install_github("jaroyle/oSCR")
library(oSCR)

load("Gobi_Ch.RData")
load("Gobi_traps.RData")

head(CR.data)
unique(CR.data$TrapID)
head(traps.tmp)

plot(traps.tmp$posx,traps.tmp$posy)

#################
## Keep all from years 2009, 2013, 2017
Sub.years.all = subset(CR.data, Year == "2009" | Year == "2013" | Year == "2017") # same as CR.data
str(Sub.years.all)

length(unique(table(Sub.years.all$AnimalID, Sub.years.all$Year)))

edf.data<-data.frame(session=as.integer(as.factor(Sub.years.all$Year)), individual=as.factor(Sub.years.all$AnimalID),
                     occasion=as.integer(Sub.years.all$Mod_Ses),TrapID=as.factor(Sub.years.all$TrapID),sex=as.factor(Sub.years.all$Sex)) 
str(edf.data)
## to process all sessions for multi-strata model change up session, but you'll want it to be a vector of integers 1:session to work best 
traps<-data.frame(TrapID=as.factor(traps.tmp[order(traps.tmp$TrapID),]$TrapID), 
                  Xcoord=as.integer(traps.tmp[order(traps.tmp$TrapID),]$posx/1000),
                  Ycoord=as.integer(traps.tmp[order(traps.tmp$TrapID),]$posy/1000))
traps
tdf.data<-list(traps,traps,traps) ### for each year for trap list 
str(tdf.data)
head(edf.data)
levels(edf.data$TrapID)

head(tdf.data)

scr.data<-data2oscr(edf=edf.data,sess.col=1,id.col=2,occ.col=3,trap.col=4,sex.col=5,tdf=tdf.data, K = c(4,5,6), ntraps=c(13,13,13),
                    remove.extracaps = FALSE) #counts for Poisson
sf<-scr.data$scrFrame
str(sf)
plot(sf)

## second version of SSDF (after testing for sensitivity of density estimates to state space size)
plot(ssDF1)
points(tdf.data[[1]]$Xcoord,tdf.data[[1]]$Ycoord, pch=19, col="red") # putting trap location info on the ssDF figure 
out.formatted1<-list(sf=sf,ssDF1=ssDF1,res=0.1*sf$mmdm)
data1<-out.formatted1

##### Running 18 SCR candidate models ###### 
## 1. Null Model ##
mod.oSCR1<-oSCR.fit(list(D~1, p0~1, sig~1),data1$sf,data1$ssDF1, encmod="P")
mod.oSCR1
#density
exp(mod.oSCR1$coef.mle$mle[3])*nrow(ssDF1[[1]])

## 2. Sex specific baseline detection## 
m2a<-oSCR.fit(list(D~1, p0~sex, sig~1),data1$sf,data1$ssDF1, encmod="P")
m2a # m2a and m2 are same. It means the additional information with open days does not influence to the model runs
exp(m2a$coef.mle$mle[4])*nrow(ssDF1[[1]]) 

## 3. sex specific spatial scalar ## 
m3a<-oSCR.fit(list(D~1, p0~1, sig~sex),data1$sf,data1$ssDF1, encmod="P")
m3a 
exp(m3a$coef.mle$mle[4])*nrow(ssDF1[[1]]) 

## 4. Session specific density ## 
m4a<-oSCR.fit(list(D~session, p0~1, sig~1),data1$sf,data1$ssDF1, encmod="P")
m4a
exp(m4a$coef.mle$mle[3])*nrow(ssDF1[[1]]) 

## 5. Session specific baseline detection ## 
m5a<-oSCR.fit(list(D~1, p0~session, sig~1),data1$sf,data1$ssDF1, encmod="P")
m5a 
exp(m5a$coef.mle$mle[5])*nrow(ssDF1[[1]]) ## 

## 6. Session specific spatial scalar ## 
m6a<-oSCR.fit(list(D~1, p0~1, sig~session),data1$sf,data1$ssDF1, encmod="P")
m6a 
exp(m6a$coef.mle$mle[5])*nrow(ssDF1[[1]]) 

##7. Session specifc eveything
m7a<-oSCR.fit(list(D~session, p0~session, sig~session),data1$sf,data1$ssDF1, encmod="P")
m7a 

##8. session specific density and encounter rate, sex specifc sigma
m8a<-oSCR.fit(list(D~session, p0~session, sig~sex),data1$sf,data1$ssDF1, encmod="P")
m8a 

## 9. Sex specific baseline detection and Session specific spatial scalar ## 
m9a<-oSCR.fit(list(D~1, p0~sex, sig~session),data1$sf,data1$ssDF1, encmod="P")
m9a
exp(m9a$coef.mle$mle[6])*nrow(ssDF1[[1]]) ## 

## 10. Session specific baseline detection and sex specific spatial scalar ## 
m10a<-oSCR.fit(list(D~1, p0~session, sig~sex),data1$sf,data1$ssDF1, encmod="P")
m10a 
exp(m10a$coef.mle$mle[6])*nrow(ssDF1[[1]])

## 11. Session specific density and sex specific baseline detection ##  
m11a<-oSCR.fit(list(D~session, p0~sex, sig~1),data1$sf,data1$ssDF1, encmod="P")
m11a
exp(m11a$coef.mle$mle[4])*nrow(ssDF1[[1]]) # 

##12. Session specific density and sex specific spatial detection ##  
m12a<-oSCR.fit(list(D~session, p0~1, sig~sex),data1$sf,data1$ssDF1, encmod="P")
m12a
exp(m12a$coef.mle$mle[4])*nrow(ssDF1[[1]])

##13. Session specifc density, sex specifc encounter rate and sigma
m13a<-oSCR.fit(list(D~session, p0~sex, sig~sex),data1$sf,data1$ssDF1, encmod="P")
m13a 
exp(m13a$coef.mle$mle[5])*nrow(ssDF1[[1]])

## 14. 
m14a<-oSCR.fit(list(D~1, p0~sex+session, sig~1),data1$sf,data1$ssDF1, encmod="P")
m14a
exp(m14a$coef.mle$mle[6])*nrow(ssDF1[[1]])

##15. session specifc density and sigma, sex specifc encounter rate
m15a<-oSCR.fit(list(D~session, p0~sex, sig~session),data1$sf,data1$ssDF1, encmod="P")
m15a 

## 16. 
m16a<-oSCR.fit(list(D~1, p0~sex+session, sig~session),data1$sf,data1$ssDF1, encmod="P")  ###
m16a
exp(m16a$coef.mle$mle[8])*nrow(ssDF1[[1]]) # 27 individuals

# 17. 
m17a<-oSCR.fit(list(D~session, p0~session, sig~1),data1$sf,data1$ssDF1, encmod="P")
m17a

#18.
m18a<-oSCR.fit(list(D~session, p0~1, sig~session),data1$sf,data1$ssDF1, encmod="P") ##
exp(m18$coef.mle$mle[5])*nrow(ssDF1[[1]]) # d0.(Intercept) -4.9649136

#19.
m19a<-oSCR.fit(list(D~1, p0~session, sig~session),data1$sf,data1$ssDF1, encmod="P")
exp(m19$coef.mle$mle[7])*nrow(ssDF1[[1]])

#20.
m20a<-oSCR.fit(list(D~session, p0~1, sig~session+sex),data1$sf,data1$ssDF1, encmod="P")

#21 session specifc density, session and sess specifc encounter rate, constant spatail scalar
m21a<-oSCR.fit(list(D~session, p0~session + sex, sig~1),data1$sf,data1$ssDF1, encmod="P")
m21a 

#22.
m22a<-oSCR.fit(list(D~1, p0~sex, sig~session+sex),data1$sf,data1$ssDF1, encmod="P")
exp(m22$coef.mle$mle[7])*nrow(ssDF1[[1]])

#23 session specifc density and spatial scalar, session and sess specifc encounter rate
m23a<-oSCR.fit(list(D~session, p0~session + sex, sig~session),data1$sf,data1$ssDF1, encmod="P")
m23a 

#24 session specifc density, session and sess specifc encounter rate, sex specifc spatial scalar
m24a<-oSCR.fit(list(D~session, p0~session + sex, sig~sex),data1$sf,data1$ssDF1, encmod="P")
m24a 

#25 density and encounter rate constant, sex and session specifc spatial scalar
m25a<-oSCR.fit(list(D~1, p0~1, sig~session + sex),data1$sf,data1$ssDF1, encmod="P")
m25a 

#26 session specific density and encounter rate, session and sex dependent sigma
m26a<-oSCR.fit(list(D~session, p0~session, sig~session + sex),data1$sf,data1$ssDF1, encmod="P")
m26a 

#27. 
m27a<-oSCR.fit(list(D~1, p0~sex+session, sig~sex),data1$sf,data1$ssDF1, encmod="P")
exp(m27a$coef.mle$mle[7])*nrow(ssDF1[[1]])

#28.
m28a<-oSCR.fit(list(D~1, p0~sex, sig~sex),data1$sf,data1$ssDF1, encmod="P")
exp(m28a$coef.mle$mle[7])*nrow(ssDF1[[1]])

#29 constant density, session specific encounter rate, session and sex dependent sigma
m29a<-oSCR.fit(list(D~session, p0~session, sig~session + sex),data1$sf,data1$ssDF1, encmod="P")
m29a 

#30 session specific density, sex and session specific encounter rate and sigma
m30a<-oSCR.fit(list(D~session, p0~session + sex, sig~session + sex),data1$sf,data1$ssDF1, encmod="P") #"global" model
m30a 

#31 session specific density, sex specific encounter rate and sex and session dependent sigma
m31a<-oSCR.fit(list(D~session, p0~sex, sig~session + sex),data1$sf,data1$ssDF1, encmod="P") #"global" model
m31a 

### Model selection #### 

f2a = fitList.oSCR(list(mod.oSCR1, m2a, m3a, m4a, m5a, m6a, m7a, m8a, m9a, m10a, m11a, m12a, m13a, m14a, m15a, m16a, m17a, m18a, m19a, m20a,  
                        m21a, m22a, m23a, m24a, m25a, m26a, m27a, m28a, m29a, m30a, m31a), 
                   drop.asu = T, rename = TRUE) # rename = T adds sensible model names 
f2a
# ms = modSel.oSCR(f1a)
ms1 <- modSel.oSCR(f2a)
ms1

#print.oSCR.fit(m16a)

########### get estimates ########################

#betas
print.oSCR.fit(m16a)  

top.model <-m16a
pred.df <-data.frame(session =factor(c(1,2,3))) #for density and sigma
pred.df.det<-data.frame(session =factor(c(1,2,3,1,2,3)),sex=factor(c(0,0,0,1,1,1)))

#overall density
pred.dens <- get.real(model =top.model,type ="dens",newdata =pred.df, N_sex=F)
pred.dens

#density by sex
pred.dens_sex <- get.real(model =top.model,type ="dens",newdata =pred.df, N_sex=T)
pred.dens_sex

#difference in space use by session
pred.sig <-get.real(model =top.model,type ="sig",newdata =pred.df)
pred.sig

#difference in encounter rate by sex and session
pred.det<-get.real(model = top.model, type = "det", newdata = pred.df.det, N_sex=T)
pred.det


