####################
# Generate Figures #
####################

out<-"Store MCMC Samples Directory"  # Store Samples

# Results
Alpha1<-matrix(scan(paste(out,"Alpha1.txt",sep="")),ncol=p,byrow=T)
Alpha2<-matrix(scan(paste(out,"Alpha2.txt",sep="")),ncol=p,byrow=T)
Beta1<-matrix(scan(paste(out,"Beta1.txt",sep="")),ncol=k,byrow=T)
Beta2<-matrix(scan(paste(out,"Beta2.txt",sep="")),ncol=k,byrow=T)
R1<-matrix(scan(paste(out,"R1.txt",sep="")),ncol=1,byrow=T)
R2<-matrix(scan(paste(out,"R2.txt",sep="")),ncol=1,byrow=T)
Sigmab<-matrix(scan(paste(out,"Sigmab.txt",sep="")),nrow=lastit,byrow=T)
Sigma<-matrix(scan(paste(out,"Sigma.txt",sep="")),nrow=lastit,byrow=T)
Sigmab<-matrix(scan(paste(out,"Sigmab.txt",sep="")),nrow=lastit,byrow=T)
Bs1<-matrix(scan(paste(out,"B1.txt",sep="")),nrow=lastit,byrow=T)
Bs2<-matrix(scan(paste(out,"B2.txt",sep="")),nrow=lastit,byrow=T)
Ls<-matrix(scan(paste(out,"L.txt",sep="")),nrow=lastit,byrow=T)

#######################
# Diagnostics: Fig S1 #
#######################

#pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\FigS1.pdf")
par(mfrow=c(3,2))
plot(1:lastit,Beta1[,1],type="l",col="lightgreen", main=expression(beta[11]),
     xlab="Iteration",ylab=expression(beta[11]))
legend(cex=.7,"topleft",legend=c("ESS = 2000","p = 0.78"))
abline(h=mbeta1[1])
effectiveSize(Beta1[,1])
geweke.diag(Beta1[,1])
2*(1-pnorm(geweke.diag(Beta1[,1])$z))

plot(1:lastit,Beta2[,1],type="l",col="lightgreen",main=expression(beta[21]),
     xlab="Iteration",ylab=expression(beta[21]))
legend(cex=.7,"topleft",legend=c("ESS = 2000","p = 0.88"))
abline(h=mbeta2[1])
effectiveSize(Beta2[,1])
geweke.diag(Beta2[,1])
2*(pnorm(geweke.diag(Beta2[,1])$z))

plot(1:lastit,R1,type="l",col="lightgreen",main=expression(r[1]),
     xlab="Iteration",ylab=expression(r[1]))
legend(cex=.7,"topleft",legend=c("ESS = 1131","p = 0.68"))
abline(h=mr1)
effectiveSize(R1)
geweke.diag(R1)
2*pnorm(geweke.diag(R1)$z)

plot(1:lastit,R2,type="l",col="lightgreen",main=expression(r[2]),
     xlab="Iteration",ylab=expression(r[2]))
legend(cex=.7,"topleft",legend=c("ESS = 917","p = 0.96"))
abline(h=mr2)
effectiveSize(R2)
geweke.diag(R2)
2*(1-pnorm(geweke.diag(R2[,1])$z))

plot(1:lastit,Sigma[,4],type="l",col="lightgreen",main=expression(Sigma[beta[22]]),
     xlab="Iteration",ylab=expression(Sigma[beta[22]]))
legend(cex=.7,"topleft",legend=c("ESS = 2000","p = 0.10"))
abline(h=msigma[4])
effectiveSize(Sigma[,4])
2*(1-pnorm(geweke.diag(Sigma[,4])$z))

plot(1:lastit,Sigmab[,1],type="l",col="lightgreen",main=expression(Sigma[b[11]]),
     xlab="Iteration",ylab=expression(Sigma[b[11]]))
legend(cex=.7,"topleft",legend=c("ESS = 462","p = 0.86"))
abline(h=msigmab[1])
effectiveSize(Sigmab[,1])
2*(pnorm(geweke.diag(Sigmab[,1])$z))
dev.off()

###############################################
# Pop Average and County-Specific Predictions #
###############################################
Pred1<-Pred2<-array(0,dim=c(ncounty,lday,lastit))
for (j in 1:lastit){
  for (h in 1:ncounty){
    Bs21<-matrix(Bs1[j,],ncol=k,byrow=T)
    Bs22<-matrix(Bs2[j,],ncol=k,byrow=T)
    beta1<-Beta1[j,]
    beta2<-Beta2[j,]
    r1<-R1[j]
    r2<-R2[j]
    alpha1<-Alpha1[j,]
    alpha2<-Alpha2[j,]
    Pred1[h,,j]<-r1*c(exp(alpha1[1]*mean(W[id==h,1])+alpha1[2]*mean(W[id==h,2])+Xtmp%*%beta1+Xtmp%*%Bs21[h,]))     # County x Day Predictions (per capita mean per 100k); don't include alpha and W for covariate-adjusted predictions
    Pred2[h,,j]<-r2*c(exp(alpha2[1]*mean(W[id==h,1])+alpha2[2]*mean(W[id==h,2])+Xtmp%*%beta2+Xtmp%*%Bs22[h,]))     # County x Day Predictions (per capita mean per million)
  } 
  print(j)
}

pred1<-colMeans(t(colMeans(Pred1[,,1:lastit])))     # PA prediction (average across counties and iterations)
pred2<-colMeans(t(colMeans(Pred2[,,1:lastit])))     # PA prediction (average across counties and iterations)

countypred1<-apply(Pred1[,,1:lastit],c(1,2),mean)   # County-specific predictions, i.e., post mean for each county/day
countypred2<-apply(Pred2[,,1:lastit],c(1,2),mean)   # County-specific predictions, i.e., post mean for each county/day

#####################################
# Credible Intervals for Mean Trend #
#####################################
# County specific
countylower1<-apply(Pred1[,,1:lastit],c(1,2),quantile,prob=0.025)  # Lower CI for each county 
countyupper1<-apply(Pred1[,,1:lastit],c(1,2),quantile,prob=0.975)  # Average lower bounds across county 

countylower2<-apply(Pred2[,,1:lastit],c(1,2),quantile,prob=0.025)  # Lower CI for each county 
countyupper2<-apply(Pred2[,,1:lastit],c(1,2),quantile,prob=0.975)  # Average lower bounds across county

# Population average
palower1<-apply(t(colMeans(Pred1[,,1:lastit])),2,quantile,prob=0.025)
paupper1<-apply(t(colMeans(Pred1[,,1:lastit])),2,quantile,prob=0.975)

palower2<-apply(t(colMeans(Pred2[,,1:lastit])),2,quantile,prob=0.025)
paupper2<-apply(t(colMeans(Pred2[,,1:lastit])),2,quantile,prob=0.975)

############################
# Population Average Plots #
############################

plot1a_data <- tibble(
  day = 1:lday,
  orr = orr1,
  postmean = pred1,
  palower = palower1,
  paupper = paupper1
)

plot1b_data <- tibble(
  day = 1:lday,
  orr = orr2,
  postmean = pred2,
  palower = palower2,
  paupper = paupper2
)

###############################
# Fig 1a: Pop Avg Plot for y1 #
###############################
#pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\Fig1a.pdf")
ggplot(plot1a_data, aes(x = day, y = orr ))+
  geom_point(aes(color = "Daily Average"))+
  geom_line(aes(day, postmean, color="Posterior Mean"),size=1,linetype=1)+
  geom_ribbon(aes(ymin = palower, ymax = paupper, color="95% Credible Interval"),alpha=0.3, show.legend = F)+
  scale_color_manual(breaks = c("Daily Average", 
                                "Posterior Mean", 
                                "95% Credible Interval"
  ), 
  values = c("Daily Average"="darkorange",
             "Posterior Mean"="blue4",
             "95% Credible Interval" = "gray40"
  ))+
  ylab("Cases per 100K")+
  #ylim(0,6)+
  guides( color = guide_legend(title="Incidence Rate Trends",
                               override.aes = list(
                                 linetype = c(0, 1, 1),
                                 shape = c(19,NA, NA )),
                               reverse = F))+
  
  theme(legend.key = element_rect(fill = "white"), legend.position = c(0.85, 0.85),legend.text=element_text(size=11))

#dev.off()

###############################
# Fig 1b: Pop Avg Plot for y2 #
###############################
#pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\Fig1b.pdf")
ggplot(plot1b_data, aes(x = day, y = orr ))+
  geom_point(aes(color = "Daily Average"))+
  geom_line(aes(day, postmean, color="Posterior Mean"),size=1,linetype=1)+
  geom_ribbon(aes(ymin = palower, ymax = paupper, color="95% Credible Interval"),alpha=0.3, show.legend = F)+
  scale_color_manual(breaks = c("Daily Average", 
                                "Posterior Mean", 
                                "95% Credible Interval"
  ), 
  values = c("Daily Average"="darkorange",
             "Posterior Mean"="blue4",
             "95% Credible Interval" = "gray40"
  ))+
  ylab("Deaths per Million")+
  scale_y_continuous(name="Date",breaks=c(0,5,10),limits=c(0,10))+
  # ylim(0,10)+
  guides( color = guide_legend(title="Death Rate Trends",
                               override.aes = list(
                                 linetype = c(0, 1, 1),
                                 shape = c(19,NA, NA )),
                               reverse = F))+
  #ggtitle("Population Average Death Trend Across All Counties")+
  theme(legend.key = element_rect(fill = "white"), legend.position = c(0.85, 0.85),legend.text=element_text(size=11))

#dev.off()

######################################
# Figs 2a-2d: County Specific Plots  #
######################################
#############
# Y1: Cases #
#############
countyid<-sample(1:ncounty,1)  # Random select 
print(countyid) 

# County 70, 151
countyid<-151 # 70 

simcurve1<-exp(W[id==countyid,]%*%truealpha1+Xtmp%*%truebeta1+Xtmp%*%trueB1[countyid,]) # True, Simulated Curve for Y1

# County inc data
county_inc <- tibble(
  day = 1:lday,
  orr1 = incidence1[id==countyid],
  postmean = countypred1[countyid, ],
  countyl = countylower1[countyid, ],
  countyu = countyupper1[countyid, ],
  simcurve = simcurve1
)


# Incidence Plot

# pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\Fig2c.pdf")
ggplot(county_inc, aes(x = day, y = orr1 ))+
  geom_point(aes(color = "Daily Incidence Rate"), shape = 16)+
  geom_line(aes(day, simcurve, color="Simulated Curve"),size=1,linetype=2)+
  geom_line(aes(day, postmean, color="Posterior Mean"),size=1,linetype=1)+
  geom_ribbon(aes(ymin = countyl, ymax = countyu, color="95% Credible Interval"),alpha=0.3, show.legend = F,
              linetype=1)+
  scale_color_manual(breaks = c("Daily Incidence Rate",
                                "Simulated Curve",
                                "Posterior Mean", 
                                "95% Credible Interval"), 
                     values = c("Daily Incidence Rate"="darkorange",
                                "Simulated Curve"="darkgreen",
                                "Posterior Mean"="blue4",
                                "95% Credible Interval" = "gray40"
                     ))+
  ylab("Cases per 100K")+
  xlab("Day")+
  guides( color = guide_legend(title="Incidence Rate Trends",
                               override.aes = list(
                                 linetype = c(0, 2,1, 1),
                                 shape = c(16,NA,NA, NA )),
                               reverse = F))+
  # ylim(0, 15)+
  # scale_y_continuous(name="Date",breaks=c(0,10,20))+
  # ggtitle("Cases")+
  theme(legend.key = element_rect(fill = "white"),
        legend.title = element_text(size = 11),
        plot.title=element_text(hjust=.5),
        legend.position = c(0.85, 0.85))

# dev.off()

###############
# Y2: Deaths  #
###############
simcurve2<-exp(W[id==countyid,]%*%truealpha2+Xtmp%*%truebeta2+Xtmp%*%trueB2[countyid,]) # True, simulated Curve for y2

county_death <- tibble(
  day = 1:lday,
  orr1 = incidence2[id==countyid],
  postmean = countypred2[countyid, ],
  countyl = countylower2[countyid, ],
  countyu = countyupper2[countyid, ],
  simcurve = simcurve2
)

# Deaths 

# pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\Fig2d.pdf")
ggplot(county_death, aes(x = day, y = orr1 ))+
  geom_point(aes(color = "Daily Death Rate"), shape = 16)+
  geom_line(aes(day, simcurve, color="Simulated Curve"),size=1,linetype=2)+
  geom_line(aes(day, postmean, color="Posterior Mean"),size=1,linetype=1)+
  geom_ribbon(aes(ymin = countyl, ymax = countyu, color="95% Credible Interval"),alpha=0.3, show.legend = F,linetype=1)+
  scale_color_manual(breaks = c("Daily Death Rate", 
                                "Simulated Curve",
                                "Posterior Mean", 
                                "95% Credible Interval"), 
                     values = c("Daily Death Rate"="darkorange",
                                "Simulated Curve"="darkgreen",
                                "Posterior Mean"="blue4",
                                "95% Credible Interval" = "gray40"
                     ))+
  ylab("Deaths per Million")+
  xlab("Day")+
  guides( color = guide_legend(title="Death Rate Trends",
                               override.aes = list(
                                 linetype = c(0, 2,1, 1),
                                 shape = c(16,NA,NA, NA )),
                               reverse = F))+
  
  #scale_y_continuous(breaks=c(0,10,20,30))+
  #ylim(0, 20)+
  
  #ggtitle("Deaths")+
  theme(legend.key = element_rect(fill = "white"),
        legend.title = element_text(size = 11),
        plot.title=element_text(hjust=.5),
        legend.position = c(0.85, 0.85))

# dev.off()

# Save files and workspaces
#save.image("C:\\Brian\\Covid\\R\\Workspaces\\Simulation\\Simulation_07-30-21.RData")
#load("C:\\Brian\\Covi\\R\\Workspaces\\Simulation\\Simulation_07-30-21.RData")
