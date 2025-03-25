
var_dist <- function(continuous_array) { 
     continuous_array <- continuous_array[!is.na(continuous_array)]
     return( c(summary(continuous_array),length(continuous_array),sd(continuous_array)) ) }

baseline_lab_distribution <- function(mtab){
     bl <- data.frame()
     for(pheno in c("UPCR","logUPCR","eGFR","age")){
          bl <- rbind(bl,c(pheno,var_dist(mtab[,pheno])) ) }
     bl <- rbind(bl,c("UACR",var_dist(mtab[!mtab$ID %in% pn1$ID,"UPCR"])) ) 
     colnames(bl) <- c("pheno","Min","1stQuartile","Median","Mean","3rdQuartile","Max","N","SD")
     return(bl) } 

output_baseline_table <- function(mtab){
     bl <- data.frame()
     bl <- rbind(bl,c("N",nrow(mtab)))
     c1 <- paste0(nrow(mtab[mtab$Class=="AFR",])," (",round(nrow(mtab[mtab$Class=="AFR",])/nrow(mtab)*100,digits=1),"%)")
     bl <- rbind(bl,c("Ancestry",c1))
     c1 <- paste0(nrow(mtab[mtab$GENDER=="Female",])," (",round(nrow(mtab[mtab$GENDER=="Female",])/nrow(mtab)*100,digits=1),"%)")
     bl <- rbind(bl,c("GENDER",c1))

     for(pheno in c("age","UPCR","logUPCR","eGFR")){
          c1 <- paste0(round(mean(mtab[,pheno]),digits=2)," (",round(sd(mtab[,pheno]),digits=2),")" )
          bl <- rbind(bl,c(pheno,c1)) }
     for(pheno in c("CVD","HF","HTN","ACEi","ARB","CCB","TD","t2d","Insulin","SGLT2","ACEiandARB","ACEiorARB","ACEionly","ARBonly","SGLT2only","ACEi_ARB_SGLT2") ){
          c1 <- paste0(sum(mtab[,pheno])," (",round(sum(mtab[,pheno])/nrow(mtab)*100,digits=1),"%)")
          bl <- rbind(bl,c(pheno,c1)) }
     for(ckds in c(2:4) ){
          c1 <- paste0(nrow(mtab[mtab$CKD==ckds,])," (",round(nrow(mtab[mtab$CKD==ckds,])/nrow(mtab)*100,digits=1),"%)")
          bl <- rbind(bl,c(paste0("CKD stage",ckds),c1)) }
     colnames(bl) <- c("pheno","ss")
     return(bl) }


PCR_change_bymonths <- function(PCRdata,baseline,measument){
     PCRdata$date <- as.Date(PCRdata$date,format = "%m/%d/%Y")
     baseline$upcr_date <- as.Date(baseline$upcr_date)
     PCRdata <- inner_join(PCRdata,baseline,by="ID",relationship = "many-to-many")
     PCRdata$gap <- PCRdata$date-PCRdata$upcr_date

     pcr0 <- PCRdata[PCRdata$gap<30 & PCRdata$gap > -30,] %>% group_by(ID,upcr_date) %>% summarise(PCR0 = mean(UPCR))
     mpcr <- pcr0
     for(nmonth in c(6,9,12,24)){
          pcrt <- PCRdata[PCRdata$gap>(365*nmonth/12-61) & PCRdata$gap<(365*nmonth/12+61),] %>% group_by(ID,upcr_date) %>% summarise( !! paste0('PCR',nmonth) := mean(UPCR) )
          mpcr <- left_join(mpcr,pcrt,by=c("ID","upcr_date"))
          mpcr[,paste0('actualPCR',nmonth,'C')] <- mpcr[,paste0('PCR',nmonth)]-mpcr[,'PCR0']
          mpcr[,paste0('logPCR',nmonth,'C')] <- log(mpcr[,paste0('PCR',nmonth)])-log(mpcr[,'PCR0']) }
     mpcr$measument <- measument
     return(mpcr) }

GMC <- function(var_array) {
     var_array <- var_array[!is.na(var_array)]
     meanoc <- mean(var_array)
     nsam = length(var_array)
     sdoc <- sd(var_array)/sqrt(nsam)
     return(c(nsam,meanoc,sdoc,(exp(meanoc) -1)*100,(exp(meanoc-qt(0.975,nsam-1)*sdoc)-1)*100,(exp(meanoc+qt(0.975,nsam-1)*sdoc)-1)*100) )}


individual_egslope <- function(egtab){
     dat <- idd[!duplicated(idd$ID),c("ID","epi_date")]
     yn <- c("x1"="1year","x2"="2years", "x3"="3years", "x100"="All")
     for(GFRyear in c(1,2,3,100)){
          sub_eGFR <- egtab[egtab$timefromidd <= GFRyear,]
          sub_eGFR <- sub_eGFR[sub_eGFR$ID %in% sub_eGFR[duplicated(sub_eGFR$ID),]$ID,]
          for(eGFRcol in c('eGFR','eGFR2')){
               cformula <- as.formula(paste(eGFRcol,'~ timefromidd + (timefromidd|ID)'))
               lmdata <- data.frame(coef(lmer(cformula,data=sub_eGFR))$ID)
               lmdata$ID<- row.names(lmdata)
               colnames(lmdata)[2] <- paste0(eGFRcol,'_slope_',yn[paste0("x",GFRyear)])
               dat <- left_join(dat,lmdata[,2:3],by="ID")
          } }
     return(dat) }

egslope <- function(sampleIDs){
     egtab <- eg[eg$ID %in% sampleIDs,]
     res <- data.frame()
     covariates <- c("age", "eGFR0", "age + eGFR0")

     for(eGFRcol in c('eGFR','eGFR2')){
          for(GFRyear in c(1,2,3,4,100)){
               mtab <- egtab[egtab$timefromidd <= GFRyear,]
               mtab2 <- mtab[order(-mtab$ed30),]
               mtab2 <- mtab2[!duplicated(mtab2$ID),]
               mtab3 <- mtab[order(-mtab$ed40),]
               mtab3 <- mtab3[!duplicated(mtab3$ID),]

               mtab <- mtab[mtab$ID %in% mtab[duplicated(mtab$ID),]$ID,]
               cformula <- as.formula(paste(eGFRcol, '~ timefromidd + (timefromidd|ID)'))
               fit <- lmer(cformula, data = mtab)
               cis = c(NA,NA)
               tryCatch({cis = confint(fit,"timefromidd", level=0.95)},error=function(e) {
                              cat("An error occurred:", conditionMessage(e), "\n") })
               res <- rbind(res,c(GFRyear,eGFRcol,coef(summary(fit))[2,][c(1,2,5)],cis,length(unique(mtab$ID)),sum(mtab2$ed30),sum(mtab2$ed40),"None") )

               for(cov in covariates){
                    cformula <- as.formula(paste(eGFRcol, '~ timefromidd + (timefromidd|ID)','+',cov))
                    fit <- lmer(cformula, data = mtab)
                    cis = c(NA,NA)
                    tryCatch({ cis = confint(fit,"timefromidd", level=0.95)},error=function(e) {
                                   cat("An error occurred:", conditionMessage(e), "\n") })
                    res <- rbind(res,c(GFRyear,eGFRcol,coef(summary(fit))[2,][c(1,2,5)],cis,length(unique(mtab$ID)),sum(mtab2$ed30),sum(mtab2$ed40),cov) ) }
          } } 
     colnames(res) <- c("eGFR_duration","eGFRapproach","beta","se","P","CI2.5","CI97.5","Nsam","ED30","ED40","covariate")
     return(res) }

data_for_association <- function(sub_eGFR){
     sub_eGFR <- sub_eGFR[sub_eGFR$ID %in% sub_eGFR[duplicated(sub_eGFR$ID),]$ID,]
     eGFR_slope <- lmer(eGFR ~ timefromidd + (timefromidd|ID), data = sub_eGFR)
     lmdata <- data.frame(coef(eGFR_slope)$ID)
     lmdata$ID<- row.names(lmdata)
     lmdata <- inner_join(lmdata, pcrtab,by="ID")
     return(lmdata) }

UPCR_months <-  c(6,9,12,24)
eGFR_months <-  c(24,36,360)

association_test <- function(sampleIDs){
     covariates <- c("age", "logUPCR0", "age + logUPCR0")
     res <- data.frame()
     for(em in eGFR_months){
          lmdata <- data_for_association(eg[eg$ID %in% sampleIDs & eg$timefromidd <= em/12,])
          for(um in UPCR_months){
               for(xx in c('actual','log')){
                    expo = paste0(xx,"PCR",um,"C")
                    cformula <- as.formula(paste('timefromidd ~',expo))
                    fit <- lm(cformula, data = lmdata)
                    res <- rbind(res, c(lmfit_coef(fit),sum(!is.na(lmdata[,expo])),em,um,xx,"None") )

                    for(cov in covariates){
                         cformula <- as.formula(paste('timefromidd ~',expo,'+',cov))
                         fit <- lm(cformula, data = lmdata)
                         res <- rbind(res, c(lmfit_coef(fit),sum(!is.na(lmdata[,expo])),em,um,xx,cov) )
                    }
               }
          }
     }
     colnames(res) <- c("beta_intercept","CI_intercept","se_intercept","P_intercept","beta_expo","CI_expo","SE_expo","P_expo","R","Npatients","Duration_of_eGFR_slope","Duration_of_UPCR_change","UPCR used","covariate")
     return(res) 
}

    
lmfit_coef <- function(lmfit){
     x <- summary(lmfit)
     cis = confint(lmfit, level=0.95)

     beta_intercept <- signif(x$coef[1,][1], digits=2)
     se_intercept <- signif(x$coef[1,][2], digits=2)
     CI_intercept <- paste0("(",signif(cis[1,][1],digits=2),"~",signif(cis[1,][2],digits=2),")")

     beta1 <- signif(x$coef[2,][1], digits=2)
     se1 <- signif(x$coef[2,][2], digits=2)
     CI1 <- paste0("(",signif(cis[2,][1],digits=2),"~",signif(cis[2,][2],digits=2),")")
     return(c(beta_intercept,CI_intercept,se_intercept,x$coef[1,][4],beta1,CI1,se1,x$coef[2,][4],x$r.squared))  
}

association_log <- function(lmdata,UPCRmonth){
     cformula <- as.formula(paste0('timefromidd ~ PCR',UPCRmonth,"C"))
     fit <- lm(cformula, data = lmdata)
     return(summary(fit))  }

scatter_plot <- function(lmdata,eGFRmonth,UPCRmonth,filename){
     splot <- ggplot(lmdata,aes(x= lmdata[,paste0('PCR',UPCRmonth,"C")],y=timefromidd)) + geom_point() +
               labs(x = "Change of log(UPCR)",y = "eGFR slope") + 
               #scale_x_continuous(breaks = scales::pretty_breaks(n = 20),limits=c(0,15)) +
               #scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,1)) +
               ggtitle(paste("change of log(UPCR) in ",UPCRmonth,"months\n","eGFR slope in ",eGFRmonth,"months")) + 
               theme(plot.title = element_text(hjust = 0.5)) + 
               geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkgrey")
     ggsave(file = filename, print(splot),width=4,height=4)  }

#output UPCR change
GMC_output <- function(sub_pcrtab){
     PCR_GMC <- data.frame()
     #PCR_dist <- data.frame()
     for(mopcrc in c(6,9,12,24)){
          PCR_GMC <- rbind(PCR_GMC,c(mopcrc,GMC(sub_pcrtab[,paste0("PCR",mopcrc,"C")])) )
          #PCR_dist <- rbind(PCR_dist,c(mopcrc,var_dist(sub_pcrtab[,paste0("PCR",mopcrc)])) ) 
     }
     colnames(PCR_GMC) <- c("Umonth","N","mean","sd","Gmean","CI1","CI2")
     #colnames(PCR_dist) <- c("Umonth","Min","1stQuartile","Median","Mean","3rdQuartile","Max","N","SD")
     return(PCR_GMC)
}



CoxPH_coef <- function(Coxmodel,dependent){ 
     sumx <- summary(Coxmodel)
     #p.value <- signif(x$wald["pvalue"], digits=2)
     #beta <- signif(x$coef[dependent,][1], digits=2)
     #SE <- signif(x$coef[dependent,][3], digits=2)
     HR <-signif(sumx$coef[dependent,][2], digits=2)
     HR.confint.lower <- signif(sumx$conf.int[dependent,"lower .95"], 2)
     HR.confint.upper <- signif(sumx$conf.int[dependent,"upper .95"],2)
     HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
     #res<-c(beta, SE, HR, p.value, x$n, x$nevent)
     res<-c(sapply(sumx$coef[dependent,][c(1,3,5)], function (x) signif(x,digits=2)),HR,sumx$n,sumx$nevent)
     return(res) }


CoxPH <- function(dat){
     exposures1 <- c("actualPCR6C","actualPCR9C","actualPCR12C","actualPCR24C","logPCR6C","logPCR9C","logPCR12C","logPCR24C","UI30m6","UI30m9","UI30m12","UI30m24","UI30")
     exposures2 <- c("eGFR_slope_1year","eGFR_slope_2years","eGFR_slope_3years","eGFR_slope_All","ed30","ed30m12","ed30m24","ed30m36","ed40","ed40m12","ed40m24","ed40m36")

     res <- data.frame()
     for(expo in exposures1){
          for(status in c("ESKD","Clinical Endpoint")){

               cformula <- as.formula(paste('Surv(time,event) ~',expo))
               fit <- coxph(cformula,data=dat[dat$type==status,])
               tryCatch({res <- rbind(res,c(CoxPH_coef(fit,expo),status,expo,"None",rep(NA,8)))},error=function(e) {
                              cat("An error occurred:", conditionMessage(e), "\n") })

               cformula <- as.formula(paste('Surv(time,event) ~',expo,'+ age'))
               fit <- coxph(cformula,data=dat[dat$type==status,])
               tryCatch({res <- rbind(res,c(CoxPH_coef(fit,expo),status,expo,'age',CoxPH_coef(fit,'age')[1:4],rep(NA,4)))},error=function(e) {
                         cat("An error occurred:", conditionMessage(e), "\n") })
               
               cformula <- as.formula(paste('Surv(time,event) ~',expo,'+ logUPCR0'))
               fit <- coxph(cformula,data=dat[dat$type==status,])
               tryCatch({res <- rbind(res,c(CoxPH_coef(fit,expo),status,expo,'logUPCR0',rep(NA,4),CoxPH_coef(fit,'logUPCR0')[1:4]))},error=function(e) {
                         cat("An error occurred:", conditionMessage(e), "\n") })

               cformula <- as.formula(paste('Surv(time,event) ~',expo,'+ age + logUPCR0'))
               fit <- coxph(cformula,data=dat[dat$type==status,])
               tryCatch({res <- rbind(res,c(CoxPH_coef(fit,expo),status,expo,'age + logUPCR0',CoxPH_coef(fit,'age')[1:4],CoxPH_coef(fit,'logUPCR0')[1:4]))},error=function(e) {
                         cat("An error occurred:", conditionMessage(e), "\n") })
          }
     }

     for(expo in exposures2){

          cformula <- as.formula(paste('Surv(time, event) ~',expo))
          fit <- coxph(cformula,data=dat[dat$type=="ESKD",])
          tryCatch({res <- rbind(res,c(CoxPH_coef(fit,expo),"ESKD",expo,"None",rep(NA,8)))},error=function(e) {
                              cat("An error occurred:", conditionMessage(e), "\n") })

          cformula <- as.formula(paste('Surv(time,event) ~',expo,'+ age'))
          fit <- coxph(cformula,data=dat[dat$type=="ESKD",])
          tryCatch({res <- rbind(res,c(CoxPH_coef(fit,expo),"ESKD",expo,'age',CoxPH_coef(fit,'age')[1:4],rep(NA,4)))},error=function(e) {
                    cat("An error occurred:", conditionMessage(e), "\n") })
          
          cformula <- as.formula(paste('Surv(time,event) ~',expo,'+ eGFR0'))
          fit <- coxph(cformula,data=dat[dat$type=="ESKD",])
          tryCatch({res <- rbind(res,c(CoxPH_coef(fit,expo),"ESKD",expo,'eGFR0',rep(NA,4),CoxPH_coef(fit,'eGFR0')[1:4]))},error=function(e) {
                    cat("An error occurred:", conditionMessage(e), "\n") })

          cformula <- as.formula(paste('Surv(time,event) ~',expo,'+ age + eGFR0'))
          fit <- coxph(cformula,data=dat[dat$type=="ESKD",])
          tryCatch({res <- rbind(res,c(CoxPH_coef(fit,expo),"ESKD",expo,'age + eGFR0',CoxPH_coef(fit,'age')[1:4],CoxPH_coef(fit,'eGFR0')[1:4]))},error=function(e) {
                    cat("An error occurred:", conditionMessage(e), "\n") })
     }
     
     colnames(res)<-c("beta","SE","p_value","HR (95% CI)","N","Nevent","status","exposure","covariate","beta_age","SE_age","p_value_age","HR_age (95% CI)", "beta_cov2","SE_cov2", "p_value_cov2", "HR_cov2 (95% CI)")
     return(res) 
}



#----------------------------------------------------------------------------------------------------------------------------------------

Predict_CoxPH <- function(dat){
     exposures1 <- c("logPCR6C","logPCR9C","logPCR12C","logPCR24C")
     chg_lg <- c(0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0)
     #chg_lg0 <- c(0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0)
     res <- data.frame()
     
     for(status in c("ESKD","Clinical Endpoint")){
          sub_data = dat[dat$type==status,]

          for(expo in exposures1){
               #chg_lg <- chg_lg0+mean(sub_data[,expo],na.rm=T)

               cformula <- as.formula(paste('Surv(time,event) ~',expo))
               fit <- coxph(cformula,data=sub_data)
               newdata <- data.frame(expo=chg_lg)
               colnames(newdata)[1] <- expo
               tryCatch({res <- rbind(res,c(predict(fit,newdata, type="risk"),status,expo,"None",CoxPH_coef(fit,expo)))},error=function(e) {
                              cat("An error occurred:", conditionMessage(e), "\n") })

               cformula <- as.formula(paste('Surv(time,event) ~',expo,'+ age'))
               fit <- coxph(cformula,data=sub_data)
               newdata <- data.frame(expo=chg_lg,age=mean(sub_data[!is.na(sub_data[,expo]),]$age))
               colnames(newdata)[1] <- expo
               tryCatch({res <- rbind(res,c(predict(fit,newdata, type="risk"),status,expo,'age',CoxPH_coef(fit,expo)))},error=function(e) {
                         cat("An error occurred:", conditionMessage(e), "\n") })
               
               cformula <- as.formula(paste('Surv(time,event) ~',expo,'+ logUPCR0'))
               fit <- coxph(cformula,data=sub_data)
               newdata <- data.frame(expo=chg_lg,logUPCR0=mean(sub_data[!is.na(sub_data[,expo]),]$logUPCR0))
               colnames(newdata)[1] <- expo
               tryCatch({res <- rbind(res,c(predict(fit,newdata, type="risk"),status,expo,'logUPCR0',CoxPH_coef(fit,expo)))},error=function(e) {
                         cat("An error occurred:", conditionMessage(e), "\n") })

               cformula <- as.formula(paste('Surv(time,event) ~',expo,'+ age + logUPCR0'))
               fit <- coxph(cformula,data=sub_data)
               newdata <- data.frame(expo=chg_lg,age=mean(sub_data[!is.na(sub_data[,expo]),]$age),logUPCR0=mean(sub_data[!is.na(sub_data[,expo]),]$logUPCR0))
               colnames(newdata)[1] <- expo
               tryCatch({res <- rbind(res,c(predict(fit,newdata, type="risk"),status,expo,'age + logUPCR0',CoxPH_coef(fit,expo)))},error=function(e) {
                         cat("An error occurred:", conditionMessage(e), "\n") })
          }
     }
     
     colnames(res)<-c(paste0('HR',c(1:length(chg_lg))),"status","exposure","covariate","beta","SE","p_value","HR","N","Nevent")
     return(res) 
}
