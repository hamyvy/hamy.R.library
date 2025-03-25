library(ggplot2)

##################### Forest plot
# method #1
library(forestplot)

tabletext <- cbind(c( "Low", "Intermediate", "High"))
clrs <- fpColors(box = "darkred",line = "darkred", summary = "red")

mydat %>%
  forestplot(labeltext = tabletext,fn.ci_norm = fpDrawCircleCI,
             boxsize = 0.2, zero = 1,
             col = clrs,grid = TRUE,
             txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex = 1),
                              xlab  = gpar(fontfamily = "", cex = 1.3)),
             xlab = "Adjusted Odds ratio (95% Confidence Interval)",vertices = TRUE)

# method #2
p = ggplot(data=cph,aes( x = strata,y = HR, ymin = CI1, ymax = CI2 ))+
    geom_pointrange()+ geom_hline(yintercept =1, linetype=2)+ 
    xlab('s_het interval')+ ylab("Hazard Ratio (95% Confidence Interval)")+ 
    geom_text(aes(x=tt+0.2, y=HR, label=paste0("p=",pval)))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
      axis.line = element_line(colour = "black"), plot.title=element_text(size=16,face="bold"),
      axis.text=element_text(face="bold"),axis.title=element_text(size=14,face="bold"))+
    coord_flip()

    
##################### Kaplan Meiyer plot
survp <- survfit2(Surv(time, event) ~ T2Dn, data = df[df$type==pheno & df$severity==sv,]) %>% ggsurvfit(linewidth = 0.8) +
                    labs(x = "Years",y = "Overall probability" ) + add_confidence_interval() + 
                    add_risktable() + scale_x_continuous(breaks = scales::pretty_breaks(n = 20),limits=c(0,15)) +
                    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,1)) +
                    ggtitle(paste0(pheno,'-',ptitle2[[sv]])) + theme(plot.title = element_text(hjust = 0.5))+
                    scale_color_manual(values = c('purple', 'cyan3')) + scale_fill_manual(values = c('purple', 'cyan3'))
ggsave(file = paste0('PCRgap',nmonth,'/Surv2_',pheno,'_',sv,'_non-diabetics_VS_diabetics_without_DNdef1.pdf'), print(survp)) 


##################### confidence interval plot as ribbon
p <- ggplot(mydat, aes(x =bin , y = OR)) + geom_ribbon(aes(ymin = OR1,
     ymax = OR2, fill = cohort), alpha = 0.2) + geom_point(aes(colour = cohort),
     size = 1,show.legend = FALSE) + labs(title="")+
     xlab('Percentile of PRS')+ ylab("Predicted odds ratio")+
     scale_fill_discrete(name="",labels=c((expression("Bio"*italic("Me")*"-AA")),(expression("Bio"*italic("Me")*"-HA")),"UKBB-AA"))+
     theme_bw()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.15, 0.85))+ theme(legend.title = element_blank(),legend.text = element_text(size=12))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
       axis.line = element_line(colour = "black"), plot.title=element_text(size=16,face="bold"),
       axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold"))


##################### add error bar
p1 <- ggplot(d1,aes(x=centre,y=all))+geom_errorbar(aes(ymin=all-se,ymax=all+se),width=.1,colour="grey")+
      geom_point(size=2)+ labs(title="Genome-wide load score",x="",y="Load score (normalized)")+ 
      scale_y_continuous(expand = c(0,0),limits=c(-0.1,0.2),breaks = c(-0.1,0,0.1,0.2))+
      theme(plot.title = element_text(size=16,face="bold",hjust=0.5), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),axis.text = element_text(face="bold",size=11,colour="black"),axis.text.x=element_text(angle=60, hjust = 1),
            axis.title.x = element_text(size=14, face="bold",vjust = -1.5),axis.title.y = element_text(size=14, face="bold"))+
      geom_text(x=11, y=0.19, label=expression(paste("p-ANOVA = 4.7x10"^{-25})),size=4.5,family="Times")


p2 <- ggplot(xx[xx$Disease=="daly",],aes(ylab,beta))+geom_pointrange(aes(ymin =lower, ymax =upper))+ 
    geom_hline(yintercept =0, linetype=2,col="grey")+ xlab('')+ ylab("Mean year difference") +
    geom_text(aes(x=ylab, y=beta, label=paste0("p=",P)),position = position_dodge(0.5),vjust = -0.6,size=4)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
      axis.line = element_line(colour = "black"), plot.title=element_text(size=14,face="bold"),
      axis.text=element_text(face="bold",size=12),axis.title=element_text(size=12,face="bold"))+


##################### combine plots
library(gridExtra)
library(Rmisc)

p1 <- ggplot(ae, aes(x=V2)) + geom_histogram(fill="grey", color="black")+
    labs(title="AGE,SEX,chip,PCs,pop-density\nbirth-location,Townsend (pval=0.09)",x="number of pheno with pvalue<0.05")+
    geom_vline(aes(xintercept=30),linetype="dashed")
p2 <- ggplot(c4x, aes(x=V2)) + geom_histogram(fill="grey", color="black")+
    labs(title="AGE,SEX,chip,PCs,Townsend (pval=0.08)",x="number of pheno with pvalue<0.05")+
    geom_vline(aes(xintercept=32),linetype="dashed")
p3 <- ggplot(c5x, aes(x=V2)) + geom_histogram(fill="grey", color="black")+
    labs(title="AGE,SEX,chip,PCs (pval=0.08)",x="number of pheno with pvalue<0.05")+
    geom_vline(aes(xintercept=31),linetype="dashed")

pdf("my_phewas/coding_permut.pdf",12,4)
grid.arrange(p1,p2,p3, ncol=3)
dev.off()


##################### Manhattan plots
library("CMplot")
library(qqman)

png('qqplot.png')
qq(a$P)
dev.off()

CMplot(ss,type="p",plot.type="m",LOG10=TRUE,threshold=c(5E-8),threshold.lty=1,file="jpg",memo="",dpi=300,
    file.output=TRUE,verbose=TRUE,width=14,height=6,chr.labels.angle=45)

############################################################################

pdf("/sc/arion/projects/GENECAD/hamy/ukbb/Figure-S2.pdf",width=10, height=12)
layout(mat = matrix(c(1,4,7,10,2,5,8,11,3,6,9,12), nrow = 4, ncol = 3), heights = c(3,3,3,3),widths = c(5,4,4))
par(mar=c(5, 10, 2, 2))
plot(c(0,3.2), c(0,3.2), col="red", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,3.2), ylim=c(0,4.4), las=1, xaxs="i", yaxs="i", bty="l",main="load score (genome-wide)")
points(load[[4]], load[[1]], pch=20, cex=.4, bg="black") 
text(-1.5, 2, "A", xpd=NA,cex=2.5,font=2)
par(mar=c(5,5,2,2))
plot(c(0,3.2), c(0,3.2), col="red", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,3.2), ylim=c(0,4.4), las=1, xaxs="i", yaxs="i", bty="l",main="load score (non-coding)")
points(load[[5]], load[[2]], pch=20, cex=.4, bg="black") 
plot(c(0,3.2), c(0,3.2), col="red", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,3.2), ylim=c(0,9.3), las=1, xaxs="i", yaxs="i", bty="l",main="load score (coding)")
points(load[[6]], load[[3]], pch=20, cex=.4, bg="black") 

par(mar=c(5, 10, 2, 2))
plot(c(0,3.2), c(0,3.2), col="red", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,3.2), ylim=c(0,4.4), las=1, xaxs="i", yaxs="i", bty="l", main="burden score (genome-wide)")
points(burden[[4]], burden[[1]], pch=20, cex=.4, bg="black") 
text(-1.5, 2, "B", xpd=NA,cex=2.5,font=2)
par(mar=c(5,5,2,2))
for(i in 2:3){
  plot(c(0,3.2), c(0,3.2), col="red", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,3.2), ylim=c(0,4.4), las=1, xaxs="i", yaxs="i", bty="l",main=paste0("burden score (",scores2[i],")"))
  points(burden[[i+3]], burden[[i]], pch=20, cex=.4, bg="black") }

par(mar=c(5, 10, 2, 2))
plot(c(0,3.2), c(0,3.2), col="red", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,3.2), ylim=c(0,4.2), las=1, xaxs="i", yaxs="i", bty="l",main="burden score\n(genome-wide,rare variant)")
points(xx,brare[[1]], pch=20, cex=.4, bg="black") 
text(-1.5, 2, "C", xpd=NA,cex=2.5,font=2)
par(mar=c(5,5,2,2))
for(i in 2:3){
  plot(c(0,3.2), c(0,3.2), col="red", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,3.2), ylim=c(0,4.2), las=1, xaxs="i", yaxs="i", bty="l",main=paste0("burden score\n(",scores2[i],",rare variant,)") )
  points(xx,brare[[i]], pch=20, cex=.4, bg="black") }

par(mar=c(5, 10, 2, 2))
plot(c(0,3.2), c(0,3.2), col="red", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,3.2), ylim=c(0,4.2), las=1, xaxs="i", yaxs="i", bty="l", main="burden score\n(genome-wide,common variant)")
points(xx,bcommon[[i]], pch=20, cex=.4, bg="black") 
text(-1.5, 2, "D", xpd=NA,cex=2.5,font=2)
par(mar=c(5,5,2,2))
for(i in 2:3){
  plot(c(0,3.2), c(0,3.2), col="red", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,3.2), ylim=c(0,4.2), las=1, xaxs="i", yaxs="i", bty="l",main=paste0("burden score\n(",scores2[i],",common variant)") )
  points(xx,bcommon[[i]], pch=20, cex=.4, bg="black") }

dev.off()


