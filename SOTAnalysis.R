library(tidyverse)
library(lme4) #linear mixed model
library(rstatix) #chi-squared tests for mixed model
library(effectsize) # effect size for mixed model

##SOT Analysis##
SOTDataPath="/cloud/project/SOT32_21Fall.csv"
SOTInfoPath="/cloud/project/SOTTrialInfo.csv"


SOT_df <-read.csv(SOTDataPath)
SOT_df$Subject <- as.factor(SOT_df$Subject)
SOT_df$Trial <- as.factor(SOT_df$Trial)

SOTInfo <-read.csv(SOTInfoPath)
SOTInfo$PS180 <- ifelse(SOTInfo$PerspectiveShift<180
                        ,SOTInfo$PerspectiveShift,
                        360-SOTInfo$PerspectiveShift)
SOTInfo$PS90 <- ifelse(SOTInfo$PS180<90
                        ,SOTInfo$PS180,
                        180-SOTInfo$PS180)
SOTInfo$PS45 <- ifelse(SOTInfo$PS90<45
                       ,SOTInfo$PS90,
                       90-SOTInfo$PS90)
SOTInfo$PDir180 <- ifelse(SOTInfo$PointingAngle<180,
                          SOTInfo$PointingAngle,
                          360-SOTInfo$PointingAngle)
SOTInfo$PDir90 <- ifelse(SOTInfo$PDir180<90,
                         SOTInfo$PDir180,
                         180-SOTInfo$PDir180)
SOTInfo$PDir45 <- ifelse(SOTInfo$PDir90<45,
                         SOTInfo$PDir90,
                         90-SOTInfo$PDir90)

#SOTInfo$PS <- cut(SOTInfo$PS180, breaks = c(-0.01, 45, 90, 135,180))
#SOTInfo$PQ <- ifelse(SOTInfo$PointingQuadrant == 'B','B',ifelse(SOTInfo$PointingQuadrant == 'F','F','LR'))
#SOTInfo$PQ <-factor(SOTInfo$PQ,levels=c('F','LR','B'))

SOTtrial <- SOT_df %>%
  group_by(Trial)%>%
  dplyr::summarise(
    n=n(),
    stand = unique(Center),
    face =  unique(Top),
    target = unique(Target),
    ang_error = mean(AngularError),
    se = sd(AngularError)/sqrt(n)
  )

SOT_sum <- merge(SOTtrial,SOTInfo,by.x = "Trial", by.y = "TrialID")%>%
  select(Trial,ang_error,se,PS180,PS90,PS45,PDir180,PDir90,PDir45)

psych::describe(SOT_sum)

chart.Correlation(SOT_sum[,c("ang_error","PS180","PS90","PS45","PDir180","PDir90","PDir45")], histogram=TRUE, pch=19)

ggplot(SOT_sum,aes(x=PS180,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Perspective Shift",y="Angular Error")+
  ylim(0,35)+
  xlim(0,180)+
  theme_bw(base_size=25)

ggplot(SOT_sum,aes(x=PS90,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Perspective Shift(90)",y="Angular Error")+
  ylim(0,35)+
  xlim(0,90)+
  theme_bw(base_size=30)

ggplot(SOT_sum,aes(x=PS45,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Perspective Shift(45)",y="Angular Error")+
  ylim(0,45)+
  xlim(0,45)+
  theme_bw(base_size=30)

ggplot(SOT_sum,aes(x=PDir45,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Pointing Direction(45)",y="Angular Error")+
  ylim(0,45)+
  xlim(0,45)+
  theme_bw(base_size=30)

ggplot(SOT_sum,aes(x=PDir180,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Pointing Direction",y="Angular Error")+
  ylim(0,45)+
  xlim(0,180)+
  theme_bw(base_size=25)

ggplot(SOT_sum,aes(x=PDir90,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Pointing Direction (90)",y="Angular Error")+
  ylim(0,35)+
  xlim(0,90)+
  theme_bw(base_size=25)

names(SOT_sum)


sot_mod_data <- SOT_sum%>%
  mutate(PS180_s = scale(PS180,center=T,scale=T))%>%
  mutate(PS45_s = scale(PS45,center = T, scale = T))%>%
  mutate(PD45_s = scale(PDir45,center = T, scale = T))%>%
  mutate(PD90_s = scale(PDir90,center = T, scale = T))

SOT.lm1 <- lm(ang_error~PS45_s+PD45_s+PS180_s+PS180_PD45_s,data=sot_mod_data)
summary(SOT.lm1)
confint(SOT.lm1)
effectsize::eta_squared(SOT.lm1)



SOTtrial_ps <- SOT_sum%>%
  group_by(PS)%>%
  dplyr::summarise(
    n=n(),
    ps_anchor = mean(PS180),
    ang_m = mean(ang_error),
    ang_se = sd(ang_error)/sqrt(n)
  )

ggplot(SOTtrial_ps ,aes(x=ps_anchor,y=ang_m))+
  geom_point(size=3)+
  geom_line()+
  geom_errorbar(aes(ymin=ang_m-ang_se, ymax=ang_m+ang_se),position=position_dodge(0),width=5)+
  ylim(0,35)+
  ylab('Angular Error')+
  scale_x_continuous(name ="Perspective Shift",breaks=c(22.5,67.5,112.5,157.5),limits=c(0,180),labels = c('[0,45]','(45,90]','(90,135]','(135,180]'))+
  theme_classic(base_size = 20)

SOTtrial_pq <- SOT_sum%>%
  group_by(PQ)%>%
  dplyr::summarise(
    n=n(),
    ang_m = mean(ang_error),
    ang_se = sd(ang_error)/sqrt(n)
  )

ggplot(SOTtrial_pq ,aes(x=PQ,y=ang_m))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=ang_m-ang_se, ymax=ang_m+ang_se),position=position_dodge(0),width=0.3)+
  ylim(0,35)+
  ylab('Angular Error')+
  xlab('Pointing Quadrant')+
  theme_classic(base_size = 20)

## interaction

SOT_combine <- merge(SOT_df[,c("Subject","Trial","AngularError")],SOT_sum[,c("Trial","PerspectiveShift","PointingAngle","PointingQuadrant","PS180","PS","PQ")],by="Trial")
SOT.lmm <- lmer(AngularError~ 
                     PQ*PS180+
                     (1|Subject)+(1|Trial), 
                   data = SOT_combine, 
                   REML = F)

summary(SOT.lmm)
effectsize::eta_squared(SOT.lmm)
Anova(SOT.lmm)

SOT_sum$PS_right <- SOT_sum$PS180 < 90

SOTtrial_pqps <- SOT_sum %>% 
  group_by(PQ,PS_right)%>%
  summarise(
    n = n(),
    ang_m = mean(ang_error),
    ang_se = sd(ang_error)/sqrt(n)
  )


ggplot(SOTtrial_pqps,aes(x=PQ,y=ang_m,color=PS_right,shape=PS_right))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=ang_m-ang_se, ymax=ang_m+ang_se),position=position_dodge(0))+
  ylim(0,45)+
  labs(x= 'Pointing Quadrant', 
       y = 'SOT anglular error',
       color = "Perspective Shift",
       shape = "Perspective Shift")+
  theme_classic(base_size = 20)+
  scale_color_manual(values=c('#999999','#E69F00'),labels = c(">90", "<90"))+
  scale_shape_discrete(labels = c(">90", "<90"))
