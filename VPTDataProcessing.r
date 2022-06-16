## ---------------------------
##
## Script name: VPTDataProcessing
##
## Purpose of script:
##
## Author: Chuanxiuyue (Carol) He
##         
##
## Date Created: 2021-8-23
##
## Email: carol.hcxy@gmail.com 
##        scottmatsubara@gmail.com
##
## ---------------------------


## ---------------------------

## set working directory


thePath = "/cloud/project/iVTTData/"
trialPath="/cloud/project/ViewPointTransformationTrialInfo.csv"
SOTPath="/cloud/project/SOT32_21Fall.csv"
QPath="/cloud/project/iVTT-Exp1_Qualtrics.csv"
## ---------------------------

## load up the packages
## use install.packages() to install new packages
library(plyr)
library(dplyr)
library(tidyverse)
library(reshape2)
library(lme4) #linear mixed model
library(rstatix) #chi-squared tests for mixed model
library(effectsize) # effect size for mixed model

## ---------------------------

#####----------Load Experiment Data----------#####
## read csv files
csv_files_ls <-  list.files(
  path=thePath, 
  pattern = "*.csv")

csv_files_df <- lapply(csv_files_ls, 
                       function(x) 
                         { 
                         tmp <- try(read.csv(
                           file = x, 
                           header = TRUE))
                         if (!inherits(tmp, 
                                       'try-error')) 
                           tmp
                         }
                       )

## combine all csv files together
combined_df <- do.call("rbind", 
                       lapply(csv_files_df, 
                              as.data.frame)
                       )

#read in the trial data
trials_df <-read.csv(trialPath)

#####----------Experiment Data Cleaning----------#####

#Get the last location of the player when they press the trigger and note the start time for each trial
output_data<-combined_df %>% 
  group_by(Participant.ID, Level) %>%
  mutate(
    startTime=first(Time)
  ) %>%
  summarise_all(last)

#remove NA rows
output_data <- na.omit(output_data)

#find distance traveled from the origin
output_data$PosX<-as.numeric(output_data$PosX)
output_data$PosZ<-as.numeric(output_data$PosZ)
output_data$distTraveled <- sqrt((output_data$PosX)^2+(output_data$PosZ)^2)


#find elapsed time
output_data$Time<-as.numeric(output_data$Time)
output_data$startTime<-as.numeric(output_data$startTime)
output_data$elapsedTime<-output_data$Time-output_data$startTime

#rename columns for consistency
names(output_data)[1] <- "ID"
names(output_data)[8] <- "FacingAngle"

#merge experiment output data with the trial data
names(trials_df)[1] <- "Level"
output_data <- merge(output_data,trials_df,by="Level")

#find the distance from the target
output_data$distFromTarget <- sqrt((output_data$PosX-output_data$TargetX)^2+(output_data$PosZ-output_data$TargetZ)^2)

#find angular error 
output_data$angularError <- acos(((output_data$PosX*output_data$TargetX)+(output_data$PosZ*output_data$TargetZ))
                                 /(sqrt((output_data$PosX^2+output_data$PosZ^2))*sqrt((output_data$TargetX^2+output_data$TargetZ^2))))

#convert radians to degrees
output_data$angularError <- output_data$angularError * 180 / pi

#positive angular error means that the subject went further than they needed to (overshot)
#negative angular error means that they undershot
output_data$yCrossProduct<-(output_data$PosZ*output_data$TargetX)-(output_data$PosX*output_data$TargetZ)

#make the angularError signed
output_data$angularError <- output_data$angularError * sign(output_data$yCrossProduct)

#abs travel distance without direction
output_data$distError <- output_data$distTraveled - output_data$TravelDistance

#select useful columns
output_data<-output_data %>%
  dplyr::select(ID,Level,gender,PosX,PosZ,TargetX,TargetZ,PerspectiveShift,PointingAngle,TravelDistance,FacingAngle,elapsedTime,distTraveled,distFromTarget,angularError,distError)

#remove practice trials
output_data <- output_data[!(output_data$Level %in% c(-1,0)),]
#remove invalid participants
#output_data <- output_data[!(output_data$ID %in% c(4)),]

#results by trials
trial_error <- output_data %>%
  group_by(Level)%>%
  dplyr::summarise(
    n=n(),
    PS = unique(PerspectiveShift),
    PS_180 = ifelse(PS<180, PS,360-PS),
    PS_90 = ifelse(PS_180<90, PS_180,180-PS_180),
    PS_45 = ifelse(PS_90<45, PS_90,90-PS_90),
    TDis = unique(TravelDistance),
    TDir = unique(PointingAngle),
    TDir_180 = ifelse(TDir<180,TDir,360-TDir),
    TDir_90 = ifelse(TDir_180<90,TDir_180,180-TDir_180),
    TDir_45 = ifelse(TDir_90<45,TDir_90,90-TDir_90),
    fromT_error = mean(distFromTarget),
    signed_ang = mean(angularError),
    signed_dist = mean(distError),
    ang_error = mean(abs(angularError)),
    ang_error_sd = sd(abs(angularError)),
    dist_error = mean(abs(distError)),
    dist_error_sd = sd(abs(distError))
  )

chart.Correlation(trial_error[,c("signed_ang","ang_error","signed_dist","dist_error","PS_180","PS_90","PS_45","TDir_180","TDir_90","TDir_45","TDis")], histogram=TRUE, pch=19)

psych::describe(trial_error[,c("ang_error","signed_ang","dist_error","signed_dist")])

chart.Correlation(trial_error[,c("ang_error","dist_error","signed_dist","PS_180","PS_90","PS_45","TDir_180","TDir_90","TDir_45","TDis")], histogram=TRUE, pch=19)
#write.csv(trial_error,file = "iVPT-Trial.csv",row.names=FALSE)

#visualize results
hist(trial_error$dist_error)
hist(trial_error$ang_error)
hist(trial_error$fromT_error)
plot(trial_error$dist_error,trial_error$ang_error)

cor.test(trial_error$dist_error,trial_error$ang_error)

#descriptive stats
psych::describe(trial_error$dist_error)
psych::describe(trial_error$ang_error)

#trial_error$PQ <- factor(trial_error$PQ, levels = c("F","LR","B"))

#trial error by features
##perspective shift
#plot(trial_error$PS,trial_error$ang_error)

plot(trial_error$PS_180,trial_error$ang_error)
plot(trial_error$TDir_180,trial_error$ang_error)
plot(trial_error$TDir_90,trial_error$ang_error)
plot(trial_error$TDir_45,trial_error$ang_error)
cor.test(trial_error$PS_180,trial_error$ang_error)
cor.test(trial_error$TDir_45,trial_error$ang_error)
hist(trial_error$ang_error)

ggplot(trial_error,aes(x=PS_180,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Perspective Shift",y="Angular Error")+
  ylim(0,45)+
  xlim(0,180)+
  theme_bw(base_size=30)

ggplot(trial_error,aes(x=TDir_180,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Travel Direction",y="Angular Error")+
  ylim(0,45)+
  xlim(0,180)+
  theme_bw(base_size=30)

ggplot(trial_error,aes(x=TDir_90,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Travel Direction (90)",y="Angular Error")+
  ylim(0,45)+
  xlim(0,90)+
  theme_bw(base_size=30)

ggplot(trial_error,aes(x=TDir_45,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Travel Direction (45)",y="Angular Error")+
  ylim(0,45)+
  xlim(0,45)+
  theme_bw(base_size=30)

ggplot(trial_error,aes(x=TDis,y=ang_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Travel Distance",y="Angular Error")+
  ylim(0,45)+
  xlim(1,3)+
  theme_bw(base_size=30)

psych::describe(trial_error)


library(plotly)
# install.packages("plotly")
fig <- plot_ly(trial_error, x = ~PS_180, y = ~TDir_180, z = ~ang_error)
fig <- fig %>% add_markers(opacity = 0.5)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Perspective Shift'),
                                   yaxis = list(title = 'Travel Direction'),
                                   zaxis = list(title = 'Abs Angluar Error')))


trial_mod_data <- trial_error%>%
  mutate(PS_s = scale(PS_180,center = T, scale = T))%>%
  mutate(T180_s = scale(TDir_180,center = T, scale = T))%>%
  mutate(T90_s = scale(TDir_90,center=T,scale=T))%>%
  mutate(T45_s = scale(TDir_45,center=T,scale=T))%>%
  mutate(TD_s = scale(TDis,center=T,scale=T))

trial.lm1 <- lm(ang_error~PS_s*T45_s+T90_s,data=trial_mod_data)
summary(trial.lm1)
confint(trial.lm1)
effectsize::eta_squared(trial.lm1)


ggplot(trial_error,aes(x=PS_180,y=dist_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Perspective Shift",y="Distance Error")+
  xlim(0,180)+
  ylim(0.15,0.45)+
  theme_bw(base_size=30)

ggplot(trial_error,aes(x=PS_90,y=dist_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Perspective Shift",y="Distance Error")+
  xlim(0,90)+
  ylim(0.15,0.45)+
  theme_bw(base_size=30)


ggplot(trial_error,aes(x=TDir_180,y=dist_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Travel Direction",y="Distance Error")+
  xlim(0,180)+
  ylim(0.15,0.45)+
  theme_bw(base_size=30)

ggplot(trial_error,aes(x=TDir_45,y=dist_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Travel Direction (45)",y="Distance Error")+
  xlim(0,45)+
  ylim(0.15,0.45)+
  theme_bw(base_size=30)

ggplot(trial_error,aes(x=TDis,y=dist_error))+
  geom_point(size=5,alpha=0.5)+
  labs(x="Travel Distance",y="Distance Error")+
  ylim(0.15,0.45)+
  xlim(1.4,2.9)+
  theme_bw(base_size=30)

trial.lm2 <- lm(dist_error~PS_s+T180_s+TD_s,data=trial_mod_data)
summary(trial.lm2)
confint(trial.lm2)
effectsize::eta_squared(trial.lm2)

trial.lm2b <- lm(dist_error~PS_s*T45_s+T180_s,data=trial_mod_data)
summary(trial.lm2b)
confint(trial.lm2b)
effectsize::eta_squared(trial.lm2b)



#plot(trial_error$PS,trial_error$dist_error)
#plot(trial_error$PS_180,trial_error$dist_error)
#cor.test(trial_error$PS_180,trial_error$dist_error)

##pointing Quatriant

# trial_qua <- trial_error %>% 
#   group_by(PQ)%>%
#   dplyr::summarise(
#     n = n(),
#     ang_m = mean(ang_error),
#     ang_se = sd(ang_error)/sqrt(n),
#     dist_m = mean(dist_error),
#     dis_se = sd(dist_error)/sqrt(n))
  
# ggplot(trial_error, aes(x = PQ, y = ang_error)) +
#   geom_point() +
#   geom_bar(data = trial_qua, aes(y = ang_m),stat = "identity", alpha = .3)+
#   ylim(0,45)

# trial_quaInter <- trial_error %>% 
#   group_by(PQ,PS_right)%>%
#   dplyr::summarise(
#     n = n(),
#     ang_m = mean(ang_error),
#     ang_se = sd(ang_error)/sqrt(n),
#     dist_m = mean(dist_error),
#     dis_se = sd(dist_error)/sqrt(n))
# 
# ggplot(trial_quaInter,aes(x=PQ,y=ang_m,color=PS_right,shape=PS_right))+
#   geom_point(size=3)+
#   geom_errorbar(aes(ymin=ang_m-ang_se, ymax=ang_m+ang_se),position=position_dodge(0))+
#   ylim(0,45)+
#   labs(x= 'Pointing Quadrant', 
#        y = 'iVPT anglular error',
#        color = "Perspective Shift",
#        shape = "Perspective Shift")+
#   theme_classic(base_size = 20)+
#   scale_color_manual(values=c('#999999','#E69F00'),labels = c(">90", "<90"))+
#   scale_shape_discrete(labels = c(">90", "<90"))

# trial_disInter <- trial_error %>% 
#   group_by(PQ,TD_S)%>%
#   dplyr::summarise(
#     n = n(),
#     ang_m = mean(ang_error),
#     ang_se = sd(ang_error)/sqrt(n),
#     dist_m = mean(dist_error),
#     dis_se = sd(dist_error)/sqrt(n))

# ggplot(trial_disInter,aes(x=PQ,y=dist_m,color=TD_S,shape=TD_S))+
#   geom_point(size=3)+
#   geom_errorbar(aes(ymin=dist_m-dis_se, ymax=dist_m+dis_se),position=position_dodge(0))+
#   ylim(0.1,0.5)+
#   labs(x= 'Pointing Quadrant', 
#        y = 'iVPT distance error',
#        color = "Perspective Shift",
#        shape = "Perspective Shift")+
#   theme_classic(base_size = 20)+
#   scale_color_manual(values=c('#999999','#E69F00'),labels = c(">90", "<90"))+
#   scale_shape_discrete(labels = c(">90", "<90"))

output_data$Level <- as.factor(output_data$Level)
#output_data$PointingQuadrant <- as.factor(output_data$PointingQuadrant)
output_data$ID <- as.factor(output_data$ID)


model_data <- output_data
model_data$PerspectiveShift <- ifelse(output_data$PerspectiveShift<=180,output_data$PerspectiveShift,360-output_data$PerspectiveShift)
model_data$PointingAngle <- ifelse(output_data$PointingAngle<180,output_data$PointingAngle,360-output_data$PointingAngle)
model_data$PointingAngle_90 <- ifelse(model_data$PointingAngle<90,model_data$PointingAngle,180-model_data$PointingAngle)
model_data$PointingAngle_45 <- ifelse(model_data$PointingAngle_90<45,model_data$PointingAngle_90,90-model_data$PointingAngle_45)
model_data$AbsAngError <- abs(model_data$angularError)
model_data$AbsDistError <- abs(model_data$distError)

library(MASS)
fitdistr(model_data$AbsAngError, "gamma")
ks.test(model_data$AbsAngError,0.998045842,0.045048728)

model_data <- model_data%>%
  mutate(PS_s = scale(PerspectiveShift,center = T,scale=T)) %>%
  mutate(PA_s = scale(PointingAngle_45,center=T, scale=T))%>%
  mutate(TD_s = scale(TravelDistance,center = T, scale=T))%>%
  mutate(AbsAngError_trans = log(AbsAngError+1))

hist(model_data$AbsAngError_trans)
hist(model_data$PA_s)


library(lme4) #linear mixed model
library(rstatix) #chi-squared tests for mixed model
trial.lmm1 <- lmer(AbsAngError_trans ~ 
                    PS_s+PA_s+
                    (1|ID)+(1|Level), 
                  data = model_data, 
                  REML = F)

summary(trial.lmm1)
effectsize::eta_squared(trial.lmm1)
Anova(trial.lmm1)
confint(trial.lmm1)


# ggplot(trial_error, aes(x = PointingQuadrant, y = dist_error)) +
#   geom_point() +
#   geom_bar(data = trial_qua,aes(y=dist_m), stat = "identity", alpha = .3)

trial.lmm2 <- lmer(AbsDistError ~ 
                     PointingAngle+PerspectiveShift+TravelDistance+
                     (1|ID)+(1|Level), 
                   data = model_data, 
                   REML = F)

summary(trial.lmm2)
Anova(trial.lmm2)
confint(trial.lmm2)
effectsize::eta_squared(trial.lmm2)


#bayesian linear model

install.packages('brms')
install.packages('bayestestR')
install.packages('emmeans')

library(brms)
library(bayestestR)
library(emmeans)

set.seed(123)

library(MASS)
fitdistr(model_data$AbsAngError, "gamma")
ks.test(model_data$AbsAngError,0.998045842,0.045048728)

fit1 <- brm(formula = AbsAngError ~ PerspectiveShift*PointingAngle + (1 | ID) + (1 | Level), 
            data = model_data, family = Gamma(link = "log"),
            chains=2,iter=1000)

brm_glm_reg <- brm(y~x, data=fake_data_a, family=Gamma(link="log"),
                   prior=c(prior(normal(0,2),class="Intercept"),
                           prior(normal(0,2),class="b"),
                           prior(gamma(0.01,0.01),class="shape")),
                   chains=2,iter=1000, cores=4)



##reliability
##install.packages("splithalf")
library("splithalf")

output_data['abs_ang_error'] <- abs(output_data$angularError)
output_data['log_ang_error'] <- log(output_data$abs_ang_error)
output_data['abs_dist_error'] <- abs(output_data$distError)


trial_45 <- c('1','2','3','4','33','34','29','20','31','32','47','48')
trial_90 <- c('5','6','7','8','35','36','25','26','27','28','45','46')

data_36 <- output_data[!output_data$Level%in%trial_45,]
data_24 <- data_36[!data_36$Level%in%trial_90,]

set.seed(123)
splithalf(data=output_data,
          outcome = "accuracy",
          score = "average",
          halftype = "random",
          permutations = 5000,
          var.ACC = "abs_ang_error",
          var.trialnum = "Level",
          var.participant = "ID",
          average="mean")

splithalf(data=data_36,
          outcome = "accuracy",
          score = "average",
          halftype = "random",
          permutations = 5000,
          var.ACC = "abs_ang_error",
          var.trialnum = "Level",
          var.participant = "ID",
          average="mean")

splithalf(data=data_24,
          outcome = "accuracy",
          score = "average",
          halftype = "random",
          permutations = 5000,
          var.ACC = "abs_ang_error",
          var.trialnum = "Level",
          var.participant = "ID",
          average="mean")

splithalf(data=output_data,
          outcome = "accuracy",
          score = "average",
          halftype = "random",
          permutations = 5000,
          var.ACC = "log_ang_error",
          var.trialnum = "Level",
          var.participant = "ID",
          average="mean")

splithalf(data=output_data,
          outcome = "accuracy",
          score = "average",
          halftype = "random",
          permutations = 5000,
          var.ACC = "abs_dist_error",
          var.trialnum = "Level",
          var.participant = "ID",
          average="mean")

#performance by participants
ID_error <- output_data %>%
  group_by(ID)%>%
  dplyr::summarise(
    gender = unique(gender),
    n=n(),
    fromT_error = mean(distFromTarget),
    ang_error = mean(abs(angularError)),
    dist_error = mean(abs(distError)),
    log_error = mean(log_ang_error)
  )

#visualize results
hist(ID_error$dist_error)
hist(ID_error$ang_error)
hist(ID_error$fromT_error)
plot(ID_error$dist_error,ID_error$ang_error)

SOT_df <-read.csv(SOTPath)
SOT_df$Subject <- as.factor(SOT_df$Subject)
SOT_df$Trial <- as.factor(SOT_df$Trial)

splithalf(data=SOT_df,
          outcome = "accuracy",
          score = "average",
          halftype = "random",
          permutations = 5000,
          var.ACC = "AngularError",
          var.trialnum = "Trial",
          var.participant = "Subject",
          average="mean")

SOT_error <- SOT_df%>%
  group_by(Subject)%>%
  dplyr::summarise(
    age = unique(Age),
    n=n(),
    sot_error = mean(AngularError),
    sot_sd = sd(AngularError)
  )

sub <- merge(ID_error,SOT_error,by.x="ID",by.y="Subject",all.x=TRUE)

plot(sub$ang_error,sub$sot_error)
plot(sub$log_error,sub$sot_error)
cor.test(sub$log_error,sub$sot_error)
0.6881631/sqrt(0.85*0.80)

ggplot(sub,aes(x=ang_error,y=sot_error))+
  geom_point(size=5,alpha=0.5)+
  geom_abline(intercept = 0, slope = 1,color='red')+
  labs(x="iVTT Angular Error",y="SOT Angular Error")+
  xlim(5,80)+
  ylim(5,80)+
  theme_bw(base_size=22)

t.test(sub$ang_error,sub$sot_error,paired=TRUE)

#Qualtrics
Q_df <-read.csv(QPath)
Q_df$Subject <- as.factor(Q_df$Subject)
Q_sum <- Q_df%>%
  mutate(SA = rowMeans(.[grepl('SA', colnames(Q_df))]))%>%
  mutate(GPS = rowMeans(.[grepl('GPS', colnames(Q_df))]))%>%
  mutate(SBSOD = rowMeans(.[grepl('SBSOD', colnames(Q_df))]))%>%
  select(Subject,Order,MS,VRExp,SA,GPS,SBSOD,Age,Gender)

sub_sum <- merge(sub,Q_sum,by.x="ID",by.y="Subject")%>%
  select(ID,gender,Gender,age,Age,Order,ang_error, log_error,sot_error,dist_error,MS,VRExp,SA,GPS,SBSOD)

psych::describe(sub_sum[,7:14])
#install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
chart.Correlation(sub_sum[,c("log_error","sot_error","dist_error","SBSOD","GPS","SA")], histogram=TRUE, pch=19)

cor.test(sub_sum$sot_error,sub_sum$dist_error)

cor.test(sub_sum$log_error,sub_sum$sot_error)
cor.test(sub_sum$log_error,sub_sum$dist_error)
cor.test(sub_sum$sot_error,sub_sum$SBSOD)
psych::r2d(0.3868849)


psych::cohen.d(sub_sum$log_error,sub_sum$dist_error)


ggplot(data = sub_sum, aes(x=ang_error,y=sot_error))+
  geom_point()+
  ylim(0,50)+
  xlim(0,80)+
  labs(x= 'iVPT Angluar Error', 
       y = 'SOT Anglular Error')+
  theme_classic(base_size = 20)
  

SBSOD_raw <-  Q_df[grepl('SBSOD|Subject',colnames(Q_df))]

SBSOD_long <-  reshape(SBSOD_raw, 
                       varying = list(names(SBSOD_raw)[2:16]),
                       direction = "long",
                       v.names = "SBSOD",
                       idvar = c("Subject"),
                       timevar = "Item",
                       times = 1:15)

splithalf(data=SBSOD_long,
          outcome = "accuracy",
          score = "average",
          halftype = "random",
          permutations = 5000,
          var.ACC = "SBSOD",
          var.trialnum = "Item",
          var.participant = "Subject",
          average="mean")

SA_raw <-  Q_df[grepl('SA|Subject',colnames(Q_df))]

SA_long <-  reshape(SA_raw, 
                       varying = list(names(SA_raw)[2:14]),
                       direction = "long",
                       v.names = "SA",
                       idvar = c("Subject"),
                       timevar = "Item",
                       times = 1:13)

splithalf(data=SA_long,
          outcome = "accuracy",
          score = "average",
          halftype = "random",
          permutations = 5000,
          var.ACC = "SA",
          var.trialnum = "Item",
          var.participant = "Subject",
          average="mean")

GPS_raw <-  Q_df[grepl('GPS|Subject',colnames(Q_df))]

GPS_long <-  reshape(GPS_raw, 
                       varying = list(names(GPS_raw)[2:9]),
                       direction = "long",
                       v.names = "GPS",
                       idvar = c("Subject"),
                       timevar = "Item",
                       times = 1:8)

splithalf(data=GPS_long,
          outcome = "accuracy",
          score = "average",
          halftype = "random",
          permutations = 5000,
          var.ACC = "GPS",
          var.trialnum = "Item",
          var.participant = "Subject",
          average="mean")

0.69/sqrt(0.85*0.80)
0.39/sqrt(0.85*0.91)
0.13/sqrt(0.85*0.89)
0.17/sqrt(0.85*0.86)
0.12/sqrt(0.85*0.75)

0.1/sqrt(0.80*0.91)
0.32/sqrt(0.80*0.89)
0.22/sqrt(0.80*0.86)
0.21/sqrt(0.80*0.75)

0.16/sqrt(0.91*0.89)
0.07/sqrt(0.91*0.86)
0.11/sqrt(0.91*0.75)

sub_test <- sub_sum

sub_test <- sub_test%>%
  mutate(log_s = scale(log_error,center = T,scale = T))%>%
  mutate(sot_s = scale(sot_error,center = T,scale = T))

trial.lm5 <- lm(SBSOD~log_s*sot_s,data=sub_test)
summary(trial.lm5)
confint(trial.lm5)
effectsize::eta_squared(trial.lm5)

trial.lm5 <- lm(SA~log_s*sot_s,data=sub_test)
summary(trial.lm5)
confint(trial.lm5)
effectsize::eta_squared(trial.lm5)

trial.lm5 <- lm(GPS~log_s*sot_s,data=sub_test)
summary(trial.lm5)
confint(trial.lm5)
effectsize::eta_squared(trial.lm5)

trial.lm5 <- lm(log_error~log(VRExp),data=sub_test)
summary(trial.lm5)
confint(trial.lm5)
effectsize::eta_squared(trial.lm5)

hist(log(sub_test$VRExp))
hist(sub_test$VRExp)

View(Q_df[Q_df$VRExp>1,])
sub_test[sub_test$VRExp>1,]
psych::describe(sub_test$VRExp)

t.test(log_error~Order,data = sub_sum, paired = FALSE)
t.test(sot_error~Order,data = sub_sum, paired = FALSE)


#Todo: sex differences
t.test(log_error~Gender,data=sub_test,paired=FALSE)
t.test(sot_error~Gender,data=sub_test,paired=FALSE)


#round the data values to 2 decimal places
output_data <- output_data %>% mutate_if(is.numeric, ~round(., 2))

psych::describe(trial_mod_data)
library(PerformanceAnalytics)
chart.Correlation(trial_mod_data[,c("ang_error","dist_error","PS_180","TDir_180","TDir_90","TDir_45","TDis")], histogram=TRUE, pch=19)



