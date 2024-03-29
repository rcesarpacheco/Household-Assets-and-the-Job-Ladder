library(here)
library(ggplot2)
library(data.table)
library(ggrepel)
library(RColorBrewer)
library(kableExtra)
library(latex2exp)
library(forcats)
library(ggthemes)
library(cowplot)
library(dplyr)
library(haven)
library(Hmisc)
library(kableExtra)
library(cowplot)




# histogram delta wage introduction ---------------------------------------
data <- data.table(read_dta(here('Dropbox/micro_data/firm_wage_assets_clean.dta')))

data <- data[order(id,interview_year),]

data[,delta_w:=lwage-shift(lwage),by=id]
data[,wage_shift:=shift(lwage),by=id]

data[,same_job:=(firm_num==shift(firm_num)),by=id]
data[,same_occ:=(occup==shift(occup)),by=id]
data[,consec_years:=(interview_year==shift(interview_year)+1),by=id]
data_hist <- data[consec_years & same_job & same_occ,]

#comparing with madeline old values
nrow(data_hist[delta_w <= 0.05,])/nrow(data_hist)
# plot

plt <- ggplot(data_hist, aes(x=delta_w)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=0.05, fill="#007bc8", color="#e9ecef", alpha=0.8) +
  # ggtitle('Change in Real Wages for Workers Holding the Same Job') +
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlim(c(-1,1))+
  scale_y_continuous(labels = scales::percent)+
  xlab(TeX("$\\Delta \\ln(wage)$"))+
  ylab('Frequency')
plt

ggsave(here('Results/Figures//hist_change_wages.png'),plot = plt,height = 5,width=8,scale=1)


factorx <- factor(cut(data_hist$delta_w, breaks=c(-10,-5,0,5,10)/100))
#Tabulate and turn into data.frame
xout <- as.data.frame(table(factorx))
#Add cumFreq and proportions
xout <- transform(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))

quantile(data_hist$delta_w)

breaks <- c(-0.10,-0.05,0,0.05,0.10)
nrow(data_hist[delta_w< -0.10])/nrow(data_hist)*100
nrow(data_hist[delta_w> 0.10])/nrow(data_hist)*100

# histogram delta wages residuals introduction --------------------------------------

data <- data.table(read_dta(here('Dropbox/micro_data/firm_wage_assets_clean.dta')))

data <- data[order(id,interview_year),]

data[,delta_w:=lwage_resid-shift(lwage_resid),by=id]
data[,wage_shift:=shift(lwage_resid),by=id]

data[,same_job:=(firm_num==shift(firm_num)),by=id]
data[,same_occ:=(occup==shift(occup)),by=id]
data[,consec_years:=(interview_year==shift(interview_year)+1),by=id]
data_hist <- data[consec_years & same_job & same_occ,]

#comparing with madeline old values
nrow(data_hist[delta_w <= 0.05,])/nrow(data_hist)
# plot

plt <- ggplot(data_hist, aes(x=delta_w)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=0.05, fill="#007bc8", color="#e9ecef", alpha=0.8) +
  # ggtitle('Change in Real Wages for Workers Holding the Same Job') +
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlim(c(-1,1))+
  scale_y_continuous(labels = scales::percent)+
  xlab(TeX("$\\Delta \\ln(wage)$"))+
  ylab('Frequency')
plt

ggsave(here('Results/Figures/hist_change_wages_resids.png'),plot = plt,height = 5,width=8,scale=1)


factorx <- factor(cut(data_hist$delta_w, breaks=c(-10,-5,0,5,10)/100))
#Tabulate and turn into data.frame
xout <- as.data.frame(table(factorx))
#Add cumFreq and proportions
xout <- transform(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))

quantile(data_hist$delta_w)

breaks <- c(-0.10,-0.05,0,0.05,0.10)
nrow(data_hist[delta_w< -0.10])/nrow(data_hist)*100
nrow(data_hist[delta_w> 0.10])/nrow(data_hist)*100




# histogram delta wages residuals introduction --------------------------------------

data <- data.table(read_dta(here('Dropbox/micro_data/firm_wage_assets_clean.dta')))

data <- data[order(id,interview_year),]

data[,delta_w:=lwage_resid-shift(lwage_resid),by=id]
data[,wage_shift:=shift(lwage_resid),by=id]

data[,same_job:=(firm_num==shift(firm_num)),by=id]
data[,same_occ:=(occup==shift(occup)),by=id]
data[,consec_years:=(interview_year==shift(interview_year)+1),by=id]
data_hist <- data[consec_years & same_job & same_occ,]

#comparing with madeline old values
nrow(data_hist[delta_w <= 0.05,])/nrow(data_hist)
# plot

plt <- ggplot(data_hist, aes(x=delta_w)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=0.05, fill="#007bc8", color="#e9ecef", alpha=0.8) +
  # ggtitle('Change in Real Wages for Workers Holding the Same Job') +
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  xlim(c(-1,1))+
  scale_y_continuous(labels = scales::percent)+
  xlab(TeX("$\\Delta \\ln(wage)$"))+
  ylab('Frequency')
plt

ggsave(here('Results/Figures//hist_change_wages.png'),plot = plt,height = 5,width=8,scale=1)


factorx <- factor(cut(data_hist$delta_w, breaks=c(-10,-5,0,5,10)/100))
#Tabulate and turn into data.frame
xout <- as.data.frame(table(factorx))
#Add cumFreq and proportions
xout <- transform(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))

quantile(data_hist$delta_w)

breaks <- c(-0.10,-0.05,0,0.05,0.10)
nrow(data_hist[delta_w< -0.10])/nrow(data_hist)*100
nrow(data_hist[delta_w> 0.10])/nrow(data_hist)*100

# preparing firm output data for matlab Log case --------------------------
data_cluster <- fread(here('Data/firms_logs_sep28.csv'))
tails_cuttof <- 1/100
data_cluster <- rename(data_cluster,"mean_wage"="clust_mean_lwage_resid",
                       "sigma_wage"="clust_std_lwage_resid",
                       "theta"="theta1")
data_cluster[,sigma_ou_matlab:=sqrt(2*theta)*sigma_wage]
data_cluster[,lower_w:=qlnorm(tails_cuttof,mean=mean_wage,sd = sigma_wage )]
data_cluster[,upper_w:=qlnorm(1-tails_cuttof,mean=mean_wage,sd = sigma_wage )]
# data_cluster[,lower_w:=6.4]
# data_cluster[,upper_w:=100]
data_cluster[,min_w:=min(data_cluster$lower_w)]
data_cluster[,max_w:=max(data_cluster$upper_w)]
data_cluster[,firm_share:=firm_size/sum(firm_size)]
fwrite(data_cluster,file = here('Data/firms_out.csv'))


# generate fake dataset ---------------------------------------------------
tails_cuttof <- 2.5/100
n_firms <- 25
mean_wage <- rep(-0.5,n_firms)
sigma_wage <- seq(from=1,to=3,length.out=n_firms)
firm_rank <- 1:n_firms
firm_size <- rep(1,n_firms)
rank_fun <- 0
thetas <- rep(0.15,n_firms)
data_fake <- data.table(firm_cluster=firm_rank,
                        mean_wage=mean_wage,
                        sigma_wage=sigma_wage,
                        firm_size=firm_size,
                        theta=thetas,
                        rank_fun=rank_fun,
                        firm_rank=firm_rank)
data_fake[,firm_share:=firm_size/sum(firm_size)]
data_fake[,sigma_ou_matlab:=sqrt(2*theta)*sigma_wage]
data_fake[,lower_w:=qlnorm(tails_cuttof,mean=mean_wage,sd = sigma_wage )]
data_fake[,upper_w:=qlnorm(1-tails_cuttof,mean=mean_wage,sd = sigma_wage )]
data_fake[,min_w:=min(data_fake$lower_w)]
data_fake[,max_w:=max(data_fake$upper_w)]
fwrite(data_fake,file = here('Data/firms_out_fake_feb_2024.csv'))




# plot sigma exercise -----------------------------------------------------
theta <- 0.15
data_ex_sigmas <- fread(here('Data/firms_out_fake_feb_2024.csv'))
data_plot_sigma <- fread(here('Results/sigma_exercise_output_matlab2.csv'))

data_plot_sigma <- data_plot_sigma[,c(1,2,8,501,1001)]
colnames(data_plot_sigma) <-c('sigma_wage','No assets','Poor','Middle Wealth','Rich')
data_plot_sigma[,sigma_wage:=sigma_wage/sqrt(2*theta)]

data_plot_sigma <- melt(data_plot_sigma,id.vars = 'sigma_wage',variable.name = 'Wealth',value.name = 'Value Job Offer')


plt1 <- ggplot(data_plot_sigma[Wealth=='No assets',],aes(x=sigma_wage,y=`Value Job Offer`))+
  geom_line(color='#007bc8')+
  geom_point(color='#007bc8',size=0.5)+
  scale_y_continuous(labels = scales::label_number(accuracy = 1))+
  theme_classic()+
  ggtitle('No assets')+
  ylab('')+
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust=0.5,size = 9))

plt2 <- ggplot(data_plot_sigma[Wealth=='Poor',],aes(x=sigma_wage,y=`Value Job Offer`))+
  geom_line(color='#007bc8')+
  geom_point(color='#007bc8',size=0.5)+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01))+
  theme_classic()+
  ggtitle('Poor')+
  ylab('')+
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust=0.5,size = 9))

plt3 <- ggplot(data_plot_sigma[Wealth=='Middle Wealth',],aes(x=sigma_wage,y=`Value Job Offer`))+
  geom_line(color='#007bc8')+
  geom_point(color='#007bc8',size=0.5)+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
  ggtitle('Middle Wealth')+
  theme_classic()+
  ylab('')+
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust=0.5,size = 9))


plt4 <- ggplot(data_plot_sigma[Wealth=='Rich',],aes(x=sigma_wage,y=`Value Job Offer`))+
  geom_line(color='#007bc8')+
  geom_point(color='#007bc8',size=0.5)+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
  xlab(TeX(r"($\sigma_f$)"))+
  ggtitle('Rich')+
  ylab('')+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust=0.5,size = 9))

plt1
plt2
plt3

figure <- plot_grid(plt1,plt2,plt3,plt4,
                    ncol = 1, nrow = 4, align = "v")
figure
teste <- grid.arrange(figure,left='Value Job Offer')

ggsave(here('Results/Figures/sigma_exercise2.png'),plot = teste,height = 4,width=4,scale=1.5)




# plot firm characteristic -------------------------------------------------------------------
data_cluster <- fread(here('wage in levels/firms_levels_out.csv'))

plt <- ggplot(data_cluster,aes(x=mean_wage_resid,y=sd_diff_wage_resid,label=firm_rank))+
  geom_point(aes(size=firm_size),color='#007bc8')+
  theme_classic()+
  geom_text_repel(nudge_y=0.5,show.legend=FALSE,min.segment.length = Inf)+
  ylab(TeX(r"($\sigma_f$)"))+
  labs(size="Firm Size")+
  xlab(TeX(r"($\bar {w}_f$)"))
plt
# ggsave(here('Figures/firm_characteristics.png'),plot = plt,height = 5,width=8,scale=1.2)


# indifference curve data -------------------------------------------------

# data_firms <- fread(here('wage in levels/firms_levels_out.csv'))
theta <- 0.15
tails_cuttof <- 2.5/100
# min_sigma <- min(data_firms$sigma_ou_matlab)
# max_sigma <- max(data_firms$sigma_ou_matlab)
# min_mean <- min(data_firms$mean_wage_resid)
# max_mean <- max(data_firms$mean_wage_resid)
min_sigma <- 0.05
max_sigma <- 0.5
min_mean <- -0.35
max_mean <- 0.35
n_sigma <- 10
n_mean <- 10
grid_mean <- seq(from=min_mean,to=max_mean,length.out =n_mean)
grid_sigma <- seq(from=min_sigma,to=max_sigma,length.out =n_sigma)
grid <- expand.grid(grid_mean,grid_sigma)
data_indifference <- data.table(mean_wage=grid$Var1,sigma_ou_matlab=grid$Var2)
data_indifference[,sd_diff_wage_resid:=sigma_ou_matlab/sqrt(2*theta)]
data_indifference[,firm_share:=0]
data_indifference[,rank_fun:=sd_diff_wage_resid-mean_wage]
data_indifference[,firm_rank:=rank(rank_fun)]

data_indifference[,lower_w:=qlnorm(tails_cuttof,mean=mean_wage,sd = sd_diff_wage_resid)]
data_indifference[,upper_w:=qlnorm(1-tails_cuttof,mean=mean_wage,sd = sd_diff_wage_resid)]
data_indifference[,min_w:=min(data_indifference$lower_w)]
data_indifference[,max_w:=max(data_indifference$upper_w)]
data_indifference[,theta:=theta]

fwrite(data_indifference,file = here('Data/data_indifference.csv'))

# indifference curve plots ------------------------------------------------
pallete <- brewer.pal(6,name = 'RdBu')

data_plot_indifference_poor <- fread(here("Results/indifference_curve_poor.csv"),col.names = c('Mean_log_wage','Sd_log_wage','V'))
data_plot_indifference_rich <- fread(here("Results/indifference_curve_rich.csv"),col.names = c('Mean_log_wage','Sd_log_wage','V'))



data_plot_indifference_poor[,'Wealth':='Poor']
data_plot_indifference_rich[,'Wealth':='Rich']

data_plot_indifference <- rbind(data_plot_indifference_poor,data_plot_indifference_rich)
# data_plot_indifference[,V:=round(V,4)]


data_cluster <- fread('Data/cluster_data_monthly.csv')

line_width <- 0.7
plt_poor <- ggplot()+
  # geom_contour(data = data_plot_indifference[Wealth=='Poor',], 
  #              aes(x = Mean_log_wage, y = Sd_log_wage, z = V, colour = after_stat(level)),linewidth=0.2,
  #              breaks = quantile(data_plot_indifference[Wealth=='Poor',V],seq(0,1,length.out=15)))+
  
  geom_contour(data = data_plot_indifference[Wealth=='Poor',], 
               aes(x = Mean_log_wage, y = Sd_log_wage, z = V, colour = after_stat(level)),linewidth=line_width,
               bins = 15)+
  scale_color_gradient(low = pallete[6],high = pallete[1],labels = function(x) sprintf("%.2f", x))+
  theme_bw()+
  ylab(TeX(r"($\sigma_f$)"))+
  xlab(TeX(r"($\bar {w}_f$)"))+
  labs(colour="Value Job\nOffer ")+
  geom_point(data=data_cluster,aes(x=mean_lwage_resid,y=std_lwage_resid))+
  # geom_text_repel(data=data_cluster,aes(x=mean_lwage_resid,y=std_lwage_resid,label=firm_rank),nudge_y=0.01,show.legend=FALSE,min.segment.length = Inf)+
  theme(panel.grid = element_blank(),legend.title.align = 0.5)

plt_rich <- ggplot()+
  geom_contour(data = data_plot_indifference[Wealth=='Rich',], 
               aes(x = Mean_log_wage, y = Sd_log_wage, z = V, colour = after_stat(level)),linewidth=line_width,
               bins = 15)+
  scale_color_gradient(low = pallete[6],high = pallete[1],labels = function(x) sprintf("%.3f", x))+
               theme_bw()+
  ylab(TeX(r"($\sigma_f$)"))+
  xlab(TeX(r"($\bar {w}_f$)"))+
  labs(colour="Value Job\nOffer ")+
  geom_point(data=data_cluster,aes(x=mean_lwage_resid,y=std_lwage_resid))+
  # geom_text_repel(data=data_cluster,aes(x=mean_wage_resid,y=sigma_ou_matlab,label=firm_rank),nudge_y=0.01,show.legend=FALSE,min.segment.length = Inf)+
  theme(panel.grid = element_blank(),legend.title.align = 0.5)


plt_rich
plt_poor
ggsave(here('Results/Figures/indifference_curve_levels_rich.png'),plot = plt_rich,height = 5,width=8,scale=1.2)
ggsave(here('Results/Figures/indifference_curve_levels_poor.png'),plot = plt_poor,height = 5,width=8,scale=1.2)



# data firm in logs -------------------------------------------------------

data_cluster <- fread(here('firms_logs.csv'))
tails_cuttof <- 2.5/100
theta <- 0.16
data_cluster[,sigma_ou_matlab:=sqrt(2*theta)*sd_diff_wage_resid]
data_cluster[,lower_w:=qlnorm(tails_cuttof,mean=mean_wage_resid,sd = sd_diff_wage_resid)]
data_cluster[,upper_w:=qlnorm(1-tails_cuttof,mean=mean_wage_resid,sd = sd_diff_wage_resid)]
# data_cluster[,lower_w:=6.4]
# data_cluster[,upper_w:=100]
data_cluster[,min_w:=min(data_cluster$lower_w)]
data_cluster[,max_w:=max(data_cluster$upper_w)]
data_cluster[,firm_share:=firm_size/sum(firm_size)]
fwrite(data_cluster,file = here('firms_logs_out.csv'))




# log normal plots  -------------------------------------------------------

firm_rank_plot <- 20
mean_firm <- data_cluster[firm_rank==firm_rank_plot,mean_lwage_resid]
sigma_firm <- data_cluster[firm_rank==firm_rank_plot,sd_diff_lwage_resid]
w_min_plot <- data_cluster[firm_rank==firm_rank_plot,min_w]
w_max_plot <- data_cluster[firm_rank==firm_rank_plot,max_w]
grid_w <- seq(from=w_min_plot,to=w_max_plot,length.out=1000)
pdf <- dlnorm(grid_w,meanlog = mean_firm,sdlog = sigma_firm)
plot(grid_w,pdf)





# indifference curve plot -------------------------------------------------

data_plot_indifference <- fread(here('indifference_curve.csv'),col.names = c('Mean_log_wage','Sd_log_wage','V'))

ggplot(data_plot_indifference,aes(x=Mean_log_wage,y=Sd_log_wage,z=V))+
  geom_contour(bins=10)+
  theme_bw()
pallete <- brewer.pal(6,name = 'RdBu')


plt <- ggplot()+stat_contour(data = data_plot_indifference, aes(x = Mean_log_wage, y = Sd_log_wage, z = V, colour = ..level..), 
             breaks = round(quantile(data_plot_indifference$V, seq(0, 1, 0.07)), 5))+
  scale_color_gradient(low = pallete[6],high = pallete[1],guide = )+
  theme_bw()+
  ylab('Std. Dev. log wage')+
  labs(colour="Value Job\nOffer ")+
  xlab('Mean log wage')+
  geom_point(data=data_cluster,aes(x=mean_lwage_resid,y=sigma_ou_matlab))+
  geom_text_repel(data=data_cluster,aes(x=mean_lwage_resid,y=sigma_ou_matlab,label=firm_rank),nudge_y=0.01,show.legend=FALSE,min.segment.length = Inf)+
  theme(panel.grid = element_blank(),legend.title.align = 0.5)
plt
ggsave(here('Figures/indifference_curve.png'),plot = plt,height = 5,width=8,scale=1.5)


# wage in levels ----------------------------------------------------------

data_cluster <- fread(here('wage in levels/firms_levels.csv'))
tails_cuttof <- 2.5/100
theta <- 0.12
data_cluster[,sigma_ou_matlab:=sqrt(2*theta)*sd_diff_wage]
data_cluster[,lower_w:=qnorm(tails_cuttof,mean=mean_w,sd = sd_diff_wage)]
data_cluster[,upper_w:=qnorm(1-tails_cuttof,mean=mean_w,sd = sd_diff_wage)]
# data_cluster[,min_w:=min(data_cluster$lower_w)]
data_cluster[,min_w:=1]
data_cluster[,max_w:=max(data_cluster$upper_w)]
fwrite(data_cluster,file = here('wage in levels/firms_levels_out.csv'))

tails_cuttof <- 2.5/100
theta <- 0.16
min_sigma <- min(data_cluster$sigma_ou_matlab)
max_sigma <- max(data_cluster$sigma_ou_matlab)
min_mean <- min(data_cluster$mean_w)
max_mean <- max(data_cluster$mean_w)
n_sigma <- 20
n_mean <- 20
grid_mean <- seq(from=min_mean,to=max_mean,length.out =n_mean)
grid_sigma <- seq(from=min_sigma,to=max_sigma,length.out =n_sigma)
grid <- expand.grid(grid_mean,grid_sigma)
data_indifference <- data.table(mean_w=grid$Var1,sigma_ou_matlab=grid$Var2)
data_indifference[,sd_diff_lwage_resid:=sigma_ou_matlab/sqrt(2*theta)]
data_indifference[,firm_share:=0]
data_indifference[,rank_fun:=sd_diff_lwage_resid-mean_w]
data_indifference[,firm_rank:=rank(rank_fun)]

data_indifference[,lower_w:=qnorm(tails_cuttof,mean=mean_w,sd = sd_diff_lwage_resid)]
data_indifference[,upper_w:=qnorm(1-tails_cuttof,mean=mean_w,sd = sd_diff_lwage_resid)]
# data_indifference[,min_w:=min(data_indifference$lower_w)]
data_indifference[,min_w:=1]
data_indifference[,max_w:=max(data_indifference$upper_w)]
fwrite(data_indifference,file = here('wage in levels/data_indifference_levels.csv'))

# fake dataset for test in logs -------------------------------------------

mean_lwage_resid <- c(2,2.5,3)
sd_diff_lwage_resid <- c(0.5,0.5,0.5)
data_cluster <- data.table(mean_lwage_resid,sd_diff_lwage_resid)
tails_cuttof <- 1/100
theta <- 0.16
data_cluster[,sigma_ou_matlab:=sqrt(2*theta)*sd_diff_lwage_resid]
data_cluster[,lower_w:=qlnorm(tails_cuttof,mean=mean_lwage_resid,sd = sd_diff_lwage_resid)]
data_cluster[,upper_w:=qlnorm(1-tails_cuttof,mean=mean_lwage_resid,sd = sd_diff_lwage_resid)]
data_cluster[,min_w:=min(data_cluster$lower_w)]
data_cluster[,max_w:=max(data_cluster$upper_w)]
data_cluster[,firm_share:=0]
data_cluster[,rank_fun:=sd_diff_lwage_resid-mean_lwage_resid]
data_cluster[,firm_rank:=rank(rank_fun)]

fwrite(data_cluster,file = here('data_cluster_test.csv'))


# table with results calibration ------------------------------------------

data_cluster <- fread(here('wage in levels/firms_levels_out.csv'))
data_calibration <- fread(here('wage in levels/result_calibration.csv'))
n_digits <- 2
theta <- 0.15
# data_table_calibration <- merge(data_cluster[,.(firm_rank,mean_wage_resid,sd_diff_wage_resid,firm_share)],data_calibration)
data_calibration[,prob_offer:=format(prob_offer*100,digits=n_digits,nsmall=n_digits)]
data_calibration[,frac_employed_per_firm:=format(frac_employed_per_firm*100,digits=n_digits,nsmall=n_digits)]
data_calibration[,firm_size_share:=format(firm_size_share*100,digits=n_digits,nsmall=n_digits)]
data_calibration[,sigma_vec:=format(sigma_vec/sqrt(2*theta),digits=n_digits,nsmall=n_digits)]
data_calibration[,w_mean_vec:=format(w_mean_vec,digits=n_digits,nsmall=n_digits)]

table_calibration <- kbl(data_calibration,format = 'latex',escape =F,
                             caption = 'Calibration Results',booktabs = T,align = 'cccccc',
                             centering = T,linesep = "",row.names = F,vline='',label = 'dist_fit',
                             col.names = c('Firm Rank','$\\bar w_f$','$\\sigma_f$','Size (\\%)','$p_f$(\\%)','Size Model (\\%)'))  %>%
  kable_styling(latex_options = "HOLD_position")
table_calibration
save_kable(table_calibration,file = here('Tables','table_calibration.tex'))


# table with results unemployment sensitivity analysis --------------------
n_digits <- 2
options(scipen = 0)

data_sensitivity <- fread(here('Results/unemployment_sensitivity_analysis.csv'))
data_sensitivity[,Lower_Bound:=delta/(lambda_0+delta)]
data_sensitivity[,frac_unemployed:=format(frac_unemployed*100,digits=n_digits,nsmall=n_digits)]
data_sensitivity[,Lower_Bound:=format(Lower_Bound*100,digits=n_digits,nsmall=n_digits)]
data_sensitivity[,prob_not_accept_job_offer:=format(prob_not_accept_job_offer*100,digits=n_digits,nsmall=n_digits)]
data_sensitivity[,b:=format(b,digits=n_digits,nsmall=5)]

table_calibration <- kbl(data_sensitivity,format = 'latex',escape =F,
                         caption = 'Calibration Results',booktabs = T,align = 'cccccc',
                         centering = T,linesep = "",row.names = F,vline='',label = 'sensitivity_analysis',
                         col.names = c('$\\lambda_0$','$\\lambda_1$','$\\delta$','b','Unemployment (\\%)','Prob. Rejecting offer(\\%)','L.B. Unemployment (\\%)'))  %>%
  kable_styling(latex_options = "HOLD_position")
table_calibration
save_kable(table_calibration,file = here('Tables','unemployment_analysis.tex'))


# plot firm characteristics from model ------------------------------------

data_firms <- fread(here('wage in levels/firm_characteristics_model.csv'))
data_cluster <- fread(here('wage in levels/firms_levels_out.csv'))
data_gini_nlsy <- fread(here('wage in levels/gini_data_output.csv'))
colnames(data_gini_nlsy) <- c('firm_rank','gini_income_data')


data_plot <- merge(data_firms[,.(firm_rank,ginis_income)],data_gini_nlsy)
data_plot <- merge(data_plot,data_cluster[,.(firm_rank,sd_diff_wage_resid)])
data_plot <- melt(data_plot,id.vars = c('firm_rank','sd_diff_wage_resid'),variable.name = 'Source',value.name = 'Gini Income')
data_plot[Source=='ginis_income',Source:='Model']
data_plot[Source=='gini_income_data',Source:='Data']
plt <- ggplot(data_plot,aes(x=sd_diff_wage_resid,y=`Gini Income`,color=Source,shape=Source))+
  scale_color_manual(values=c('black','red'),)+
  labs(colour="",shape='')+
  geom_point()+
  xlab(TeX(r"($\sigma_f$)"))+
  theme_classic()
ggsave(here('Figures/gini_income_comparison.png'),plot = plt,height = 5,width=8,scale=1)

# generate dateset for same mean grid sigmas ------------------------------

data_cluster <- fread(here('wage in levels/firms_levels_out.csv'))
theta <- 0.15
tails_cuttof <- 1/100
# wage <- median(data_cluster$mean_wage_resid)
wage <- 50
sigma_min <- 1
sigma_max <- 25
n_sigmas <- 50
grid_sigmas <- seq(from=sigma_min,to=sigma_max,length.out=n_sigmas)
data_ex_sigmas <- data.table(mean_wage_resid=wage,sd_diff_wage_resid=grid_sigmas)
data_ex_sigmas[,sigma_ou_matlab:=sqrt(2*theta)*sd_diff_wage_resid]
data_ex_sigmas[,firm_share:=0]
data_ex_sigmas[,rank_fun:=sd_diff_wage_resid-mean_wage_resid]
data_ex_sigmas[,firm_rank:=rank(rank_fun)]
data_ex_sigmas[,lower_w:=qnorm(tails_cuttof,mean=mean_wage_resid,sd = sd_diff_wage_resid)]
data_ex_sigmas[,upper_w:=qnorm(1-tails_cuttof,mean=mean_wage_resid,sd = sd_diff_wage_resid)]
data_ex_sigmas[,lower_w:=6.4]
data_ex_sigmas[,upper_w:=100]
data_ex_sigmas[,min_w:=min(data_ex_sigmas$lower_w)]
data_ex_sigmas[,max_w:=max(data_ex_sigmas$upper_w)]
fwrite(data_ex_sigmas,here('wage in levels/data_sigmas_exs.csv'))

# plot sigma exercise -----------------------------------------------------

data_ex_sigmas <- fread(here('wage in levels/data_sigmas_exs.csv'))
data_plot_sigma <- fread(here('wage in levels/sigma_exercise_output_matlab.csv'))

data_plot_sigma <- merge(data_ex_sigmas[,.(firm_rank,sd_diff_wage_resid)],data_plot_sigma)
data_plot_sigma <- melt(data_plot_sigma,id.vars = c('firm_rank','sd_diff_wage_resid'),variable.name = 'Wealth',value.name = 'Value Job Offer')
names_labels <- as_labeller(
  c('Poor' = "Bootom", 'Middle' = "Middle",'Rich' = "Top"))
plt <- ggplot(data_plot_sigma,aes(x=sd_diff_wage_resid,y=`Value Job Offer`))+
  geom_line(color='#007bc8')+
  geom_point(color='#007bc8',size=0.5)+
    facet_grid(rows = 'Wealth',scales="free",labeller = names_labels)+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
  xlab(TeX(r"($\sigma_f$)"))+
  theme_bw()
   # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  # annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+theme(panel.grid = element_blank(),axis.ticks = element_line())
plt

plt1 <- ggplot(data_plot_sigma[Wealth=='Poor',],aes(x=sd_diff_wage_resid,y=`Value Job Offer`))+
  geom_line(color='#007bc8')+
  geom_point(color='#007bc8',size=0.5)+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
  theme_classic()+
  ggtitle('Poor')+
  ylab('')+
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust=0.5,size = 9))

plt2 <- ggplot(data_plot_sigma[Wealth=='Middle',],aes(x=sd_diff_wage_resid,y=`Value Job Offer`))+
  geom_line(color='#007bc8')+
  geom_point(color='#007bc8',size=0.5)+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
  ggtitle('Middle Wealth')+
  theme_classic()+
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust=0.5,size = 9))


plt3 <- ggplot(data_plot_sigma[Wealth=='Rich',],aes(x=sd_diff_wage_resid,y=`Value Job Offer`))+
  geom_line(color='#007bc8')+
  geom_point(color='#007bc8',size=0.5)+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
  xlab(TeX(r"($\sigma_f$)"))+
  ggtitle('Rich')+
  ylab('')+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust=0.5,size = 9))

plt
plt1
plt2
plt3

figure <- plot_grid(plt1,plt2,plt3,
                    ncol = 1, nrow = 3, align = "v"
)

figure
# table grids -------------------------------------------------------------

data_grid <- fread(here('wage in levels/grid_parameters.csv'))
table_grid <- kbl(data_grid,format = 'latex',escape =F,
                         caption = 'Grid Parameters',booktabs = T,align = 'lcc',
                         centering = T,linesep = "",row.names = F,vline='',label = 'grid_parameters',digits = 2) %>% 
                         add_header_above(c(" ", "Grid" = 2)) %>%
  kable_styling(latex_options = "HOLD_position")
table_grid
save_kable(table_grid,file = here('Tables','table_grid_parameters.tex'))
