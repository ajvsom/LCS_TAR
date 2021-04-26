# TAR & LCS Representations -----------------------------------------------
#AR(1) Difference Equation
T  = 10
a = seq(-1.5, 1.5, .1) 
y_ts  = matrix(NA, length(a), T)  
y_ts[,1] <- 4

for (i in 1:length(a)){
  for (j in 2:T){
    y_ts[i,j] = a[i]*y_ts[i,(j-1)] + 4
  }
}
y_ts = as.data.frame(y_ts)
names(y_ts) = c(paste0('t', seq(1,10)))
y_ts$a = seq(-1.5, 1.5, .1) 
y_ts_long = reshape(y_ts, 
                    direction = "long",
                    varying = list(names(y_ts[1:10])),
                    v.names = "Y",
                    idvar = c("cond"),
                    timevar = "Time")

y_ts_long$MODEL = 'TAR'

#UNIVARIATE CHANGE SCORE MODEL
T  = 10
p = seq(-2.5,.5,.1)
x_ts  = matrix(NA, length(p), T)  
x_ts[,1] <- 4

for (i in 1:length(p)){
  for (j in 2:T){
    x_ts[i,j] = (1+p[i])*x_ts[i,(j-1)] + 4
  }
}
x_ts = as.data.frame(x_ts)
names(x_ts) = c(paste0('t', seq(1,10)))
x_ts$p = seq(-2.5, .5, .1) 
x_ts_long = reshape(x_ts, 
                    direction = "long",
                    varying = list(names(y_ts[1:10])),
                    v.names = "Y",
                    idvar = c("cond"),
                    timevar = "Time")
x_ts_long$MODEL = 'LCS'

ts_long_merged = rbind(y_ts_long, x_ts_long)


# Figure 1: TAR Trajectories ----------------------------------------------

library(ggplot2)
library(gridExtra)
theme_update(plot.title = element_text(hjust = 0.5))

y_polarize = subset(y_ts_long, a >= -1.5 & a <= -1.1)
y_polarize$a = as.factor(y_polarize$a)

y_ntrend = subset(y_ts_long, a == -1)
y_ntrend$a = as.factor(y_ntrend$a)

y_nconv = subset(y_ts_long, a > -1 & a < 0)
y_nconv$a = as.factor(y_nconv$a)
levels(y_nconv$a)[levels(y_nconv$a)=="-0.0999999999999999"] <- "-0.1"

y_none = subset(y_ts_long, a == 0)
y_none$a = as.factor(y_none$a)

y_pconv = subset(y_ts_long, a > 0 & a < 1)
y_pconv$a = as.factor(y_pconv$a)

y_trend = subset(y_ts_long, a == 1)
y_trend$a = as.factor(y_trend$a)

y_exp = subset(y_ts_long, a >= 1.1 & a <= 1.5)
y_exp$a = as.factor(y_exp$a)

c1 = ggplot(data = y_polarize) + 
  aes(x = Time, y = Y, group = a) + 
  geom_line(aes(linetype=a)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Volatility",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) +theme(legend.key.size = unit(0.5, "cm"))

c2 = ggplot(data = y_ntrend) + 
  aes(x = Time, y = Y, group = a) + 
  geom_line(aes(linetype=a)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Periodic Change",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) +theme(legend.key.size = unit(0.5, "cm"))

c3 = ggplot(data = y_nconv) + 
  aes(x = Time, y = Y, group = a) + 
  geom_line(aes(linetype=a)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Oscillatory Convergence",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) +theme(legend.key.size = unit(0.5, "cm"))

c4 = ggplot(data = y_none) + 
  aes(x = Time, y = Y, group = a) + 
  geom_line(aes(linetype=a)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Stasis",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10))+theme(legend.key.size = unit(0.5, "cm"))

c5 = ggplot(data = y_pconv) + 
  aes(x = Time, y = Y, group = a) + 
  geom_line(aes(linetype=a)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Smooth Convergence",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) +theme(legend.key.size = unit(0.5, "cm"))

c6 = ggplot(data = y_trend) + 
  aes(x = Time, y = Y, group = a) + 
  geom_line(aes(linetype=a)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Constant Growth",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) +theme(legend.key.size = unit(0.5, "cm"))

c7 = ggplot(data = y_exp) + 
  aes(x = Time, y = Y, group = a) + 
  geom_line(aes(linetype=a)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Explosive Growth", x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) +theme(legend.key.size = unit(0.5, "cm"))

grid.arrange(c1, c2, c3, c4, c5, c6, c7, nrow = 3)


# Figure 2: LCS Trajectories ----------------------------------------------

x_polarize = subset(x_ts_long, p >= -2.5 & p <= -2.1)
x_polarize$p = as.factor(x_polarize$p)

x_ntrend = subset(x_ts_long, p == -2)
x_ntrend$p = as.factor(x_ntrend$p)

x_nconv = subset(x_ts_long, p >= -1.9 & p < -1)
x_nconv$p = as.factor(x_nconv$p)

x_none = subset(x_ts_long, p == -1)
x_none$p = as.factor(x_none$p)

x_pconv = subset(x_ts_long, p >= -.9 & p < 0)
x_pconv$p = as.factor(x_pconv$p)
levels(x_pconv$p)[levels(x_pconv$p)=="-0.0999999999999996"] <- "-0.1"

x_trend = subset(x_ts_long, p == 0)
x_trend$p = as.factor(x_trend$p)

x_exp = subset(x_ts_long, p >= .1 & p <= .5)
x_exp$p = as.factor(x_exp$p)

c8 = ggplot(data = x_polarize) + 
  aes(x = Time, y = Y, group = p) + 
  geom_line(aes(linetype=p)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Volatility",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10))+theme(legend.key.size = unit(0.5, "cm"))

c9 = ggplot(data = x_ntrend) + 
  aes(x = Time, y = Y, group = p) + 
  geom_line(aes(linetype=p)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Periodic Change",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10))+theme(legend.key.size = unit(0.5, "cm"))

c10 = ggplot(data = x_nconv) + 
  aes(x = Time, y = Y, group = p) + 
  geom_line(aes(linetype=p)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Oscillatory Convergence",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10))+theme(legend.key.size = unit(0.5, "cm"))

c11 = ggplot(data = x_none) + 
  aes(x = Time, y = Y, group = p) + 
  geom_line(aes(linetype=p)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Stasis",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10))+theme(legend.key.size = unit(0.5, "cm"))

c12 = ggplot(data = x_pconv) + 
  aes(x = Time, y = Y, group = p) + 
  geom_line(aes(linetype=p)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Smooth Convergence",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10))+theme(legend.key.size = unit(0.5, "cm"))

c13 = ggplot(data = x_trend) + 
  aes(x = Time, y = Y, group = p) + 
  geom_line(aes(linetype=p)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Constant Growth",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10))+theme(legend.key.size = unit(0.5, "cm"))

c14 = ggplot(data = x_exp) + 
  aes(x = Time, y = Y, group = p) + 
  geom_line(aes(linetype=p)) + 
  #geom_point(aes(shape=a)) + 
  labs(title="Explosive Growth", x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) +theme(legend.key.size = unit(0.5, "cm"))

grid.arrange(c8, c9, c10, c11, c12, c13, c14, nrow = 3)



#Import data
library(readxl)
table_for_graphs <- read_excel("OneDrive - Michigan State University/Research/Current_Projects/LCSM_Sims/figures_tables/table_for_graphs.xlsx")
View(table_for_graphs)

df = table_for_graphs

#Load packages
library(ggplot2)
library(gridExtra)
theme_update(plot.title = element_text(hjust = 0.5))

# Figure 3: Bias  ----------------------------------------------------------
# T = 5 ==================================================================
df$N = as.factor(df$N)
bias_t5 = ggplot(data = subset(df, T == 5)) + 
  aes(x = CHANGE, y = as.numeric(BIAS)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="T = 5",x="Beta", y = "Bias") + 
  scale_x_continuous(breaks=seq(-1.5, 1.5, .1)) + 
  scale_color_manual(values=c("red", "blue"))
## T = 10 ==================================================================
bias_t10 = ggplot(data = subset(df, T == 10)) + 
  aes(x = CHANGE, y = as.numeric(BIAS)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="T = 10",x="Beta", y = "Bias") + 
  scale_x_continuous(breaks=seq(-1.5, 1.5, .1)) + 
  scale_color_manual(values=c("red", "blue"))
## T = 30 ==================================================================
bias_t30 = ggplot(data = subset(df, T == 30)) + 
  aes(x = CHANGE, y = as.numeric(BIAS)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="T = 30",x="Beta", y = "Bias") + 
  scale_x_continuous(breaks=seq(-1.5, 1.5, .1)) + 
  scale_color_manual(values=c("red", "blue"))

grid.arrange(bias_t5, bias_t10, bias_t30)

# Predicted vs Actual  -----------------------------------------------------
# Prediction Functions -----------------------------------------------------
project_tar = function(a, T){
  y_ts = matrix(NA, length(a), T)
  y_ts[,1] <- 4
  for (i in 1:length(a)){
    for (j in 2:T){
      y_ts[i,j] = a[i]*y_ts[i,(j-1)] + 4
    }
  }
  y_ts = as.data.frame(y_ts)
  names(y_ts) = c(paste0('t', seq(1,T)))
  y_ts_long = reshape(y_ts, 
                      direction = "long",
                      varying = list(names(y_ts[1:T])),
                      v.names = "Y",
                      idvar = c("cond"), #note: cond = coefficient
                      timevar = "Time")
  a = list(y_ts, y_ts_long)
  return(a)
}
project_lcsm = function(p, T){
  x_ts = matrix(NA, length(p), T)
  x_ts[,1] <- 4
  for (i in 1:length(p)){
    for (j in 2:T){
      x_ts[i,j] = (1+p[i])*x_ts[i,(j-1)] + 4
    }
  }
  x_ts = as.data.frame(x_ts)
  names(x_ts) = c(paste0('t', seq(1,T)))
  x_ts_long = reshape(x_ts, 
                      direction = "long",
                      varying = list(names(x_ts[1:T])),
                      v.names = "Y",
                      idvar = c("cond"), #note: cond = coefficient
                      timevar = "Time")
  a = list(x_ts, x_ts_long)
  return(a)
}


# T = 5 -------------------------------------------------------------------
# Generate Data -----------------------------------------------------------
# Actual
b = seq(-1.5, 1.5, .1)
temp_t5 = as.data.frame(project_tar(b, 5)[2])
COEF = seq(-1.5, 1.5, .1)
N = NA
actualt5 = cbind(temp_t5, COEF, N)
actualt5$N = 200
actualt5$MODEL = 'b'

# TAR
temp_tarT5 = project_tar(subset(df, T == 5 & MODEL == "TAR")$EST, 5)
temp_tarT5 = as.data.frame(temp_tarT5[2])
COEF = rep(seq(-1.5, 1.5, .1), 5*2) #T = 5 * levels(N) = 2
N = rep(c(rep(200,31), rep(500,31)), 5) #levels(b) = 31; T = 5
tarT5 = cbind(temp_tarT5, COEF, N)
tarT5$MODEL = 'TAR'

# LCS
temp_lcsT5 = project_lcsm(subset(df, T == 5 & MODEL == "LCS")$EST, 5)
temp_lcsT5 = as.data.frame(temp_lcsT5[2])
COEF = rep(seq(-1.5, 1.5, .1), 5*2) #T = 5 * levels(N) = 2
N = rep(c(rep(200,31), rep(500,31)), 5) #levels(b) = 31; T = 5
lcsT5 = cbind(temp_lcsT5, COEF, N)
lcsT5$MODEL = 'LCS'

predT5 = rbind(actualt5, tarT5, lcsT5)
predT5 = subset(predT5, select = -c(cond))
predT5$N = as.factor(predT5$N)
predT5$COEF = as.factor(predT5$COEF)
levels(predT5$COEF)[levels(predT5$COEF)=="-0.0999999999999999"] <- "-0.1"

# Graphing Functions --------------------------------------
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
# Graph Nonstationary Cases for T = 5 --------------------------------------
# Negative Nonstationary Cases
n15t5 = ggplot(data = subset(predT5, COEF == -1.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n14t5 = ggplot(data = subset(predT5, COEF == -1.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n13t5 = ggplot(data = subset(predT5, COEF == -1.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n12t5 = ggplot(data = subset(predT5, COEF == -1.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n11t5 = ggplot(data = subset(predT5, COEF == -1.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n10t5 = ggplot(data = subset(predT5, COEF == -1.0)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.0",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

# Positive Nonstationary Cases
p10t5 = ggplot(data = subset(predT5, COEF == 1.0)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.0",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

p11t5 = ggplot(data = subset(predT5, COEF == 1.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

p12t5 = ggplot(data = subset(predT5, COEF == 1.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

p13t5 = ggplot(data = subset(predT5, COEF == 1.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

p14t5 = ggplot(data = subset(predT5, COEF == 1.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

p15t5 = ggplot(data = subset(predT5, COEF == 1.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

grid.arrange(n15t5+theme(legend.position='hidden'), 
             n14t5+theme(legend.position='hidden'), 
             n13t5+theme(legend.position='hidden'), 
             n12t5+theme(legend.position='hidden'), 
             n11t5+theme(legend.position='hidden'), 
             n10t5+theme(legend.position='hidden'), 
             p10t5+theme(legend.position='hidden'), 
             p11t5+theme(legend.position='hidden'), 
             p12t5+theme(legend.position='hidden'), 
             p13t5+theme(legend.position='hidden'), 
             p14t5+theme(legend.position='hidden'), 
             p15t5+theme(legend.position='hidden'), 
             nrow = 3)

# Graph Stationary Cases for T = 5 --------------------------------------
# Negative Stationary Cases
n09t5 = ggplot(data = subset(predT5, COEF == -0.9)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.9",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n08t5 = ggplot(data = subset(predT5, COEF == -0.8)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.8",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n07t5 = ggplot(data = subset(predT5, COEF == -0.7)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.7",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n06t5 = ggplot(data = subset(predT5, COEF == -0.6)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.6",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n05t5 = ggplot(data = subset(predT5, COEF == -0.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n04t5 = ggplot(data = subset(predT5, COEF == -0.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n03t5 = ggplot(data = subset(predT5, COEF == -0.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n02t5 = ggplot(data = subset(predT5, COEF == -0.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

n01t5 = ggplot(data = subset(predT5, COEF == -0.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 

# Positive Stationary Cases
p09t5 = ggplot(data = subset(predT5, COEF == 0.9)) + 
 aes(x = as.numeric(Time), y = as.numeric(Y)) + 
 geom_line(aes(color = MODEL, linetype=N)) + 
 labs(title="b = 0.9",x="Time", y = "Y values") + 
 scale_x_continuous(breaks=seq(0, 5)) 

p08t5 = ggplot(data = subset(predT5, COEF == 0.8)) + 
 aes(x = as.numeric(Time), y = as.numeric(Y)) + 
 geom_line(aes(color = MODEL, linetype=N)) + 
 labs(title="b = 0.8",x="Time", y = "Y values") + 
 scale_x_continuous(breaks=seq(0, 5)) 

p07t5 = ggplot(data = subset(predT5, COEF == 0.7)) + 
 aes(x = as.numeric(Time), y = as.numeric(Y)) + 
 geom_line(aes(color = MODEL, linetype=N)) + 
 labs(title="b = 0.7",x="Time", y = "Y values") + 
 scale_x_continuous(breaks=seq(0, 5)) 

p06t5 = ggplot(data = subset(predT5, COEF == 0.6)) + 
 aes(x = as.numeric(Time), y = as.numeric(Y)) + 
 geom_line(aes(color = MODEL, linetype=N)) + 
 labs(title="b = 0.6",x="Time", y = "Y values") + 
 scale_x_continuous(breaks=seq(0, 5)) 

p05t5 = ggplot(data = subset(predT5, COEF == 0.5)) + 
 aes(x = as.numeric(Time), y = as.numeric(Y)) + 
 geom_line(aes(color = MODEL, linetype=N)) + 
 labs(title="b = 0.5",x="Time", y = "Y values") + 
 scale_x_continuous(breaks=seq(0, 5)) 

p04t5 = ggplot(data = subset(predT5, COEF == 0.4)) + 
 aes(x = as.numeric(Time), y = as.numeric(Y)) + 
 geom_line(aes(color = MODEL, linetype=N)) + 
 labs(title="b = 0.4",x="Time", y = "Y values") + 
 scale_x_continuous(breaks=seq(0, 5)) 

p03t5 = ggplot(data = subset(predT5, COEF == 0.3)) + 
 aes(x = as.numeric(Time), y = as.numeric(Y)) + 
 geom_line(aes(color = MODEL, linetype=N)) + 
 labs(title="b = 0.3",x="Time", y = "Y values") + 
 scale_x_continuous(breaks=seq(0, 5)) 

p02t5 = ggplot(data = subset(predT5, COEF == 0.2)) + 
 aes(x = as.numeric(Time), y = as.numeric(Y)) + 
 geom_line(aes(color = MODEL, linetype=N)) + 
 labs(title="b = 0.2",x="Time", y = "Y values") + 
 scale_x_continuous(breaks=seq(0, 5)) 

p01t5 = ggplot(data = subset(predT5, COEF == 0.1)) + 
 aes(x = as.numeric(Time), y = as.numeric(Y)) + 
 geom_line(aes(color = MODEL, linetype=N)) + 
 labs(title="b = 0.1",x="Time", y = "Y values") + 
 scale_x_continuous(breaks=seq(0, 5)) 

p0t5 = ggplot(data = subset(predT5, COEF == 0)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 5)) 


grid.arrange(n09t5+theme(legend.position='hidden'), 
             n08t5+theme(legend.position='hidden'), 
             n07t5+theme(legend.position='hidden'), 
             n06t5+theme(legend.position='hidden'), 
             n05t5+theme(legend.position='hidden'), 
             n04t5+theme(legend.position='hidden'), 
             n03t5+theme(legend.position='hidden'), 
             n02t5+theme(legend.position='hidden'), 
             n01t5+theme(legend.position='hidden'), 
             p0t5+theme(legend.position='hidden'), 
             p01t5+theme(legend.position='hidden'), 
             p02t5+theme(legend.position='hidden'),
             p03t5+theme(legend.position='hidden'),
             p04t5+theme(legend.position='hidden'),
             p05t5+theme(legend.position='hidden'),
             p06t5+theme(legend.position='hidden'),
             p07t5+theme(legend.position='hidden'),
             p08t5+theme(legend.position='hidden'),
             p09t5+theme(legend.position='hidden'),
             nrow = 4)

grid.arrange(n15t5+theme(legend.position='hidden'), 
             n14t5+theme(legend.position='hidden'), 
             n13t5+theme(legend.position='hidden'), 
             n12t5+theme(legend.position='hidden'), 
             n11t5+theme(legend.position='hidden'), 
             n10t5+theme(legend.position='hidden'), 
             n09t5+theme(legend.position='hidden'), 
             n08t5+theme(legend.position='hidden'), 
             n07t5+theme(legend.position='hidden'), 
             n06t5+theme(legend.position='hidden'), 
             n05t5+theme(legend.position='hidden'), 
             n04t5+theme(legend.position='hidden'), 
             n03t5+theme(legend.position='hidden'), 
             n02t5+theme(legend.position='hidden'), 
             n01t5+theme(legend.position='hidden'), 
             p0t5+theme(legend.position='hidden'), 
             p01t5+theme(legend.position='hidden'), 
             p02t5+theme(legend.position='hidden'),
             p03t5+theme(legend.position='hidden'),
             p04t5+theme(legend.position='hidden'),
             p05t5+theme(legend.position='hidden'),
             p06t5+theme(legend.position='hidden'),
             p07t5+theme(legend.position='hidden'),
             p08t5+theme(legend.position='hidden'),
             p09t5+theme(legend.position='hidden'),
             p10t5+theme(legend.position='hidden'), 
             p11t5+theme(legend.position='hidden'), 
             p12t5+theme(legend.position='hidden'), 
             p13t5+theme(legend.position='hidden'), 
             p14t5+theme(legend.position='hidden'), 
             p15t5+theme(legend.position='hidden'), 
             nrow = 4)
# T = 10 -------------------------------------------------------------------
# Generate Data -----------------------------------------------------------
# Actual
b = seq(-1.5, 1.5, .1)
temp_t10 = as.data.frame(project_tar(b, 10)[2])
COEF = seq(-1.5, 1.5, .1)
N = NA
actualt10 = cbind(temp_t10, COEF, N)
actualt10$N = 200
actualt10$MODEL = 'b'

# TAR
temp_tarT10 = project_tar(subset(df, T == 10 & MODEL == "TAR")$EST, 10)
temp_tarT10 = as.data.frame(temp_tarT10[2])
COEF = rep(seq(-1.5, 1.5, .1), 10*2) #T = 10 * levels(N) = 2
N = rep(c(rep(200,31), rep(500,31)), 10) #levels(b) = 31; T = 10
tarT10 = cbind(temp_tarT10, COEF, N)
tarT10$MODEL = 'TAR'

# LCS
temp_lcsT10 = project_lcsm(subset(df, T == 10 & MODEL == "LCS")$EST, 10)
temp_lcsT10 = as.data.frame(temp_lcsT10[2])
COEF = rep(seq(-1.5, 1.5, .1), 10*2) #T = 10 * levels(N) = 2
N = rep(c(rep(200,31), rep(500,31)), 10) #levels(b) = 31; T = 10
lcsT10 = cbind(temp_lcsT10, COEF, N)
lcsT10$MODEL = 'LCS'

predT10 = rbind(actualt10, tarT10, lcsT10)
predT10 = subset(predT10, select = -c(cond))
predT10$N = as.factor(predT10$N)
predT10$COEF = as.factor(predT10$COEF)
levels(predT10$COEF)[levels(predT10$COEF)=="-0.0999999999999999"] <- "-0.1"

# Graph Nonstationary Cases for T = 10 --------------------------------------
# Negative Nonstationary Cases
n15T10 = ggplot(data = subset(predT10, COEF == -1.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n14T10 = ggplot(data = subset(predT10, COEF == -1.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n13T10 = ggplot(data = subset(predT10, COEF == -1.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n12T10 = ggplot(data = subset(predT10, COEF == -1.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n11T10 = ggplot(data = subset(predT10, COEF == -1.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n10T10 = ggplot(data = subset(predT10, COEF == -1.0)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.0",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

# Positive Nonstationary Cases
p10T10 = ggplot(data = subset(predT10, COEF == 1.0)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.0",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p11T10 = ggplot(data = subset(predT10, COEF == 1.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p12T10 = ggplot(data = subset(predT10, COEF == 1.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p13T10 = ggplot(data = subset(predT10, COEF == 1.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p14T10 = ggplot(data = subset(predT10, COEF == 1.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p15T10 = ggplot(data = subset(predT10, COEF == 1.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

grid.arrange(n15T10+theme(legend.position='hidden'), 
             n14T10+theme(legend.position='hidden'), 
             n13T10+theme(legend.position='hidden'), 
             n12T10+theme(legend.position='hidden'), 
             n11T10+theme(legend.position='hidden'), 
             n10T10+theme(legend.position='hidden'), 
             p10T10+theme(legend.position='hidden'), 
             p11T10+theme(legend.position='hidden'), 
             p12T10+theme(legend.position='hidden'), 
             p13T10+theme(legend.position='hidden'), 
             p14T10+theme(legend.position='hidden'), 
             p15T10+theme(legend.position='hidden'), 
             nrow = 3)

# Graph Stationary Cases for T = 10 --------------------------------------
# Negative Stationary Cases
n09T10 = ggplot(data = subset(predT10, COEF == -0.9)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.9",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n08T10 = ggplot(data = subset(predT10, COEF == -0.8)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.8",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n07T10 = ggplot(data = subset(predT10, COEF == -0.7)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.7",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n06T10 = ggplot(data = subset(predT10, COEF == -0.6)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.6",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n05T10 = ggplot(data = subset(predT10, COEF == -0.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n04T10 = ggplot(data = subset(predT10, COEF == -0.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n03T10 = ggplot(data = subset(predT10, COEF == -0.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n02T10 = ggplot(data = subset(predT10, COEF == -0.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

n01T10 = ggplot(data = subset(predT10, COEF == -0.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

# Positive Stationary Cases
p09T10 = ggplot(data = subset(predT10, COEF == 0.9)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.9",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p08T10 = ggplot(data = subset(predT10, COEF == 0.8)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.8",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p07T10 = ggplot(data = subset(predT10, COEF == 0.7)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.7",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p06T10 = ggplot(data = subset(predT10, COEF == 0.6)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.6",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p05T10 = ggplot(data = subset(predT10, COEF == 0.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p04T10 = ggplot(data = subset(predT10, COEF == 0.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p03T10 = ggplot(data = subset(predT10, COEF == 0.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p02T10 = ggplot(data = subset(predT10, COEF == 0.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p01T10 = ggplot(data = subset(predT10, COEF == 0.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 

p0T10 = ggplot(data = subset(predT10, COEF == 0)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 10)) 


grid.arrange(n09T10+theme(legend.position='hidden'), 
             n08T10+theme(legend.position='hidden'), 
             n07T10+theme(legend.position='hidden'), 
             n06T10+theme(legend.position='hidden'), 
             n05T10+theme(legend.position='hidden'), 
             n04T10+theme(legend.position='hidden'), 
             n03T10+theme(legend.position='hidden'), 
             n02T10+theme(legend.position='hidden'), 
             n01T10+theme(legend.position='hidden'), 
             p0T10+theme(legend.position='hidden'), 
             p01T10+theme(legend.position='hidden'), 
             p02T10+theme(legend.position='hidden'),
             p03T10+theme(legend.position='hidden'),
             p04T10+theme(legend.position='hidden'),
             p05T10+theme(legend.position='hidden'),
             p06T10+theme(legend.position='hidden'),
             p07T10+theme(legend.position='hidden'),
             p08T10+theme(legend.position='hidden'),
             p09T10+theme(legend.position='hidden'),
             nrow = 4)

grid.arrange(n15T10+theme(legend.position='hidden'), 
             n14T10+theme(legend.position='hidden'), 
             n13T10+theme(legend.position='hidden'), 
             n12T10+theme(legend.position='hidden'), 
             n11T10+theme(legend.position='hidden'), 
             n10T10+theme(legend.position='hidden'), 
             n09T10+theme(legend.position='hidden'), 
             n08T10+theme(legend.position='hidden'), 
             n07T10+theme(legend.position='hidden'), 
             n06T10+theme(legend.position='hidden'), 
             n05T10+theme(legend.position='hidden'), 
             n04T10+theme(legend.position='hidden'), 
             n03T10+theme(legend.position='hidden'), 
             n02T10+theme(legend.position='hidden'), 
             n01T10+theme(legend.position='hidden'), 
             p0T10+theme(legend.position='hidden'), 
             p01T10+theme(legend.position='hidden'), 
             p02T10+theme(legend.position='hidden'),
             p03T10+theme(legend.position='hidden'),
             p04T10+theme(legend.position='hidden'),
             p05T10+theme(legend.position='hidden'),
             p06T10+theme(legend.position='hidden'),
             p07T10+theme(legend.position='hidden'),
             p08T10+theme(legend.position='hidden'),
             p09T10+theme(legend.position='hidden'),
             p10T10+theme(legend.position='hidden'), 
             p11T10+theme(legend.position='hidden'), 
             p12T10+theme(legend.position='hidden'), 
             p13T10+theme(legend.position='hidden'), 
             p14T10+theme(legend.position='hidden'), 
             p15T10+theme(legend.position='hidden'), 
             nrow = 4)
# T = 30 -------------------------------------------------------------------
# Generate Data -----------------------------------------------------------
# Actual
b = seq(-1.5, 1.5, .1)
temp_T30 = as.data.frame(project_tar(b, 30)[2])
COEF = seq(-1.5, 1.5, .1)
N = NA
actualT30 = cbind(temp_T30, COEF, N)
actualT30$N = 200
actualT30$MODEL = 'b'

# TAR
temp_tarT30 = project_tar(subset(df, T == 30 & MODEL == "TAR")$EST, 30)
temp_tarT30 = as.data.frame(temp_tarT30[2])
COEF = rep(seq(-1.5, 1.5, .1), 30*2) #T = 30 * levels(N) = 2
N = rep(c(rep(200,31), rep(500,31)), 30) #levels(b) = 31; T = 30
tarT30 = cbind(temp_tarT30, COEF, N)
tarT30$MODEL = 'TAR'

# LCS
temp_lcsT30 = project_lcsm(subset(df, T == 30 & MODEL == "LCS")$EST, 30)
temp_lcsT30 = as.data.frame(temp_lcsT30[2])
COEF = rep(seq(-1.5, 1.5, .1), 30*2) #T = 30 * levels(N) = 2
N = rep(c(rep(200,31), rep(500,31)), 30) #levels(b) = 31; T = 30
lcsT30 = cbind(temp_lcsT30, COEF, N)
lcsT30$MODEL = 'LCS'

predT30 = rbind(actualT30, tarT30, lcsT30)
predT30 = subset(predT30, select = -c(cond))
predT30$N = as.factor(predT30$N)
predT30$COEF = as.factor(predT30$COEF)
levels(predT30$COEF)[levels(predT30$COEF)=="-0.0999999999999999"] <- "-0.1"

# Graph Nonstationary Cases for T = 30 --------------------------------------
# Negative Nonstationary Cases
n15T30 = ggplot(data = subset(predT30, COEF == -1.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n14T30 = ggplot(data = subset(predT30, COEF == -1.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n13T30 = ggplot(data = subset(predT30, COEF == -1.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n12T30 = ggplot(data = subset(predT30, COEF == -1.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n11T30 = ggplot(data = subset(predT30, COEF == -1.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n10T30 = ggplot(data = subset(predT30, COEF == -1.0)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -1.0",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

# Positive Nonstationary Cases
p10T30 = ggplot(data = subset(predT30, COEF == 1.0)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.0",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p11T30 = ggplot(data = subset(predT30, COEF == 1.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p12T30 = ggplot(data = subset(predT30, COEF == 1.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p13T30 = ggplot(data = subset(predT30, COEF == 1.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p14T30 = ggplot(data = subset(predT30, COEF == 1.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p15T30 = ggplot(data = subset(predT30, COEF == 1.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 1.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

grid.arrange(n15T30+theme(legend.position='hidden'), 
             n14T30+theme(legend.position='hidden'), 
             n13T30+theme(legend.position='hidden'), 
             n12T30+theme(legend.position='hidden'), 
             n11T30+theme(legend.position='hidden'), 
             n10T30+theme(legend.position='hidden'), 
             p10T30+theme(legend.position='hidden'), 
             p11T30+theme(legend.position='hidden'), 
             p12T30+theme(legend.position='hidden'), 
             p13T30+theme(legend.position='hidden'), 
             p14T30+theme(legend.position='hidden'), 
             p15T30+theme(legend.position='hidden'), 
             nrow = 3)

# Graph Stationary Cases for T = 30 --------------------------------------
# Negative Stationary Cases
n09T30 = ggplot(data = subset(predT30, COEF == -0.9)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.9",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n08T30 = ggplot(data = subset(predT30, COEF == -0.8)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.8",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n07T30 = ggplot(data = subset(predT30, COEF == -0.7)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.7",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n06T30 = ggplot(data = subset(predT30, COEF == -0.6)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.6",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n05T30 = ggplot(data = subset(predT30, COEF == -0.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n04T30 = ggplot(data = subset(predT30, COEF == -0.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n03T30 = ggplot(data = subset(predT30, COEF == -0.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n02T30 = ggplot(data = subset(predT30, COEF == -0.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

n01T30 = ggplot(data = subset(predT30, COEF == -0.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = -0.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

# Positive Stationary Cases
p09T30 = ggplot(data = subset(predT30, COEF == 0.9)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.9",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p08T30 = ggplot(data = subset(predT30, COEF == 0.8)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.8",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p07T30 = ggplot(data = subset(predT30, COEF == 0.7)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.7",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p06T30 = ggplot(data = subset(predT30, COEF == 0.6)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.6",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p05T30 = ggplot(data = subset(predT30, COEF == 0.5)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.5",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p04T30 = ggplot(data = subset(predT30, COEF == 0.4)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.4",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p03T30 = ggplot(data = subset(predT30, COEF == 0.3)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.3",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p02T30 = ggplot(data = subset(predT30, COEF == 0.2)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.2",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p01T30 = ggplot(data = subset(predT30, COEF == 0.1)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0.1",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 

p0T30 = ggplot(data = subset(predT30, COEF == 0)) + 
  aes(x = as.numeric(Time), y = as.numeric(Y)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="b = 0",x="Time", y = "Y values") + 
  scale_x_continuous(breaks=seq(0, 30, 5)) 


grid.arrange(n09T30+theme(legend.position='hidden'), 
             n08T30+theme(legend.position='hidden'), 
             n07T30+theme(legend.position='hidden'), 
             n06T30+theme(legend.position='hidden'), 
             n05T30+theme(legend.position='hidden'), 
             n04T30+theme(legend.position='hidden'), 
             n03T30+theme(legend.position='hidden'), 
             n02T30+theme(legend.position='hidden'), 
             n01T30+theme(legend.position='hidden'), 
             p0T30+theme(legend.position='hidden'), 
             p01T30+theme(legend.position='hidden'), 
             p02T30+theme(legend.position='hidden'),
             p03T30+theme(legend.position='hidden'),
             p04T30+theme(legend.position='hidden'),
             p05T30+theme(legend.position='hidden'),
             p06T30+theme(legend.position='hidden'),
             p07T30+theme(legend.position='hidden'),
             p08T30+theme(legend.position='hidden'),
             p09T30+theme(legend.position='hidden'),
             nrow = 4)


grid.arrange(n15T30+theme(legend.position='hidden'), 
             n14T30+theme(legend.position='hidden'), 
             n13T30+theme(legend.position='hidden'), 
             n12T30+theme(legend.position='hidden'), 
             n11T30+theme(legend.position='hidden'), 
             n10T30+theme(legend.position='hidden'), 
             n09T30+theme(legend.position='hidden'), 
             n08T30+theme(legend.position='hidden'), 
             n07T30+theme(legend.position='hidden'), 
             n06T30+theme(legend.position='hidden'), 
             n05T30+theme(legend.position='hidden'), 
             n04T30+theme(legend.position='hidden'), 
             n03T30+theme(legend.position='hidden'), 
             n02T30+theme(legend.position='hidden'), 
             n01T30+theme(legend.position='hidden'), 
             p0T30+theme(legend.position='hidden'), 
             p01T30+theme(legend.position='hidden'), 
             p02T30+theme(legend.position='hidden'),
             p03T30+theme(legend.position='hidden'),
             p04T30+theme(legend.position='hidden'),
             p05T30+theme(legend.position='hidden'),
             p06T30+theme(legend.position='hidden'),
             p07T30+theme(legend.position='hidden'),
             p08T30+theme(legend.position='hidden'),
             p09T30+theme(legend.position='hidden'),
             p10T30+theme(legend.position='hidden'), 
             p11T30+theme(legend.position='hidden'), 
             p12T30+theme(legend.position='hidden'), 
             p13T30+theme(legend.position='hidden'), 
             p14T30+theme(legend.position='hidden'), 
             p15T30+theme(legend.position='hidden'), 
             nrow = 4)

# Figure 7: T1 Error  ----------------------------------------------------------

# T = 5 ==================================================================
df$N = as.factor(df$N)
type1_t5 = ggplot(data = subset(df, T == 5)) + 
  aes(x = CHANGE, y = as.numeric(T1ERROR)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="T = 5",x="Beta", y = "Type 1 Error") + 
  scale_x_continuous(breaks=seq(-1.5, 1.5, .1)) + 
  scale_color_manual(values=c("red", "blue")) +
  geom_hline(yintercept = .05)
## T = 10 ==================================================================
type1_t10 = ggplot(data = subset(df, T == 10)) + 
  aes(x = CHANGE, y = as.numeric(T1ERROR)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="T = 10",x="Beta", y = "Type 1 Error") + 
  scale_x_continuous(breaks=seq(-1.5, 1.5, .1)) + 
  scale_color_manual(values=c("red", "blue")) +
  geom_hline(yintercept = .05)
## T = 30 ==================================================================
type1_t30 = ggplot(data = subset(df, T == 30)) + 
  aes(x = CHANGE, y = as.numeric(T1ERROR)) + 
  geom_line(aes(color = MODEL, linetype=N)) + 
  labs(title="T = 30",x="Beta", y = "Type 1 Error") + 
  scale_x_continuous(breaks=seq(-1.5, 1.5, .1)) + 
  scale_color_manual(values=c("red", "blue")) +
  geom_hline(yintercept = .05)

grid.arrange(type1_t5, type1_t10, type1_t30)





