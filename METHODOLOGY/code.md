# Extensions of threshold analysis in network meta-analysis

This study provides a framework to derive how changes in treatment effects within invariant intervals impact on heterogeneity of pairwise and network meta-analysis, local and global methods of inconsistency of network meta-analysis estimates.

## 1 Network meta-analysis with netmeta package

### 1.1 Libraries:

```{r}
library(netmeta)
library(ggplot2)
library(dplyr)
```

### 1.2 Data:

For the shake of this analysis dummy data was generated to indicate a generic example of the methods that will be used.

```{r}
# Column - > study_id
study_id <- c(1,1, #Study 1
              2,2, #Study 2
              3,3, #Study 3
              4,4, #Study 4
              5,5, #Study 5
              6,6, #Study 6
              7,7) #Study 7
# Group of patiens in each study arm
n <- c(90,85, # n.1, n.2 Study 1
       95,95, # n.1, n.2 Study 2
       99,97, # ...
       101,111,
       105,102,
       109,100,
       194,201)
# Patients treated in each study arm
r<- c(13,18, # r.1, r.2 Study 1
      14,20, # r.1, r.2 Study 2
      17,22, #...
      18,24,
      20,24,
      22,30,
      42,54)
# Treatment_id
t <-c('treatment-0','treatment-1', #Study 1
      'treatment-0','treatment-2', #Study 2
      'treatment-0','treatment-1', #Study 3
      'treatment-0','treatment-2', #Study 4
      'treatment-0','treatment-2', #Study 5
      'treatment-0','treatment-2', #Study 6
      'treatment-1','treatment-2') 
data <- cbind(study_id,r,n,t)
data
```

### 1.3 Network meta-analysis with netmeta package:

**Pairwise** function transform meta-analysis data from two arm-based formats in contrast based format that is needed as input to netmata function.

**Pairwise** is used to perform pairwise comparisons between two specific **treatments** or **interventions** in a network meta-analysis. This method allows the researcher to estimate and compare the effectiveness of multiple treatments even when **direct** head-to-head comparisons are limited or unavailable.

In this example for **summary measure** we use OR (Odds Ratios because they can interpret easier for the reader). Other options include RR(Risk Ratio), MD(Mean Difference), SMD(Standard Mean Difference), RD(Risk Difference), HR(Hazard Ratio).

```{r}
my_data_l=pairwise(treat=t,event=r,n,studlab=study_id,sm="OR")
my_data_l
```

**Netconnection** provide information on network connectivity: (these information are also included in netmeta)

```{r}
netconnection(treat1 = treat1, treat2= treat2, studlab=study_id, data = my_data_l)
```

**Netmeta** is used for conduction network meta analysis to compare and rank multiple treatments or interventions when there are both direct and indirect evidence from a network of studies.

To run this function we need information about the studies, treatments and their effect sizes. So we will run the network meta analysis using **netmata** and we will take the outputs of pairwise as an input for this function. **TE** column, contains the effect size of all comparisons and **seTE** the respective standard error of each comparison.

We use as a **reference group** the treatment-0, **95%** the level to calculate **confidence intervals** and we state comb.random for random effects model.

```{r}
results.NMA=netmeta(TE, seTE, treat1, treat2, studlab, data=my_data_l, sm="OR",level.ma=0.95, ref="treatment-0", comb.random=TRUE)
summary(results.NMA)
```

**Data:**

-   This network meta-analysis include **k =** **6 studies** in total, where **5** of them are **two-arm studies** and there is **one three-arm study**.

-   The number of pairwise comparisons: m = 6

-   The number of treatments: n = 3

**Results:\
**The next table for the common effect model shows the fitted values for each comparison in our network meta-analysis. The ***Q*** column in this table indicates which comparison contributes substantially to the overall inconsistency in our network. For our example all studies have approximately the same value of Q.

The treatment effect estimate for **random effects model** show that there is **not a statistically significant difference** between the reference treatment compared to the treatment-1 while 1 is included on the confidence interval of the effect estimates and therefore the p-value is above 0.05 indicating for a not statistically significant difference. However, there is a **statistical significant difference** between the reference treatment compared with the treatment-2, while 1 is not included on the confidence interval and the p-value is \< 0.05 (p-value = 0.01).

**Heterogeneity / Inconsistency:\
**The heterogeneity in network meta analysis is composed from two components: The **within-study** heterogeneity, and inconsistency **between** study designs.

-   **Within study** heterogeneity refers to the variations in the results within a single study, such as differences in effect sizes or outcomes among participants or subgroups. It reflects the diversity of data within a particular research investigation.

-   **Between** study heterogeneity pertains to differences in the results across multiple studies, it captures the variability in effect sizes, study designs or populations across various research studies.

In network meta-analysis, there are also several metrics that are commonly used to assess and quantify heterogeneity.

-   **Q-statistic** is a measure of heterogeneity that tests the null hypothesis that there is no variation in treatment effects among studies. It compares the observed differences between study estimates and the expected differences under the assumption of heterogeneity. If the Q-statistic is significant, it suggest that there is substantial heterogeneity, indicating that treatment effects are not consistent across studies.

-   **I squared** is a measure that quantifies the proportion of the total variation in the total treatment effects across studies that is due to heterogeneity rather than chance. It is expressed as a percentage and the interpretation is given with the rule of thumb as low, moderate and high with the values of 25%, 50%, 75% respectively.

-   **Tau-squared** is a metric that estimates the between-study variance or the amount of true heterogeneity present in the NMA. It provides a numerical measure of the extent of heterogeneity and is used to calculate the confidence interval for the summary treatment effect. A larger tau-squared suggests greater variation in treatment effects among studies.

**Inconsistency Factor\
Netsplit** to examine the contribution of direct and indirect evidence and to test inconsistency in network meta - analysis.

```{r}
net = netsplit(results.NMA, digits = 2)
forest(net)
```

**Network Graph**\
**Net-graph** shows the direct connection between the studies, for this example:

```{r}
netgraph(results.NMA, number.of.studies = TRUE)
```

-   2 studies compare **treatment-0** with **treatment-1**

-   4 study compare **treatment-0** with **treatment-2**

-   1 study compare **treatment-1** with **treatment-2**

**Effect Estimate Table**\
**Net-league** function is a square matrix showing all pairwise comparisons in network meta-analysis. We evaluate the **inconsistency** between our studies for direct and mixed(both direct and indirect comparisons) comparisons. On the upper right triangle we can see the direct comparisons and on the lower left triangle we can see the mixed comparisons for either fixed or random effects model.

```{r}
netleague(results.NMA)
```

**Treatment Ranking\
**Then we use **net-rank** function to rank the treatments and indicate which treatment has the highest effect. It allows us to generate a ranking of treatments that suggest which treatment is more or less likely to produce the larger benefits. For our example we consider that small values are not in favor of the treatment.

```{r}
netrank(results.NMA,small.values = 'bad')
```

**Probability of being best:\
**This is the probability that a particular innervation is the best among all interventions, the values range from 0 to 1, where 1 indicates high confidence that the corresponding innervation is the best.

**Surface Under the Cumulative Ranking Curve:**\
SUCRA is a summary statistic that combines the probabilities of each innervation being at a various possible ranks. **Rankogram** is a graphical representation of the cumulative ranking probabilities for each innervation. The x-axis represents the possible ranks, and the y-axis represents the cumulative probability of being at or below each rank.

```{r}
ran <- rankogram(results.NMA, small.values = "undesirable")
plot(ran,)
print(ran, cumulative.rankprob = TRUE)
```

**Results:**\
The results rank treatment-2 the treatment with the highest effect size(netrank) and the highest cumulative probability (SUCRA) and therefore the "best" treatment followed by treatment-1. For this assumption we need to also investigate further for the various sources of heterogeneity in the existing studies.

**Forrest Plot Visualization\
**We can visual represent the results of the analysis with a **forest plot**, which provide a summary of the estimated treatment effects for each intervention in the network of studies. We can also interpret the results for the statistical significant difference, however this plot do not provide any information for the heterogeneity between the included studies.

```{r}
par(mfrow = c(1,2))
forest(results.NMA, ref="treatment-0",fs.axis=10,pooled = "random")
forest(results.NMA, ref="treatment-0",fs.axis=10,pooled = "common")
```

## 2 Threshold analysis in network meta analysis with nmathresh

### 2.1 Libraries:

```{r}
library(nmathresh)
```

### 2.2 Data:

The example of the included studies is the same as the previous example for the network meta-analysis but the format of the data is different we make the changes accordingly.

```{r}
###Column - > study_id
study_id <- c(1,2,3,4,5,6,7)
#plcaebo(treatment-0) group in studies
n.1 <- c(90,95,99,101,105,109,194)
#treated in placebo(treatment-0) group
r.1 <- c(13,14,17,18,20,22,42)
#intervation group in studies
n.2 <- c(85,95,97,111,102,100,201)
#treated in intervation group
r.2 <- c(18,20,22,24,24,30,54)
#treatment_id
t.1 <-c(1,1,1,1,1,1,2)
t.2 <- c(2,3,2,3,3,3,3)
#n-arms
n_arms <- c(2,2,2,2,2,2,2)

my_study <- data.frame(study_id,n.1,r.1,n.2,r.2,t.1,t.2,n_arms)
my_study
```

### 2.3 Threshold analysis with nmathresh (Study Level):

First we need to do some calculations in order to run the nma_thresh:

*Log OR for **two-arm trials**, arm 2 vs. 1:*

```{r}
my_study$lor.1 <- with(my_study, log(r.2 * (n.1 - r.1) / ((n.2 - r.2) * r.1)) )
my_study$k.1 <- my_study$t.2   # Arm 2 treatment
my_study$b.1 <- my_study$t.1   # Reference treatment
```

*LOR variances and co-variances, likelihood **covariance matrix V**:\
(Since we do not use 3 arm studies the else if argument is not used for this example)*

```{r}
n <- nrow(my_study) # n number of studies
V.diag <- as.list(rep(NA, n))
attach(my_study)
for (i in 1:n){
  if (my_study$n_arms[i] == 2){
    V.diag[[i]] <- 1/r.1[i] + 1/r.2[i] + 1/(n.1[i]-r.1[i]) + 1/(n.2[i]-r.2[i])
  }
  else if (my_study$n_arms[i] == 3){
    v1 <- 1/r.1[i] + 1/r.2[i] + 1/(n.1[i]-r.1[i]) + 1/(n.2[i]-r.2[i])
    v2 <- 1/r.1[i] + 1/r.3[i] + 1/(n.1[i]-r.1[i]) + 1/(n.3[i]-r.3[i])
    # Covariance term
    c1 <- 1/r.1[i] + 1/(n.1[i] - r.1[i])
    V.diag[[i]] <- matrix(c(v1, c1, c1, v2), nrow = 2)
  }
}
detach(my_study)

library(Matrix)
V <- bdiag(V.diag)
V
```

*Reshape the data to have one row per contrast per study:*

```{r}
my_study_l <- reshape(my_study, varying = c("lor.1", "b.1", "k.1"), 
                    timevar = "c", idvar = "study_id", direction = "long")
my_study_l
```

*Construct the **design matrix X**, with a row for each contrast and K-1 columns (parameters):*

```{r}
N <- nrow(my_study_l)   # number of data points
K <- length(unique(c(my_study_l$b, my_study_l$k)))   # Number of treatments

X <- matrix(0, nrow = N, ncol = K-1)

for (i in 1:N){
  X[i, my_study_l$k[i]-1] <- 1
  if (my_study_l$b[i] != 1){
    X[i, my_study_l$b[i]-1] <- -1
  }
}
X
```

Now we can run the nma_thresh, but we also need:

1.  the posterior **summaries of the variable\
    **We got the **estimate of treatment effect, difference between first and second treatment** from the netmeta package and that we are going to use as an input for mean.dk.

2.  the posterior **covariance matrix**\
    We got **Variance - covariance matrix for fixed-effects models** from the netmata package and that we are going to use as an input for post. We could also use the came matrix for the random effects model but we need to also change the nmatype accordingly.

3.  We select opt.max = False considering that the optimal treatment is the one which minimizes the log odds.

```{r}
thresh <- nma_thresh(mean.dk = results.NMA$seTE[1:(K-1)], 
                     lhood = V, 
                     post = results.NMA$Cov.fixed[1:2,1:2],
                     nmatype = "fixed",
                     X = X,
                     opt.max = FALSE)
thresh
```

### 2.4 Forrest Plot with thresh_forest

Display using a forest plot, along with 95% confidence intervals for LORs:\
We need:\
***Labels*** *:*

```{r}
my_study_l$lab <- rep(NA, nrow(my_study_l))
for (i in 1:nrow(my_study_l)) {
  my_study_l$lab[i] <- paste0(my_study_l$study_id[i], " (", my_study_l$k[i], " vs. ", my_study_l$b[i], ")")
}
```

***95% Confidence Intervals**:*

```{r}
my_study_l$CI2.5 <- my_study_l$lor + qnorm(.025)*sqrt(diag(V))
my_study_l$CI97.5 <- my_study_l$lor + qnorm(.975)*sqrt(diag(V))
```

Sort the studies by proportion of CI covered by invariant interval. **Coverage** \<1 means that the CI extends beyond the bias invariant threshold, and the bias is plausible within the data.

```{r}
my_study_l$coverage <- apply(cbind(thresh$thresholds$lo / (my_study_l$CI2.5 - my_study_l$lor), 
                                 thresh$thresholds$hi / (my_study_l$CI97.5 - my_study_l$lor)), 1, min, na.rm = TRUE)
```

Now we cant plot the forest plot with the invariant intervals using thresh_forest:

```{r}
thresh_forest(thresh, 
              y = lor, CI.lo = CI2.5, CI.hi = CI97.5, 
              label = lab, orderby = coverage, data = my_study_l,
              CI.title = "95% Confidence Interval", y.title = "Log OR", 
              label.title = "Study (Contrast)", xlab = "Log Odds Ratio", 
              II.title = expression("Invariant Interval and "*tilde(k)*"*"),
              II.cols = rgb(.72, .80, .93),
              xlim = c(-3, 2), refline = 0, digits = 2)
```

## Methods

### Extensions of Threshold Analysis Function

Step 1 Set an upper bound for invariant interval, if the upper boundary is infinity then we make it 3 times the logOR:

```{r}
# Invariant interval
II.lo = round(my_study_l$lor + thresh$thresholds$lo,2) 
II.hi = round(my_study_l$lor + thresh$thresholds$hi,2)

# Function to set an upper and lower bound to II
for (i in 1:length(thresh$thresholds$hi)){
  if (is.infinite(thresh$thresholds$hi[i])) {
    if (my_study_l$lor[i] < 0) {
      thresh$thresholds$hi[i] = 3*abs(my_study_l$lor[i]) + my_study_l$lor[i]
    } else {
      thresh$thresholds$hi[i] = 3*my_study_l$lor[i]
    }
  }
  if(is.infinite(thresh$thresholds$lo[i])){
    if (my_study_l$lor[i] < 0) {
      thresh$thresholds$lo[i] = 3*my_study_l$lor[i]
    } else {
      thresh$thresholds$lo[i] = my_study_l$lor[i] - 3*my_study_l$lor[i]
    }
  }
}
```

```{r}
thresh_forest(thresh, 
                y = lor, CI.lo = CI2.5, CI.hi = CI97.5, 
                label = lab, orderby = coverage, data = my_study_l,
                CI.title = "95% Confidence Interval", y.title = "Log OR", 
                label.title = "Study (Contrast)", xlab = "Log Odds Ratio", 
                II.title = expression("Invariant Interval and "*tilde(k)*"*"),
                II.cols = rgb(.72, .80, .93),
                xlim = c(-3, 2), refline = 0, digits = 2)
```

Modifying the TE and run the netmeta until we reach the boundaries (Per Study)

```{r}
#Modifying the TE and run the netmeta until we reach the boundaries (Per Study)

# Lists to store the metrics in the function
heterogeneity = list()
heterogeneity_tau = list()
decomposition = list()
inconsistency_factor = list()
inconsistency_pval = list()
# Fixed Variables
treatment_effects = my_data_l$TE
# Invariant interval
II.lo = round(my_study_l$lor + thresh$thresholds$lo,2) 
II.hi = round(my_study_l$lor + thresh$thresholds$hi,2)

for (i in 1:length(my_data_l$TE)) {
  my_data_l$TE = treatment_effects
  print('------ Next Loop ------')
  print(paste0('------ Study index ------: ', i))
  print(my_data_l$TE)
  for (j in 1:floor((II.hi[i] - II.lo[i])/0.1)){
    my_data_l$TE[i] = (II.lo[i] + 0.1*j)
    results.NMA=netmeta(TE, seTE, treat1, treat2, studlab, data=my_data_l, sm="OR",level.ma=0.95, ref="treatment-0", comb.random=TRUE)
    
    heterogeneity = append(heterogeneity, c(results.NMA$tau2))
    heterogeneity_tau = append(heterogeneity, c(results.NMA$tau))#tau 
    
    decomp.NMA = decomp.design(results.NMA)
    decomposition = append(decomposition,c(decomp.NMA$Q.decomp$pval[1]))
    
    net = netsplit(results.NMA, digits = 2)
    inconsistency_factor = append(inconsistency_factor, c(net$compare.random['statistic'][1,1]))
    inconsistency_pval = append(inconsistency_pval, c(net$compare.random['p'][1,1]))
    
    #print(paste0('------ TE Adjusted Times ------: ', j))
    #print(my_data_l$TE[i])
  }
}
```

Print the Results:

```{r}
column.names <- c("|Heterogeneity|","|Decomposition|","|Inconsistency Factor|","|Inconsistency p-value|")

for (i in 1:length(my_data_l$TE)){
  if(i>1){
    x = 0
    for(k in 1:(i-1)){
      x = (x + floor((II.hi[k] - II.lo[k])/0.1))
    }
    print(paste0('____________________________ Study index ', i, ' ____________________________'))
    print('______________________________________________________________________')
    metric_1 = heterogeneity[(x+1):(x + floor((II.hi[i] - II.lo[i])/0.1))]
    #metric_5 = heterogeneity_tau[(x+1):(x + floor((II.hi[i] - II.lo[i])/0.1))]
    metric_2 = decomposition[(x+1):(x + floor((II.hi[i] - II.lo[i])/0.1))]
    metric_3 = inconsistency_factor[(x+1):(x + floor((II.hi[i] - II.lo[i])/0.1))]
    metric_4 = inconsistency_pval[(x+1):(x + floor((II.hi[i] - II.lo[i])/0.1))]
    metrics = array(c(metric_1,
                      #metric_5,
                      metric_2,
                      metric_3,
                      metric_4),
                    dim = c(floor((II.hi[i] - II.lo[i])/0.1),4))
    colnames(metrics) <- column.names
    print(metrics)
  }
  else{
    print(paste0('____________________________ Study index ', i, ' ____________________________'))
    print('______________________________________________________________________')
    metrics = array(c(heterogeneity,decomposition,inconsistency_factor,inconsistency_pval),dim = c(floor((II.hi[i] - II.lo[i])/0.1),4))
    colnames(metrics) <- column.names
    print(metrics)
  }
}
```

**Visualization of the results:**

Functions

```{r}
color.gradient <- function(x, colors=c("green","orange","red"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
metrics_plot <- function(title,metric,i){
  if(i>1){
    x = 0
    for(k in 1:(i-1)){
      x = x + floor((II.hi[k] - II.lo[k])/0.01)
    }
    j = 1:floor((II.hi[i] - II.lo[i])/0.01)
    metric1 = metric[(x+1):(x + floor((II.hi[i] - II.lo[i])/0.01))]
    
    return(plot(x =(II.lo[i] + 0.01*j),
                y = unlist(metric)[(x+1):(x + floor((II.hi[i] - II.lo[i])/0.01))],
                col= color.gradient(unlist(metric1)),
                pch= 19,
                xlim = range(II.lo[i],II.hi[i]),
                xlab = 'Treatment Effect',
                ylab= title
                #main = paste0('Study ',my_data_l$studlab[i])
    )
    )
  }
  else{
    j = 1:floor((II.hi[i] - II.lo[i])/0.01)
    metric1 = metric[1:floor((II.hi[i] - II.lo[i])/0.01)]
    return(plot(x =(II.lo[i] + 0.01*j),
                y = unlist(metric)[1:floor((II.hi[i] - II.lo[i])/0.01)],
                col= color.gradient(unlist(metric1)),
                pch= 19,
                xlim = range(II.lo[i],II.hi[i]),
                xlab = 'Treatment Effect',
                ylab= title
                #main = paste0('Study ',my_data_l$studlab[i])
    )
    )
  }
  
}


all_metrics <- function(i){
  par(mfrow = c(2,2))
  metrics_plot('Heterogeineity',heterogeneity,i)
  metrics_plot('Decomposition',decomposition,i)
  metrics_plot('Inconsistency Factor',inconsistency_factor,i)
  metrics_plot('Inconsistency p-value',inconsistency_pval,i)
  mtext(paste0('Study ',my_data_l$studlab[i]), side = 3, line = -1.5, outer = TRUE)
}
```

Visuals

```{r}
i=1
for(i in 1:length(my_study_l$study_id))
  all_metrics(i)
```
