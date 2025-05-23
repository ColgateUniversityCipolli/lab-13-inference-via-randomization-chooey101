library(tidyverse)
library(patchwork)
library(e1071)
###1a: Compute the error
finch.data <- read.csv("zebrafinches.csv")
n = nrow(finch.data)
further.sum <- finch.data |>
  summarize(
    mean = mean(further),
    sd = sd(further),
    skew = skewness(further)
  )

close.t.stat <- t.test(x=finch.data$further, mu = 0, alternative = "less")

t.stat <- close.t.stat$statistic

fZ.at.t <- dnorm(t.stat, mean = 0, sd = 1)

error.val <- ((further.sum$skew)/(sqrt(n)))*(((2*t.stat^2)+1)/6)*fZ.at.t

###1b: compute error from -10 to 10

t.seq <- seq(-10, 10, length.out=1000)

fZ <- dnorm(t.seq, mean = 0, sd = 1)

skew.further <- further.sum$skew
# Edgeworth error
error.seq <- (skew.further / sqrt(n)) *
  ((2 * t.seq^2 + 1) / 6) *
  fZ

error.df <- tibble(
  t     = t.seq,
  error = error.seq
)


error.plot <- ggplot(error.df, aes(x = t, y = error)) +
  geom_line(size = 1) +
  labs(
    title = "Edgeworth Approximation Error vs. t",
    x     = "t statistic",
    y     = "Approximation Error"
  ) +
  theme_minimal()

##1c: compute tail probability
alpha = 0.05
assumed.t <- abs(qnorm(alpha))
assumed.fz <- dnorm(assumed.t)
skew.t <-
  tail.prob <- (((skew.further)/6*(0.10*alpha))*((2*assumed.t^2)+1)*assumed.fz)^2

##Question 2: Perform bootstrapping procedure

###2a: 
closer <- finch.data$closer


n.closer <- length(closer)        
s.closer <- sd(closer)
mu.obs      <- mean(closer)


R = 1000
set.seed(42)

closer.resamples <- tibble(xbar = rep(NA, R), sd = rep(NA, R), t.stat = rep(NA, R))

#Generate null bootstrap T-stats
 for (i in 1:R) {
  x.star   <- sample(closer, size = n.closer, replace = TRUE)
  closer.resamples$xbar[i] <- mean(x.star)
  closer.resamples$sd[i] <- sd(x.star)
  closer.resamples$t.stat[i] <- (closer.resamples$xbar[i]-mu.obs) / (closer.resamples$sd[i] / sqrt(n.closer))
 }

#Plot approximate null distribution
closer_dist <- ggplot(closer.resamples, aes(x=t.stat)) +
  geom_histogram(bins = 30,
                 fill  = "skyblue",
                 color = "black") +
  labs(title = "Approximate distribution of Bootstrapped Closer Data",
       x     = "t",
       y     = "Density") +
  theme_minimal()

##2b: Compute Bootstrap P-values

further <- finch.data$further
difference <- further - closer


further.resamples <- tibble(xbar = rep(NA, R), sd = rep(NA, R), t.stat = rep(NA, R))
diff.resamples <- tibble(xbar = rep(NA, R), sd = rep(NA, R), t.stat = rep(NA, R))

# observed T‑statistics
t_closer    <- t.test(closer, mu = 0, alternative = "less")$statistic
p_closer_tt <- t.test(closer, mu = 0, alternative = "less")$p.value

t_further    <- t.test(further, mu = 0, alternative = "less")$statistic
p_further_tt <- t.test(further, mu = 0, alternative = "less")$p.value

t_diff      <- t.test(difference, mu = 0, alternative = "less")$statistic
p_diff_tt   <- t.test(difference, mu = 0, alternative = "less")$p.value

#Bootstrap T-stats for further and difference vals 
n.further = length(further)
n.diff = length(difference)
mu.further = mean(further)
mu.diff = mean(difference)

for (i in 1:R) {
  x.star   <- sample(further, size = n.further, replace = TRUE)
  further.resamples$xbar[i] <- mean(x.star)
  further.resamples$sd[i] <- sd(x.star)
  further.resamples$t.stat[i] <- (further.resamples$xbar[i]-mu.further) / (further.resamples$sd[i] / sqrt(n.further))
}

for (i in 1:R) {
  x.star   <- sample(difference, size = n.diff, replace = TRUE)
  diff.resamples$xbar[i] <- mean(x.star)
  diff.resamples$sd[i] <- sd(x.star)
  diff.resamples$t.stat[i] <- (diff.resamples$xbar[i]-mu.diff) / (diff.resamples$sd[i] / sqrt(n.diff))
}

##Remember to fix bootstrap statistics

#Bootstrapped p-vals
p_closer_boot  <- mean(abs(closer.resamples$t.stat)  >= abs(t_closer))
p_further_boot <- mean(abs(further.resamples$t.stat) >= abs(t_further))
p_diff_boot    <- mean(abs(diff.resamples$t.stat)    >= abs(t_diff))

#Part 2c: Compute 5th percentile under shifted null dist

t05_closer_boot  <- quantile(closer.resamples$t.stat,  probs = 0.05)
t05_further_boot <- quantile(further.resamples$t.stat, probs = 0.05)
t05_diff_boot    <- quantile(diff.resamples$t.stat,    probs = 0.05)

# Compute theoretical t critical values at α = 0.05 (left tail)
t05_closer_theo  <- qt(0.05, df = n.closer - 1)
t05_further_theo <- qt(0.05, df = n.further - 1)
t05_diff_theo    <- qt(0.05, df = n.diff    - 1)

#Part 2d: Compute confidence intervals
probs <- c(0.025, 0.975)

# 1. Closer data
q_closer        <- quantile(closer.resamples$t.stat,  probs = probs)
ci_closer_boot  <- c(
  mean(closer) - q_closer[2] * (s.closer / sqrt(n.closer)),
  mean(closer) - q_closer[1] * (s.closer / sqrt(n.closer))
)
tci_closer      <- t.test(closer,  conf.level = 0.95)$conf.int

# 2. Further data
q_further        <- quantile(further.resamples$t.stat, probs = probs)
ci_further_boot  <- c(
  mean(further) - q_further[2] * (sd.further / sqrt(n.further)),
  mean(further) - q_further[1] * (sd.further / sqrt(n.further))
)
tci_further      <- t.test(further, conf.level = 0.95)$conf.int

# 3. Difference data
q_diff           <- quantile(diff.resamples$t.stat,    probs = probs)
ci_diff_boot     <- c(
  mean(difference) - q_diff[2] * (sd.diff / sqrt(n.diff)),
  mean(difference) - q_diff[1] * (sd.diff / sqrt(n.diff))
)
tci_diff         <- t.test(difference, conf.level = 0.95)$conf.int

#Part 3: Randomization Testing
#3a

rand.t.closer <- replicate(R, {
  signs <- sample(c(-1,1), n.closer, replace = TRUE)
  mean(signs * closer.null)  / (s.closer   / sqrt(n.closer))
})

rand.t.further <- replicate(R, {
  signs <- sample(c(-1,1), n.further, replace = TRUE)
  mean(signs * further.null) / (sd.further / sqrt(n.further))
})

rand.t.diff <- replicate(R, {
  signs <- sample(c(-1,1), n.diff, replace = TRUE)
  mean(signs * diff.null)    / (sd.diff    / sqrt(n.diff))
})

#3b: Compute Randomization P values
p_closer_rand  <- mean(abs(rand.t.closer)  >= abs(t_closer))
p_further_rand <- mean(abs(rand.t.further) >= abs(t_further))
p_diff_rand    <- mean(abs(rand.t.diff)    >= abs(t_diff))

#3c: Compute confidence intervals

#This function will be used within the next function to calculate the p value for each iteration of the loop
rand_p_val <- function(x, mu0) {
  n    <- length(x)
  s    <- sd(x)
  x0   <- x - mu0
  t_obs <- abs(mean(x) - mu0) / (s / sqrt(n))
  rand_t <- replicate(R, {
    signs <- sample(c(-1,1), n, replace = TRUE)
    mean(signs * x0) / (s / sqrt(n))
  })
  mean(abs(rand_t) >= t_obs)
}
#Loops until it finds the first p value that is less than .05 for both upper and lower bounds
rand_ci <- function(x) {
  m   <- mean(x)
  lower <- m
  while(rand_p_val(x, lower) < 0.05) {
    lower <- lower - step
  }
  upper <- m
  while(rand_p_val(x, upper) < 0.05) {
    upper <- upper + step
  }
  c(lower = lower, upper = upper)
}

ci_closer_rand  <- rand_ci(closer)
ci_further_rand <- rand_ci(further)
ci_diff_rand    <- rand_ci(difference)


