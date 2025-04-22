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

closer.null  <- closer - mean(closer)

R = 10000
set.seed(42)

#Generate null bootstrap T-stats
resamples.null.closer <- replicate(R, {
  x.bar <- sample(closer.null, size = n.closer, replace = TRUE)
  mean(x.bar)/(s.closer/sqrt(n.closer))
})

##2b: Compute Bootstrap P-values

further <- finch.data$further
difference <- further - closer

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
further.null = further - mean(further)
diff.null = difference - mean(difference)
sd.further = sd(further)
sd.diff = sd(difference)

resamples.null.further <- replicate(R, {
  x.bar <- sample(further.null, size = n.further, replace = TRUE)
  mean(x.bar)/(sd.further/sqrt(n.further))
})

resamples.null.diff <- replicate(R, {
  x.bar <- sample(diff.null, size = n.diff, replace = TRUE)
  mean(x.bar)/(sd.diff/sqrt(n.diff))
})

#Bootstrapped p-vals
p_closer_boot <- mean(
  abs(resamples.null.closer) >= abs(t_closer)
)

p_further_boot <- mean(abs())
