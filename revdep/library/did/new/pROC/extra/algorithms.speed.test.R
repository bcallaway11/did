# This script tests which algorithm is faster: 2 or 3?
# The intuition is that 3 is much better with few thresholds, but behaves horribly
# as the number increases. However it seems less affected by the number of data points.
# It might change over time so let's make it easier to test again in the future

library(pROC)
library(dplyr)
#devtools::install_github("xrobin/xavamess")
library(xavamess)
library(ggplot2)
library(parallel)

# Number of observations to test
ns <- as.vector(outer(c(1), c(2:7), function(i, j) i * 10^j))
# Controls how many thresholds we get
norm.factors <- as.vector(outer(c(1, 2, 5), c(0:2), function(i, j) i * 10^j))
# Number of cores to execute on
# We want the number of physical cores, remove 1 to be sure
parallel::detectCores() / 2 - 1


# Loop over all those conditions
times.by.alg <- lapply(2:3, function(algorithm) {
	times <- lapply(rev(norm.factors), function(norm.factor) {
	#times <- autoParLapply(rev(norm.factors), function(norm.factor) {
		as.data.frame(t(sapply(ns, function(n) {
			print(sprintf("a=%s, norm=%s, n=%s", algorithm, norm.factor, n))
			# Get some data
			lab <- rbinom(n, 1, 0.5)
			d <- round(rnorm(n) * norm.factor)
			# How many thresholds do we have?
			nthr <- length(unique(d))
			if (nthr > 1000 && algorithm == 3) {
				# Algorithm 3 is out here anyway, no need to waste time to test it
				return(c(n=n, norm.factor=norm.factor, nthr=nthr, algorithm=algorithm, rep(NA, 3)))
			}
			# Repeat 5 times and take the median time
			time <- apply(replicate(5,  system.time(roc(lab, d, algorithm=algorithm, levels = c(0, 1), direction = "<"))), 1, median)
			return(c(n=n, norm.factor=norm.factor, nthr=nthr, algorithm=algorithm, time[1:3]))
		})))
	}) # Physical cores, not logical !!!
	#}, .maxCores = 3) # Physical cores, not logical !!!
	times.df = bind_rows(times)
})
times.by.alg.df <- bind_rows(times.by.alg)
times.by.alg.df$algorithm <- as.factor(times.by.alg.df$algorithm)


# Plot the data
library(ggplot2)
ggplot(times.by.alg.df) + geom_point(aes(n, user.self, color=algorithm)) + scale_x_log10() + scale_y_log10()
ggplot(times.by.alg.df) + geom_point(aes(n, user.self, color=nthr)) + facet_grid(algorithm ~ .) + scale_x_log10() + scale_y_log10()
ggplot(times.by.alg.df) + geom_point(aes(n, user.self, color=log(nthr))) + facet_grid(algorithm ~ .) + scale_x_log10() + scale_y_log10()
ggplot(times.by.alg.df) + geom_point(aes(nthr, user.self, color=n)) + facet_grid(algorithm ~ .) + scale_x_log10() + scale_y_log10()
ggplot(times.by.alg.df) + geom_point(aes(nthr, n, color=user.self)) + facet_grid(algorithm ~ .)

# Algorithm 3 is linear with nthr * n?
ggplot(times.by.alg.df) + geom_point(aes(nthr * n, user.self)) + facet_grid(algorithm ~ .)
plot(nthr * n ~ user.self, na.omit(times.by.alg.df %>% filter(algorithm==3)))


# Test algorithm 2
times.by.alg.df2 <- times.by.alg.df %>% filter(algorithm == 2, n > 200)
lm.2 <- lm(user.self ~ n * nthr, times.by.alg.df2)
# nthr gives low, barely significant but negative estimates which don't make sense, so remove it...
lm.2 <- lm(user.self ~ n + 0, times.by.alg.df2)
summary(lm.2)
plot(lm.2)

times.by.alg.df2$predicted.user.self <- predict(lm.2, times.by.alg.df2)
plot(times.by.alg.df2$user.self, times.by.alg.df2$predicted.user.self)

plot(times.by.alg.df2$n, times.by.alg.df2$user.self)
grid <- expand.grid(n=ns, nthr=ns)
grid$prod = grid$n * grid$nthr
grid <- grid[order(grid$n),]
lines(grid$n, predict(lm.2, grid))

# Test algorithm 3
times.by.alg.df3 <- times.by.alg.df %>% 
	filter(algorithm == 3, n > 200) %>%
	mutate(prod = nthr * n)
lm.3 <- lm(user.self ~ n:nthr + 0, times.by.alg.df3)
summary(lm.3)
plot(lm.3)

times.by.alg.df3$predicted.user.self <- predict(lm.3, times.by.alg.df3)
plot(times.by.alg.df3$user.self, times.by.alg.df3$predicted.user.self)

plot(times.by.alg.df3$n * times.by.alg.df3$nthr, times.by.alg.df3$user.self)
grid <- expand.grid(n=ns, nthr=ns)
grid$prod = grid$n * grid$nthr
grid <- grid[order(grid$prod),]
lines(grid$prod, predict(lm.3, grid))


# Predict time on initial data
times.by.alg.df$user.self.predicted.2 <- predict(lm.2, times.by.alg.df)
times.by.alg.df$user.self.predicted.3 <- predict(lm.3, times.by.alg.df)
times.by.alg.df$predicted.best <- ifelse(times.by.alg.df$user.self.predicted.3 < times.by.alg.df$user.self.predicted.2, 3, 2)
ggplot(times.by.alg.df) + geom_point(aes(user.self.predicted.2, user.self.predicted.3, color=n))+ scale_x_log10() + scale_y_log10()
ggplot(times.by.alg.df) + geom_point(aes(user.self.predicted.2, user.self.predicted.3, color=predicted.best))+ scale_x_log10() + scale_y_log10()
ggplot(times.by.alg.df) + geom_point(aes(n, nthr, color=predicted.best))+ scale_x_log10() + scale_y_log10()


### Final formula:
# Algorithm 2: user.self = 2.959e-07 * n
# Algorithm 3: user.self = 5.378e-09 * n * nthr
# Reduction:
# 2.959e-07 * n = 5.378e-09 * n * nthr
# 2.959e-07 / 5.378e-09 = nthr
# 55 = nthr