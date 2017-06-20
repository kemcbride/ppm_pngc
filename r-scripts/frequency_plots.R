# install packages by running "install.packages('package_name')
library('ggplot2') # for plotting nicely
library('ggthemes') # for nice looking themed plots
library('plyr') # for count function, frequency data


# I got rid of scd.txt so you'd need to change this to use 
# something like complete_pairs
scd <- read.csv('/home/kelly/Dropbox/gff/scd.txt', header=FALSE)
hist_df <- count(scd, 'V1')

x_max <- 500
x_tick <- 50
freq_max <- max(hist_df['freq']) # should be 74


# I use hist_df so that we can properly apply the xlim and ylim as you like.

# Need to experiment with different shapes, geom/functions.
# plot <- ggplot(scd, aes(x=V1)) + geom_histogram(binwidth=50)
# plot <- ggplot(hist_df, aes(x=V1, y=freq)) + geom_point() + xlim(0, 750) + ylim(0, freq_max)
plot <- ggplot(hist_df, aes(x=V1, y=freq)) + geom_path() + xlim(0, x_max) + ylim(0, freq_max)
x_ticks <- seq(0, x_max, x_tick)

plot + theme_hc() +
	labs(x='Distance', y='Frequency') +
	scale_x_continuous(breaks=x_ticks, limit=c(0, x_max))
