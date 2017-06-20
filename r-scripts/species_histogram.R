# install packages by running "install.packages('package_name')
library('ggplot2') # for plotting nicely
library('ggthemes') # for nice looking themed plots
library('plyr') # for count function, frequency data


# Speicfy the path to the data file
gcf_list <- read.csv('/home/kelly/Dropbox/gff/complete_list', header=FALSE)

species_map <- read.csv('/home/kelly/Dropbox/gff/taxonomy_genus_map', header=FALSE)

gcf_list <- merge(gcf_list, species_map, by='V1')


# Experiment with using V2, V3, V4, V5 (V6 is too specific)
# If you change it on one line, you MUST change it on both lines
hist_df <- count(gcf_list, 'V2')

# Produce stats BEFORE filtering
sy <- summary(hist_df$freq)
sy


# NOTE: If you want to filter out low frequency results, use:
hist_df <- hist_df[hist_df$freq > ceiling(mean(hist_df$freq)),]
hist_df <- hist_df[hist_df$freq > ceiling(mean(hist_df$freq)),]

freq_max <- max(hist_df['freq']) # Relies on hist_df, so must be after its creation


plot <- ggplot(hist_df, aes(x=reorder(V2, freq), y=freq)) + geom_bar(stat='identity', width=0.3, fill='darkblue') + ylim(0, freq_max)

pdf(file='Rplots.pdf', width=5, height=10)

plot + theme_hc() +
	theme(axis.text.x = element_text(size=8)) +
	labs(x='Genus', y='Frequency') + 
	coord_flip()
