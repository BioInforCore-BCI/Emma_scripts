
#Take coverage histogram data produced from AGeNT pipeline
#Adapters trimmed, aligned to hg38 reference genome, duplicates purged w/ LocatIt


library(ggplot2)
library(RColorBrewer)

#Make list of all hist coverage files
print(files <- list.files(pattern="_hist.txt"))

#Loop for each file in list, produce boxplots of coverage distribution for each target probe
for (i in 1:length(files)) {

df <- read.delim(files[i], header = FALSE) 

df <- df[- grep("all",df$V1),]

filename <- files[i]

ggplot(df,
	aes(x=V5,
		y=V4)) +
	geom_boxplot(outlier.shape=".",
				outlier.colour="blue") +
	coord_flip() +
	theme(axis.text.x = element_text(angle = 90, size=5, hjust=1)) +
	ylab("Target Probe") +
	xlab("Coverage") +
	theme(axis.title=element_text(size=14,face="bold")) +
	scale_x_continuous(expand = c(0, 0)) +
	geom_vline(xintercept = mean(df$V5), colour = "Red", size = 0.5)
	
	ggsave(sub(".txt",".pdf",files[i]), width = 40, height = 7)
	
}
	
