# macaque_pathogen_boxplots.R
# last edited 11/29/2017
# Note: these counts are pulled from DESeq normalized organism counts for macaque files.
# yes, it's hard-coded.  No, I'm not thrilled about it.

library(ggplot2)

jejun_c <- c(2380.573307,2277.62073,2417.293344,2137.77424,3749.917374,1996.533656,2444.551564,2182.633859,4019.267158,4571.085076,4010.87967,7031.019227)
jejun_e <- c(18904.16813,5865.781823,9784.665069,59733.75544,18271.36416,6077.190856,18627.7718,16718.17236,11248.01198,18819.92801,4956.260934,3534.010263)
coli_c <- c(1236.756498,1283.847534,1394.007192,1671.277015,2566.061224,1192.630548,1567.000547,1105.574039,2763.792663,2624.520348,2729.50058,3958.857503)
coli_e <- c(10902.55425,3493.353307,5895.143933,29134.08702,10552.62764,2913.862564,9643.625872,9074.753754,6071.685353,9628.984983,3304.572914,3105.328545)
pylor_c <- c(1474.594286,2170.18579,2292.385175,2331.810043,2446.303383,2345.896292,2039.941428,1837.222978,1863.173699,3183.799632,3226.102943,2774.14364)
pylor_e <- c(3514.158286,2290.26649,3960.291171,7944.077738,5382.691117,2042.018426,8740.645648,3676.13096,4873.884677,4772.275927,2698.953974,5967.189772)
enter_c <- c(893.5382746,845.6664571,1092.946477,999.2505779,924.3808297,910.3110938,812.2913175,865.7953789,717.7262783,1315.670413,604.7432671,919.0730522)
enter_e <- c(1468.744892,1996.830681,1353.123073,1803.25691,1249.578623,5069.040701,1393.362742,1003.490404,1751.552306,846.4962441,946.7284227,1535.48713)
flex_c <- c(98.06234955,87.48273694,154.5338242,112.9058501,150.9447779,107.6381361,132.0549212,163.3087092,91.08201501,162.327402,111.1615995,124.3581632)
flex_e <- c(129.9534651,137.9148303,174.0942866,231.4991979,221.2647732,5259.354761,221.5176204,182.2507489,159.924341,225.8766619,186.7125587,289.7709176)

df <- data.frame(jejun_c, jejun_e, coli_c,coli_e,pylor_c,pylor_e,enter_c, enter_e,flex_c, flex_e)
colnames <- c("C. jejuni\ncontrol", "C. jejuni\nICD", "C. coli\ncontrol", "C. coli\nICD",
                  "H. pylori\ncontrol", "H. pylori\nICD", "Y. enterocolitica\ncontrol", "Y. enterocolitica\nICD",
                  "S. flexneri\ncontrol", "S. flexneri\nICD")
colnames(df) <- colnames

# graphing
fd = as.data.frame(t(df))
fd$group <- row.names(fd)
fd.m <- melt(fd, id.vars="group")
fd.m$color <- rep(c("red", "green"), times=60)
pp = ggplot(fd.m, aes(group, value)) + 
  ylab("Normalized transcript abundance (log scale)") +
  geom_boxplot() + 
  scale_y_log10() +
  geom_jitter(aes(color=color), position=position_jitter(w=0.2,h=0.1)) +
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.title.x = element_blank())
  
df1 <- data.frame(a = c(1, 1,2, 2), b = c(60000, 71000, 71000, 60000))
df2 <- data.frame(a = c(3, 3,4, 4), b = c(64000, 76000, 76000, 64000))
df3 <- data.frame(a = c(5, 5,6, 6), b = c(10000, 12000, 12000, 10000))
df4 <- data.frame(a = c(7, 7,8, 8), b = c(650, 700, 700, 650))
df5 <- data.frame(a = c(9, 9,10, 10), b = c(6700, 7300, 7300, 6700))

pp + geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 74000, label = "*", size = 8) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 3.5, y = 75000, label = "*", size = 8) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 5.5, y = 13000, label = "*", size = 8) +
  geom_line(data = df4, aes(x = a, y = b)) + annotate("text", x = 7.5, y = 750, label = "*", size = 8) +
  geom_line(data = df5, aes(x = a, y = b)) + annotate("text", x = 9.5, y = 7600, label = "*", size = 8)

