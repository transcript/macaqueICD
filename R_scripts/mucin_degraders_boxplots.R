# mucin_degraders_boxplots.R
# last edited 11/29/2017
# Note: these counts are pulled from DESeq normalized organism counts for macaque files.
# yes, it's hard-coded.  No, I'm not thrilled about it.

library(ggplot2)

a_mucinophila_c <- c(10536.57992,6671.70978,8361.641069,10237.2478,14879.91165,7578.954931,8769.36808,6941.916242,7367.44203,11695.75752,8362.49380,7989.092176)
a_mucinophila_e <- c(9119.755152,7886.087368,8711.791335,8208.64825,9356.946851,8852.175593,11603.46497,9852.653292,6666.777969,4620.007634,8076.51504,8714.038831)
b_caccae_c <- c(105.3804353,169.5937269,115.2998481,145.357831,313.1168532,235.5737493,161.9976068,63.50894247,389.1023681,86.62008424,181.2417384,129.509093)
b_caccae_e <- c(364.1404388,283.1655557,295.8187472,236.7209843,351.7542548,347.1945686,380.4691954,191.1410294,233.9029477,78.6599238,286.053215,235.9989947)
b_fragilis_c <- c(225948.0941,181921.1189,233499.8077,254686.5264,432814.8146,227687.2568,229880.7461,125314.8084,176455.0093,247051.393,174350.9275,160155.6539)
b_fragilis_e <- c(746582.6572,795321.0809,803126.6676,795171.8952,634882.3892,564643.8128,1059899.248,448121.253,472179.3365,186193.0912,430149.829,580618.7673)
b_thetaiotaomicron_c <- c(190266.5714,160157.8696,198086.7405,216828.5863,366804.5425,180909.5681,189350.634,100813.3176,144990.1807,198621.2173,146919.3862,140209.7818)
b_thetaiotaomicron_e <- c(859552.8288,877028.282,871242.8265,628542.95,557195.76,518664.9646,1121225.809,366970.7732,548214.1133,160978.6165,482157.6548,521619.0186)
b_bifidum_c <- c(1499.475778,1802.604816,2010.541102,1775.393787,2355.237526,2133.695395,2037.638144,1919.525384,1911.265003,1523.012977,1303.128099,1704.957776)
b_bifidum_e <- c(2566.580936,2674.667399,2113.193251,3234.026389,2679.289465,2350.63582,2891.565885,2096.994898,2311.831459,1800.518439,1736.666171,2483.964103)
b_breve_c <- c(1629.005896,1835.602691,2157.868685,1845.03033,1787.635262,2323.138514,2152.80232,2235.773995,1842.771328,1636.914978,1577.407263,1509.958289)
b_breve_e <- c(3117.529481,2508.876167,1988.637989,7451.489219,1726.148904,3693.121485,2057.915603,1836.954195,2282.457601,2504.8495,1868.322463,2274.851069)
b_longum_c <- c(3881.512702,3991.975417,4961.896936,4228.898758,4570.757571,5089.131075,4073.740766,4134.561765,4079.745616,3264.963333,3681.019707,4164.158849)
b_longum_e <- c(10760.41765,10970.09772,5412.492293,13898.65485,5306.099465,6336.943829,6344.535204,5415.29207,5861.71666,6416.91782,4996.954696,4902.206967)
r_gnavus_c <- c(683.5092125,613.9139434,1187.428297,744.3673114,3411.85097,944.7552974,1654.525321,1522.918518,1154.91995,311.6958937,465.7912676,873.4505308)
r_gnavus_e <- c(3588.610792,6386.630384,15913.35012,10621.11358,6976.932175,3950.302647,6210.948243,6153.185346,8685.958763,4861.03896,8570.82457,3861.122794)
r_torques_c <- c(50249.63606,75207.52764,91407.15737,68574.41601,160123.7173,87332.04795,88737.83599,74847.23285,77482.37718,49056.29574,52847.67436,63771.45476)
r_torques_e <- c(158541.8738,181855.3755,156935.3839,196911.8252,151458.574,160077.2706,148349.1667,134495.4963,223629.7126,153712.3159,160726.0005,150888.4965)


df <- data.frame(a_mucinophila_c, a_mucinophila_e, b_caccae_c, b_caccae_e, b_fragilis_c, b_fragilis_e, 
                 b_thetaiotaomicron_c, b_thetaiotaomicron_e, b_bifidum_c, b_bifidum_e, b_breve_c, b_breve_e,
                 b_longum_c, b_longum_e, r_gnavus_c, r_gnavus_e, r_torques_c, r_torques_e)
colnames <- c("A. mucinophila\ncontrol", "A. mucinophila\nICD", "B. caccae\ncontrol", "B. caccae\nICD",
                  "B. fragilis\ncontrol", "B. fragilis\nICD", "B. thetaiotaomicron\ncontrol", "B. thetaiotaomicron\nICD",
              "B. bifidum\ncontrol", "B. bifidum\nICD", "B. breve\ncontrol", "B. breve\nICD",
              "B. longum\ncontrol", "B. longum\nICD", "R. gnavus\ncontrol", "R. gnavus\nICD",
              "R. torques\ncontrol", "R. torques\nICD")
colnames(df) <- colnames

# graphing with ggplot
fd = as.data.frame(t(df))
fd$group <- row.names(fd)
fd.m <- melt(fd, id.vars="group")
fd.m$color <- rep(c("red", "green"), times=108)
pp = ggplot(fd.m, aes(group, value)) + geom_boxplot() + scale_y_log10() +
  geom_jitter(aes(color=color), position=position_jitter(w=0.2,h=0.1)) +
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.title.x = element_blank()) +
  ylab("Normalized transcript abundance (log scale)")
pp

df2 <- data.frame(a = c(3, 3,4, 4), b = c(8500, 10500, 10500, 8500))
df3 <- data.frame(a = c(5, 5,6, 6), b = c(10000, 12000, 12000, 10000))
df4 <- data.frame(a = c(7, 7,8, 8), b = c(650, 800, 800, 650))
df5 <- data.frame(a = c(9, 9,10, 10), b = c(1717000, 2223000, 2223000, 1717000))
df6 <- data.frame(a = c(11, 11,12, 12), b = c(24700, 30300, 30300, 24700))
df7 <- data.frame(a = c(13, 13,14, 14), b = c(1717000, 2223000, 2223000, 1717000))
df8 <- data.frame(a = c(15, 15,16, 16), b = c(47000, 55000, 55000, 47000))
df9 <- data.frame(a = c(17, 17,18, 18), b = c(420000, 540000, 540000, 420000))

pp + 
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 3.5, y = 12000, label = "*", size = 8) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 5.5, y = 14000, label = "*", size = 8) +
  geom_line(data = df4, aes(x = a, y = b)) + annotate("text", x = 7.5, y = 850, label = "*", size = 8) +
  geom_line(data = df5, aes(x = a, y = b)) + annotate("text", x = 9.5, y = 2400000, label = "*", size = 8) +
  geom_line(data = df6, aes(x = a, y = b)) + annotate("text", x = 11.5, y = 35000, label = "*", size = 8) +
  geom_line(data = df7, aes(x = a, y = b)) + annotate("text", x = 13.5, y = 2400000, label = "*", size = 8) +
  geom_line(data = df8, aes(x = a, y = b)) + annotate("text", x = 15.5, y = 60000, label = "*", size = 8) +
  geom_line(data = df9, aes(x = a, y = b)) + annotate("text", x = 17.5, y = 600000, label = "*", size = 8) 


