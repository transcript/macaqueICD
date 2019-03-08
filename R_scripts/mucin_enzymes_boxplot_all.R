# mucin_enzymes_boxplot_all.R
# last edited 11/29/2017
# Note: these counts are pulled from DESeq normalized counts for macaque functional annotations.
# yes, it's hard-coded.  No, I'm not thrilled about it.
# Second note: these are enzyme annotation counts across ALL the functional annotations, not just
#   those for specific mucin degrading organisms.

library(ggplot2)

a_fucos_c <- c(8446.650341,6491.853377,9402.966407,8243.382424,17640.13249,8099.400686,9480.427937,7510.500698,6746.647411,9587.05242,7542.214346,8836.429775)
a_fucos_e <- c(39932.73216,42412.95782,56194.56917,27686.09316,27273.72379,21824.94202,69381.70167,22894.22588,24913.5542,14926.06013,26001.55575,25388.79285)
n_acetyl_c <- c(3739.428371,3714.71787,9040.490271,9293.501075,25091.09704,5957.515866,9530.686715,7229.047213,5174.263679,11013.03713,5792.090218,4992.996475)
n_acetyl_e <- c(32426.66223,30547.51233,38485.33599,52413.35909,59158.52079,30235.84581,56714.86949,17471.42602,18125.75133,7700.84254,10602.24144,25589.38141)
a_mannos_c <- c(3818.714081,3046.672129,3883.89483,3892.32477,6605.790599,3805.176453,3693.001863,3095.199155,3066.650991,5828.06016,3592.586173,3490.004692)
a_mannos_e <- c(9860.464788,11322.67551,11119.04564,7275.659067,7538.340539,6658.209758,12017.79476,4902.383877,4608.533613,2992.055474,6070.99478,5357.236758)
a_galac_c <- c(16325.30071,9490.758167,12104.986,15609.27128,25418.78677,10882.93617,13542.01868,11103.83013,11085.57015,14544.5774,9906.725079,11746.37553)
a_galac_e <- c(39431.86095,35540.15323,121464.7511,50687.68021,42370.61412,27572.52984,55941.43238,51347.59729,44337.30481,21202.37046,54330.50281,29617.76924)
b_galac_c <- c(35687.65417,30816.20193,39317.17847,37546.8443,60744.5881,35259.78088,39369.06771,23629.30979,27719.69316,40244.79541,28861.42953,30271.10939)
b_galac_e <- c(154741.2109,144017.3017,143590.5895,81259.86004,92412.35107,65433.91754,178452.039,79810.71668,81881.58526,46512.78095,83972.97633,85362.39586)
a_sial_c <- c(5395.618759,2599.48304,1595.660775,3095.027138,1738.166448,2081.475369,3404.22685,857.1642005,2206.320003,2563.892024,911.2274913,2156.99984)
a_sial_e <- c(846.7825199,460.2325535,525.9041989,687.781121,882.0063559,3023.615105,175.4386833,2282.127535,617.9065325,900.7291839,1781.050905,899.1858786)
o_acet_c <- c(2353.3967,1445.917716,2885.016528,2491.109846,4651.317547,1641.075623,3291.783482,1840.65223,2103.238871,3111.633549,1611.628221,2294.333913)
o_acet_e <- c(11673.99847,14310.24349,14019.86547,10882.45967,7249.039962,4672.53008,18018.14867,6032.174465,5838.574734,2945.426126,6122.402765,8915.606444)

df <- data.frame(a_fucos_c, a_fucos_e, n_acetyl_c, n_acetyl_e, a_mannos_c, a_mannos_e,
                 a_galac_c, a_galac_e, b_galac_c, b_galac_e, a_sial_c, a_sial_e, o_acet_c, o_acet_e)
colnames <- c("Alpha-L-fucosidase\ncontrol", "Alpha-L-fucosidase\nICD", 
              "N-acetylglucosaminidase\ncontrol", "N-acetylglucosaminidase\nICD",
              "Alpha-mannosidase\ncontrol", "Alpha-mannosidase\nICD", 
              "Alpha-galactosidase\ncontrol", "Alpha-galactosidase\nICD",
              "Beta-galactosidase\ncontrol", "Beta-galactosidase\nICD", 
              "Alpha-sialidase\ncontrol", "Alpha-sialidase\nICD",
              "O-acetylesterase\ncontrol", "O-acetylesterase\nICD")
colnames(df) <- colnames

# ggplot
fd = as.data.frame(t(df))
fd$group <- row.names(fd)
fd.m <- melt(fd, id.vars="group")
fd.m$color <- rep(c("red", "green"), times=84)
pp = ggplot(fd.m, aes(group, value)) + geom_boxplot() + scale_y_log10() +
  geom_jitter(aes(color=color), position=position_jitter(w=0.2,h=0.1)) +
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.title.x = element_blank()) +
  ylab("Normalized transcript abundance (log scale)")
pp

df1 <- data.frame(a = c(1, 1,2, 2), b = c(158500, 200500, 200500, 158500))
df2 <- data.frame(a = c(3, 3,4, 4), b = c(150000, 190000, 190000, 150000))
df3 <- data.frame(a = c(5, 5,6, 6), b = c(30000, 37000, 37000, 30000))
df4 <- data.frame(a = c(7, 7,8, 8), b = c(9000, 11000, 11000, 9000))
df5 <- data.frame(a = c(9, 9,10, 10), b = c(240000, 300000, 300000, 240000))
df6 <- data.frame(a = c(11, 11,12, 12), b = c(107000, 137000, 137000, 107000))
df7 <- data.frame(a = c(13, 13,14, 14), b = c(33000, 43000, 43000, 33000))

pp + 
  geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 250000, label = "*", size = 8) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 3.5, y = 240000, label = "*", size = 8) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 5.5, y = 40000, label = "*", size = 8) +
  geom_line(data = df4, aes(x = a, y = b)) + annotate("text", x = 7.5, y = 13000, label = "*", size = 8) +
  geom_line(data = df5, aes(x = a, y = b)) + annotate("text", x = 9.5, y = 400000, label = "*", size = 8) +
  geom_line(data = df6, aes(x = a, y = b)) + annotate("text", x = 11.5, y = 150000, label = "*", size = 8) +
  geom_line(data = df7, aes(x = a, y = b)) + annotate("text", x = 13.5, y = 50000, label = "*", size = 8)
