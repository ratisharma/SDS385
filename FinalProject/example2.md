





Title
=====

Introduction
------------

$\\sum\_{i=1}^n X\_i$

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

The Data
--------

The data are produced by means of

[Link to
paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004400)

Normalizing RNA and protein counts
----------------------------------

RNAseq data must be normalized to read depth, depending on

First, read in the data from the csv files.

    # The as.is = TRUE option will preserve strings in the csv file instead of converting them to factors
    rna <- read.csv("rna.csv", header = TRUE, as.is = TRUE)
    pro <- read.csv("protein.csv", header = TRUE, as.is = TRUE)

    head(rna)

    ##     gene_id t3_1 t4_1 t5_1 t6_1 t8_1 t24_1 t48_1 t168_1 t336_1 t3_2 t4_2
    ## 1 ECB_00001  287  628  289  174   14     3     4    151    308   62   82
    ## 2 ECB_00002 3948 4251 3669 2630   19     1     1    129    670  995 1794
    ## 3 ECB_00003 1434 1669 1444 1359    7     1     1     51    183  468  620
    ## 4 ECB_00004 1483 2243 1565 1123    4     0     3     13     99  436  752
    ## 5 ECB_00005  154  248  174  149    2     0     0      8     27   33   62
    ## 6 ECB_00006   94  144   30   15    1     1     0      7     18   23   26
    ##   t5_2 t6_2 t8_2 t24_2 t48_2 t168_2 t336_2 t3_3 t4_3 t5_3 t6_3 t8_3 t24_3
    ## 1  709   88   45    20     7     65     92  196  180  180  107   93    62
    ## 2 1204 2216   26    13     4     72    147 3421 3415 3500 4236 3496    26
    ## 3  226  727   25     4     1     22     46 1558 1229 1337 1604 1131     7
    ## 4  243  796   22     4     0      4     37 1594 1602 1563 1573 1089     1
    ## 5    8   70    3     1     1      3      6  102  125  139   96   71     2
    ## 6    5   16    3     0     1      0      4   14   14   20   30   17     2
    ##   t48_3 t168_3 t336_3
    ## 1     7    157    131
    ## 2     4    289    230
    ## 3     3    115    135
    ## 4     0     30     41
    ## 5     0     28     13
    ## 6     0      9     18

    head(pro)

    ##          gene_id t3_1 t4_1 t5_1 t6_1 t8_1 t24_1 t48_1 t168_1 t336_1 t3_2
    ## 1 YP_003044966.1    0    0    0    0    0     0     1      1      0    0
    ## 2 YP_003043986.1   87  101  104  133  140   132   127    130    128  104
    ## 3 YP_003045405.1    2    4    5    5    6     3     4      5      6    8
    ## 4 YP_003044440.1    0    0    0    0    0     0     0      0      0    0
    ## 5 YP_003043527.1    2    4    0    0    7     5     6      5      5    0
    ## 6 YP_003045246.1    6    3    3    1    0     0     3      2      0    0
    ##   t4_2 t5_2 t6_2 t8_2 t24_2 t48_2 t168_2 t336_2 t3_3 t4_3 t5_3 t6_3 t8_3
    ## 1    0    2    0    0     0     0      1      0  0.0  0.0  0.0  0.5    0
    ## 2  135  119  143  168   169   139    112    123 67.0 64.0 55.0 73.5   73
    ## 3    8    5    4    3     0     0      2      2  2.5  1.0  4.0  2.0    2
    ## 4    0    0    0    0     0     0      0      0  0.0  0.0  2.0  0.0    0
    ## 5    3    3    0    8     7     8      6      0  0.0  0.0  0.5  0.5    7
    ## 6    1    1    1    1     0     2      4      0  6.0  4.5 10.5  7.5    7
    ##   t24_3 t48_3 t168_3 t336_3
    ## 1   0.5   0.0    0.5    0.0
    ## 2  89.5  56.5   55.5   39.5
    ## 3   1.0   2.5    2.0    0.0
    ## 4   0.0   0.0    0.0    0.0
    ## 5   4.5   0.0    0.0    0.0
    ## 6  11.0   4.0    6.5    7.0

Let's first store the names of our

    rna_names <- rna[, 1]

We will normalize the RNA counts by the results from DEseq, a method
described here. You can find a script for producing these scaled size
factors here in the appendix.

    ssf <- read.csv("ssf.csv", header = T)

    ## Warning in read.table(file = file, header = header, sep = sep, quote =
    ## quote, : incomplete final line found by readTableHeader on 'ssf.csv'

    ssf <- unlist(ssf[1, 2:ncol(ssf)])

    rna_norm <- sweep(as.matrix(rna[, -1]), 2, ssf, `/`)
    rna_norm[1, ]

    ##      t3_1      t4_1      t5_1      t6_1      t8_1     t24_1     t48_1 
    ##  287.0000  241.4178  211.4464  210.2930  161.1804  132.0854  219.7473 
    ##    t168_1    t336_1      t3_2      t4_2      t5_2      t6_2      t8_2 
    ##  412.6012  447.6265   62.0000  141.6106 4031.9200  195.1336  274.9117 
    ##     t24_2     t48_2    t168_2    t336_2      t3_3      t4_3      t5_3 
    ##  385.3193  213.2850  683.9320  328.2673  196.0000  171.2717  176.6082 
    ##      t6_3      t8_3     t24_3     t48_3    t168_3    t336_3 
    ##  139.9926  123.2953  420.8323  179.0655  408.5317  231.2244

Now we find the average of the adjusted counts across the three
biological replicates.

    rna_norm_av <- matrix(rep(0, nrow(rna) * 9), nrow = nrow(rna))
    rna_norm_sd <- matrix(rep(0, nrow(rna) * 9), nrow = nrow(rna))

    for (t in 1:9) {
        rna_norm_t <- rna_norm[, c(t, t + 9, t + 18)]
        rna_norm_av[, t] <- apply(rna_norm_t, 1, mean)
        rna_norm_sd[, t] <- apply(rna_norm_t, 1, sd)
    }

    head(rna_norm_av)

    ##            [,1]       [,2]       [,3]      [,4]       [,5]      [,6]
    ## [1,]  181.66667  184.76673 1473.32486  181.8064  186.46249 312.74564
    ## [2,] 2788.00000 2660.58425 4321.77889 4544.8410 1670.80902 156.98801
    ## [3,] 1153.33333  960.57407 1217.83882 1784.3704  577.58314  56.20188
    ## [4,] 1171.00000 1228.41761 1353.48745 1726.7767  541.40051  27.95049
    ## [5,]   96.33333  107.11571  103.06057  153.6331   45.16064  10.94707
    ## [6,]   43.66667   37.85967   23.33548   30.9526   17.45940  19.20123
    ##           [,7]      [,8]      [,9]
    ## [1,] 204.03259 501.68829 335.70603
    ## [2,]  93.04570 620.69462 634.73777
    ## [3,]  54.04948 223.36100 222.79272
    ## [4,]  54.93682  51.89115 116.08946
    ## [5,]  10.15643  42.09497  27.86488
    ## [6,]  10.15643  14.18207  24.06792

    head(rna_norm_sd)

    ##            [,1]      [,2]        [,3]       [,4]        [,5]       [,6]
    ## [1,]  113.18274  51.25384 2215.876880   36.99658   78.906691 157.460771
    ## [2,] 1574.98222 892.09999 2218.683078 1224.22439 2567.104635 104.585551
    ## [3,]  596.74562 280.61027  140.354647  272.53800  799.157606  18.150832
    ## [4,]  638.94366 336.57347  195.809407  351.95919  782.704044  42.668607
    ## [5,]   60.69871  11.80092   50.059994   27.27346   42.472602   9.898222
    ## [6,]   43.82161  21.88463    4.565962   11.26478    5.563508  22.546953
    ##          [,7]      [,8]       [,9]
    ## [1,] 21.86226 157.84079 108.392658
    ## [2,] 34.42100 232.29124 299.502331
    ## [3,] 23.14928  80.25275  52.651176
    ## [4,] 95.15336  22.90235  38.325458
    ## [5,] 17.59145  27.08099   9.881064
    ## [6,] 17.59145  12.46809   8.935017

Finally, we normalize the averages and standard deviation values by the
maximum value of the averages across the entire duration in order to put
averages on the scale of 0 to 1, and then save these two data frames to
csv's we can refer to.

    rna_norm_av_max <- apply(rna_norm_av, 1, max)

    rna_norm_av <- sweep(rna_norm_av, 1, rna_norm_av_max, '/') 
    rna_norm_sd <- sweep(rna_norm_sd, 1, rna_norm_av_max, '/') 

    rna_norm_av <- data.frame(rna_names, rna_norm_av)
    colnames(rna_norm_av) <- c("gene_id", "t3", "t4", "t5", "t6", 
                               "t8", "t24", "t48", "t168", "t336")
    head(rna_norm_av)                          

    ##     gene_id        t3        t4        t5        t6        t8        t24
    ## 1 ECB_00001 0.1233039 0.1254080 1.0000000 0.1233987 0.1265590 0.21227202
    ## 2 ECB_00002 0.6134428 0.5854076 0.9509197 1.0000000 0.3676276 0.03454203
    ## 3 ECB_00003 0.6463531 0.5383266 0.6825034 1.0000000 0.3236902 0.03149676
    ## 4 ECB_00004 0.6781421 0.7113934 0.7838231 1.0000000 0.3135324 0.01618651
    ## 5 ECB_00005 0.6270351 0.6972177 0.6708228 1.0000000 0.2939513 0.07125462
    ## 6 ECB_00006 1.0000000 0.8670154 0.5344002 0.7088382 0.3998335 0.43972283
    ##          t48       t168       t336
    ## 1 0.13848446 0.34051437 0.22785608
    ## 2 0.02047282 0.13657125 0.13966116
    ## 3 0.03029051 0.12517636 0.12485789
    ## 4 0.03181466 0.03005088 0.06722899
    ## 5 0.06610836 0.27399678 0.18137293
    ## 6 0.23259002 0.32478031 0.55117377

    write.csv(rna_norm_av, file = "Data2/rna_norm_av.csv", row.names = FALSE)

    rna_norm_sd <- data.frame(rna_names, rna_norm_sd)
    colnames(rna_norm_sd) <- colnames(rna_norm_av)
    head(rna_norm_sd)

    ##     gene_id         t3         t4        t5         t6         t8
    ## 1 ECB_00001 0.07682131 0.03478788 1.5039975 0.02511095 0.05355689
    ## 2 ECB_00002 0.34654287 0.19628849 0.4881762 0.26936573 0.56483926
    ## 3 ECB_00003 0.33442923 0.15726010 0.0786578 0.15273622 0.44786531
    ## 4 ECB_00004 0.37002100 0.19491429 0.1133959 0.20382438 0.45327461
    ## 5 ECB_00005 0.39508880 0.07681238 0.3258412 0.17752338 0.27645479
    ## 6 ECB_00006 1.00354828 0.50117473 0.1045640 0.25797195 0.12740859
    ##          t24         t48       t168       t336
    ## 1 0.10687444 0.014838723 0.10713238 0.07357010
    ## 2 0.02301193 0.007573643 0.05111097 0.06589941
    ## 3 0.01017212 0.012973361 0.04497539 0.02950686
    ## 4 0.02470997 0.055104610 0.01326306 0.02219480
    ## 5 0.06442767 0.114503033 0.17627054 0.06431599
    ## 6 0.51634244 0.402857724 0.28552882 0.20461872

    write.csv(rna_norm_sd, file = "Data2/rna_norm_sd.csv", row.names = FALSE)

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
