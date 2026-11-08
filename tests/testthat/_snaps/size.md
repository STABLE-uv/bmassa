# print

    Code
      print(summarize_size(f$est, breaks = c(1, 10), n_draws = 200, seed = 1))
    Output
      <bmassa_size> mass of the individuals (wet mass, mg) by site: 6 groups, 1 incomplete
      Individuals without a size: 0-88 % of each group (median 0 %)
      Approximate median and classes (individuals at the mean size of their taxon): 2 of 5 groups
      Interval: prediction 95 %, from 200 draws
      Classes: [0, 1) [1, 10) [10, Inf)
       site per mean_size median_size share_excluded approximation
         S1  m2      3.66        2.74          0.000          TRUE
         S2  m2      2.81        3.58          0.000         FALSE
         S3  m2      3.60        2.87          0.000         FALSE
         S4  m2      2.58        2.31          0.000          TRUE
         S5  m2  23329.02    20996.12          0.882         FALSE
         S6  m2        NA          NA             NA            NA
      Use as.data.frame() for all the indicators and x$classes for the size classes.

---

    Code
      print(summarize_size(f$est, size = "length"))
    Output
      <bmassa_size> length of the individuals (mm) by site: 6 groups, 1 incomplete
      Individuals without a size: 0-88 % of each group (median 20 %)
      Approximate median and classes (individuals at the mean size of their taxon): 1 of 4 groups
      Interval: prediction 95 % (lengths and abundances taken as known)
       site per mean_size median_size share_excluded approximation
         S1  m2      9.17       12.00          0.400          TRUE
         S2  m2      7.00        8.00          0.000         FALSE
         S3  m2      6.00        5.88          0.000         FALSE
         S4  m2        NA          NA          1.000            NA
         S5  m2     58.00       57.36          0.882         FALSE
         S6  m2        NA          NA             NA            NA
      Use as.data.frame() for all the indicators and x$classes for the size classes.

