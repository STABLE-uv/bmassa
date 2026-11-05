# input errors are classed and informative

    Code
      estimate_biomass(ab, limno_l)
    Condition
      Error in `estimate_biomass()`:
      ! `bias_correction` is required when lengths are used: "none" (median mass of the individuals of a given length) or "lognormal" (their mean, needed for biomass).

# print() summarises the estimate

    Code
      print(est)
    Output
      <bmassa_estimate> 40 rows: 35 estimated (13 from lengths, 22 without lengths), 5 not estimated
      Lowest level shared with the reference: species 25, genus 6, family 4
      Interval: prediction 95 %, for the mean individual mass of a population; bias_correction = "lognormal"
      Not estimated:
          1  abundance is missing
          2  amphibian larvae and juveniles have no reference masses without lengths
          2  fish need a length: the reference has no masses of fish without lengths (the mean mass of a population depends on its sizes)
      Use as.data.frame() for the table and summarize_biomass() for totals.

