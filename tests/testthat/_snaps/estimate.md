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
      <bmassa_estimate> 40 rows: 37 estimated (13 from lengths, 24 without lengths), 3 not estimated
      Lowest level shared with the reference: species 27, genus 6, family 4
      Interval: prediction 95 %, for the mean individual mass of a population; bias_correction = "lognormal"
      Not estimated:
          1  abundance is missing
          2  amphibian larvae and juveniles have no reference masses without lengths
      Use as.data.frame() for the table and summarize_biomass() for totals.

