# input errors are classed and informative

    Code
      estimate_biomass(ab, limno, ln5)
    Condition
      Error in `estimate_biomass()`:
      ! `bias_correction` is required when lengths are used (also the typical lengths of tadpoles without lengths): "none" (median mass of the individuals of a given length) or "lognormal" (their mean, needed for biomass).

# print() summarises the estimate

    Code
      print(est)
    Output
      <bmassa_estimate> 40 rows: 37 estimated (14 from lengths, 23 without lengths), 3 not estimated
      Lowest level shared with the reference: species 27, genus 6, family 4
      Interval: prediction 95 %, for the mean individual mass of a population; bias_correction = "lognormal"
      Currency: wet mass (mg); 30 rows converted from the currency of their models; habitat unknown
      Not estimated:
          1  abundance is missing
          2  no conversion from dry mass to wet mass for tadpoles (the reference has their carbon, nitrogen and phosphorus only)
      Use as.data.frame() for the table and summarize_biomass() for totals.

