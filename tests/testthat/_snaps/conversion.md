# print() reports the conversion

    Code
      print(cv)
    Output
      <bmassa_estimate> 40 rows: 37 estimated (12 from lengths, 2 from typical lengths, 23 without lengths), 3 not estimated
      Lowest level shared with the reference: species 25, genus 6, family 6
      Interval: prediction 95 %, for the mean individual mass of a population; bias_correction = "lognormal"
      Currency: carbon (ug); 38 rows converted from the currency of their models; habitat unknown
      Not estimated:
          1  abundance is missing
          2  the reference data have too few conversion factors for amphibians to carbon (by convention, they go only between wet and dry mass and to energy)
      Use as.data.frame() for the table and summarize_biomass() for totals.

---

    Code
      print(en)
    Output
      <bmassa_estimate> 40 rows: 37 estimated (14 from lengths, 23 without lengths), 3 not estimated
      Lowest level shared with the reference: species 27, genus 6, family 4
      Interval: prediction 95 %, for the mean individual mass of a population; bias_correction = "lognormal"
      Currency: energy (J); 38 rows converted from the currency of their models; habitat unknown
      Not estimated:
          1  abundance is missing
          2  no conversion from dry mass to energy for tadpoles (the reference has their carbon, nitrogen and phosphorus only)
      Use as.data.frame() for the table and summarize_biomass() for totals.

