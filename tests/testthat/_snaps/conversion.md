# print() reports the conversion

    Code
      print(cv)
    Output
      <bmassa_estimate> 40 rows: 35 estimated (12 from lengths, 23 without lengths), 5 not estimated
      Lowest level shared with the reference: species 25, genus 6, family 4
      Interval: prediction 95 %, for the mean individual mass of a population; bias_correction = "lognormal"
      Currency: carbon (ug); 36 rows converted from the currency of their models; habitat unknown
      Not estimated:
          1  abundance is missing
          2  amphibian larvae and juveniles have no reference masses without lengths
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
          2  amphibian larvae and juveniles have no reference masses without lengths
      Use as.data.frame() for the table and summarize_biomass() for totals.

