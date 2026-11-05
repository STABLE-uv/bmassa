# print() reports the conversion

    Code
      print(cv)
    Output
      <bmassa_estimate> 40 rows: 33 estimated (11 from lengths, 22 without lengths), 7 not estimated
      Lowest level shared with the reference: species 23, genus 6, family 4
      Interval: prediction 95 %, for the mean individual mass of a population; bias_correction = "lognormal"
      Converted to carbon (ug): 34 rows; habitat unknown
      Not estimated:
          1  abundance is missing
          2  amphibian larvae and juveniles have no reference masses without lengths
          2  fish need a length: the reference has no masses of fish without lengths (the mean mass of a population depends on its sizes)
          2  the reference data have too few conversion factors for amphibians to carbon (by convention, they go only between wet and dry mass and to energy)
      Use as.data.frame() for the table and summarize_biomass() for totals.

---

    Code
      print(en)
    Output
      <bmassa_estimate> 40 rows: 25 estimated (13 from lengths, 12 without lengths), 15 not estimated
      Lowest level shared with the reference: species 17, genus 4, family 4
      Interval: prediction 95 %, for the mean individual mass of a population; bias_correction = "lognormal"
      Converted to energy (J): 26 rows; habitat unknown
      Not estimated:
          1  abundance is missing
          2  amphibian larvae and juveniles have no reference masses without lengths
          2  fish need a length: the reference has no masses of fish without lengths (the mean mass of a population depends on its sizes)
         10  no conversion from biovolume to energy in the reference
      Use as.data.frame() for the table and summarize_biomass() for totals.

