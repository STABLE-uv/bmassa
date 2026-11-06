# print() summarises by broad group

    Code
      print(check_coverage(ex$taxa, ex$lengths))
    Output
      <bmassa_coverage> 20 taxa, 22 rows (14 without lengths, 8 with lengths)
        phytoplankton      typical species 4, genus 1; median factor 12
        zooplankton        typical species 3; median factor 12
        zooplankton        length  species 1; median factor 2.6
        macroinvertebrate  typical species 1, genus 1, family 1; median factor 14
        macroinvertebrate  length  species 2, family 1; median factor 3.1
        fish               typical species 1, genus 1; median factor 7.8
        fish               length  species 1, genus 2; median factor 1.4
        amphibian          typical none, not estimable 1; median factor -
        amphibian          length  species 1; median factor 1.7
      Factor: the 95 % prediction interval is the median divided and multiplied by it.

