## Submission

This is a new package (first submission to CRAN).

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

## Notes for the reviewers

* `standardize_taxa()` queries the GBIF name-matching API. Its example runs
  only when 'curl', 'jsonlite' and an internet connection are available, and
  a failure of the service gives an error of class `bmassa_error_network`,
  which the example turns into a message. The tests that need the internet
  are skipped on CRAN.
* The reference tables in `R/sysdata.rda` hold the parameters of
  hierarchical models fitted to published data sets and compilations, which
  `cite_sources()` lists with their citation, DOI and licence; the package
  does not include the data themselves.
