test_that("example_survey has the three tables with consistent keys", {
  ex <- example_survey
  expect_named(ex, c("abundance", "taxa", "lengths"))
  expect_true(all(ex$abundance$taxon_id %in% ex$taxa$taxon_id))
  expect_true(all(ex$lengths$taxon_id %in% ex$taxa$taxon_id))
  expect_false(anyDuplicated(ex$taxa$taxon_id) > 0)
  expect_true(all(ex$taxa$broad_group %in% c("phytoplankton", "zooplankton", "macroinvertebrate",
                                             "fish", "amphibian")))
  expect_true(all(ex$abundance$abundance >= 0, na.rm = TRUE))
  # Phytoplankton rows give their counting unit.
  phy <- ex$abundance$taxon_id %in% ex$taxa$taxon_id[ex$taxa$broad_group == "phytoplankton"]
  expect_true(all(ex$abundance$count_unit[phy] == "cell"))
})
