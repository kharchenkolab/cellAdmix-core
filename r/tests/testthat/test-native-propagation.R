test_that("vetoed calls require every rule of a combination to be rejected", {
  rules <- data.frame(
    factor = c(1L, 1L, 2L, 3L),
    source_cell_type = c("A", "B", "A", "A"),
    target_cell_type = c("T", "T", "T", "U"),
    keep = c(FALSE, TRUE, FALSE, NA),
    stringsAsFactors = FALSE
  )
  # Factor 1 on T keeps one surviving rule; factor 2 on T loses all of its
  # rules; the NA on factor 3 counts as kept (matching the correction).
  expect_equal(cellAdmixCore:::.celladmix_vetoed_calls(rules), "2_T")

  rules$keep <- c(FALSE, FALSE, FALSE, FALSE)
  expect_setequal(cellAdmixCore:::.celladmix_vetoed_calls(rules),
    c("1_T", "2_T", "3_U"))

  rules$keep <- TRUE
  expect_null(cellAdmixCore:::.celladmix_vetoed_calls(rules))
  expect_null(cellAdmixCore:::.celladmix_vetoed_calls(NULL))
  expect_null(cellAdmixCore:::.celladmix_vetoed_calls(rules[0, ]))
  expect_null(cellAdmixCore:::.celladmix_vetoed_calls(
    rules[, setdiff(names(rules), "keep")]))
})

test_that("heatmap accepts vetoed calls in both display branches", {
  mat <- matrix(c(4, 0.2, 3.5, 4.2), nrow = 2,
    dimnames = list(c("T", "S"), c("F1", "F2")))
  ann <- list(
    score_matrix = mat,
    source_calls = list(f_1 = "S", f_2 = "S"),
    remove_calls = c("1_T", "2_T"),
    threshold = 1
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- cellAdmixCore:::celladmix_plot_bridge_heatmap(ann,
    vetoed_calls = "2_T")
  expect_identical(out$remove_calls, ann$remove_calls)
  # And without any vetoes the caption line is dropped silently.
  expect_silent(cellAdmixCore:::celladmix_plot_bridge_heatmap(ann))
})

test_that("example selection drops factor/target combinations vetoed by the native check", {
  pairs <- data.frame(
    target_cell = c("c1", "c2"),
    target_cell_type = c("T", "T"),
    factor = c(1L, 2L),
    mean_score = c(0.6, 0.5),
    factor_count = c(5, 5),
    stringsAsFactors = FALSE
  )
  cell_data <- data.frame(
    cell_id = c("c1", "c2"),
    transcript_count = c(100, 100),
    dominant_factor = c(9L, 9L),
    stringsAsFactors = FALSE
  )
  score <- structure(list(
    pairs = function() pairs,
    fit = list(cell_factors = function() cell_data)
  ), class = "CellAdmixScore")
  rules <- data.frame(
    factor = c(1L, 2L),
    source_cell_type = c("S", "S"),
    target_cell_type = c("T", "T"),
    keep = c(TRUE, FALSE),
    native_check = c("pass", "native_median"),
    stringsAsFactors = FALSE
  )
  out <- celladmix_select_example_cells(score,
    rules = rules,
    score_annotation = list(source_calls = list(f_1 = "S", f_2 = "S")),
    cell_data = cell_data,
    min_molecules = 0, min_factor_molecules = 0)
  # Factor 2's rule is vetoed, so c2 (whose only evidence is factor 2) is
  # not presented as an admixture example; c1 survives on factor 1.
  expect_identical(sort(as.character(out$target_cell)), "c1")
  expect_identical(as.integer(out$top_admix_factor), 1L)

  rules$keep <- c(TRUE, TRUE)
  out_all <- celladmix_select_example_cells(score,
    rules = rules,
    score_annotation = list(source_calls = list(f_1 = "S", f_2 = "S")),
    cell_data = cell_data,
    min_molecules = 0, min_factor_molecules = 0)
  expect_setequal(as.character(out_all$target_cell), c("c1", "c2"))
})
