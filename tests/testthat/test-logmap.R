
test_set_mappedRowNames <- function(map, i, value) {
  # Add a new column and test whether 'map' automatically extends the rows
  old.rownames <- rownames(map)
  mappedRowNames(map, i, append = TRUE) <- value

  expect_equal(rownames(map), union(old.rownames, value))
  expect_equal_no_attr(map[, i], rownames(map) %in% value)
  expect_setequal(mappedRowNames(map, i), value)
  map
}

test_that("Create a sparse logical map from empty", {
  s1 <- sample(LETTERS, 5)
  s2 <- sample(LETTERS, 5)
  map1 <- sparseLogMap(character())

  # Add a new column
  map1 <- test_set_mappedRowNames(map1, 1, s1)
  map1 <- test_set_mappedRowNames(map1, "a", s2)
  expect_identical(rownames(map1), union(s1, s2))
  expect_identical(colnames(map1), c("", "a"))
})

test_that("Create a sparse logical map from known row names", {
  s1 <- LETTERS[5:1]
  map1 <- sparseLogMap(s1)

  # Get a non-existing column
  expect_error(mappedRowNames(map1, 1), "out-of-bounds")

  # Add a new column
  s2 <- LETTERS[3:8]
  map1 <- test_set_mappedRowNames(map1, "a", s2)

  # Get an existing column
  expect_equal(mappedRowNames(map1, 1), rownames(map1)[map1[, 1]])
  expect_equal(mappedRowNames(map1, "a"), rownames(map1)[map1[, "a"]])

  # Add a new column
  s3 <- sample(LETTERS, 5)
  map1 <- test_set_mappedRowNames(map1, 2, s3)
  expect_equal(mappedRowNames(map1, 2), rownames(map1)[map1[, 2]])
  expect_equal(mappedRowNames(map1, "a"), rownames(map1)[map1[, 1]])
  expect_identical(colnames(map1), c("a", ""))

  s4 <- sample(LETTERS, 5)
  map1 <- test_set_mappedRowNames(map1, "b", s4)
  expect_equal(mappedRowNames(map1, "b"), rownames(map1)[map1[, 3]])
  expect_identical(colnames(map1), c("a", "", "b"))
})

test_that("Create a new log map with initial columns", {
  s1 <- LETTERS[5:1]
  s2 <- LETTERS[3:8]

  map1 <- sparseLogMap(s1, 2)
  expect_equal_no_attr(map1[, 1], rep(FALSE, nrow(map1)))
  expect_equal_no_attr(map1[, 2], rep(FALSE, nrow(map1)))

  map1 <- test_set_mappedRowNames(map1, 1, s2)
  map1 <- test_set_mappedRowNames(map1, "a", sample(LETTERS, 5))
  map1 <- test_set_mappedRowNames(map1, 4, sample(LETTERS, 5))
  map1 <- test_set_mappedRowNames(map1, "a", sample(LETTERS, 5))
  expect_identical(colnames(map1), c("", "", "a", ""))

})

test_that("Get rows that are mapped to all coloumns", {
  # Get which rows are mapped to all columns
  s1 <- sample(LETTERS, 5)
  s2 <- sample(LETTERS, 5)
  s3 <- sample(LETTERS, 5)
  map1 <- sparseLogMap(character())
  map1 <- test_set_mappedRowNames(map1, "s1", s1)
  map1 <- test_set_mappedRowNames(map1, "s2", s2)
  map1 <- test_set_mappedRowNames(map1, "s3", s3)
  expect_equal_no_attr(
    intersectedRows(map1),
    which(rownames(map1) %in% Reduce(intersect, list(s1, s2, s3)))
  )
})

test_that("Drop empty rows or columns", {
  map1 <- sparseLogMap(LETTERS)
  s1 <- sample(LETTERS, 5)
  s2 <- sample(LETTERS, 5)
  s3 <- sample(LETTERS, 5)

  map1 <- test_set_mappedRowNames(map1, "s1", s1)
  map1 <- test_set_mappedRowNames(map1, "s2", s2)
  map1 <- test_set_mappedRowNames(map1, "s3", s3)

  # Drop empty rows
  map2 <- dropMapRows(map1)
  expect_equal(rownames(map2), rownames(map1)[rowSums(map1) > 0])

  # Drop empty columns
  map2 <- dropMapRows(map1)
  mappedRowNames(map2, 2) <- NULL
  map2 <- dropMapCols(map2)
  expect_equal(colnames(map2), c("s1", "s3"))

  map2 <- dropMapRows(map1)
  mappedRowNames(map2, 2) <- character()
  map2 <- dropMapCols(map2)
  expect_equal(colnames(map2), c("s1", "s3"))
})

test_that("Remove columns", {
  map1 <- sparseLogMap(LETTERS)
  s1 <- sample(LETTERS, 5)
  s2 <- sample(LETTERS, 5)
  s3 <- sample(LETTERS, 5)
  s4 <- sample(LETTERS, 5)

  map1 <- test_set_mappedRowNames(map1, "s1", s1)
  map1 <- test_set_mappedRowNames(map1, "s2", s2)
  map1 <- test_set_mappedRowNames(map1, "s3", s3)
  map1 <- test_set_mappedRowNames(map1, "s4", s4)
  map1 <- test_set_mappedRowNames(map1, "s5", s4)

  # Remove columns
  map2 <- removeMapCols(map1, c(1, 3))
  expect_error(removeMapCols(map2, c("s1", "s2")), "out of bounds")
  map2 <- removeMapCols(map2, c("s2"))
  expect_equal(colnames(map2), c("s4", "s5"))
})

test_that("sparseLogMap with append = FALSE", {
  s1 <- sample(LETTERS, 5)

  ## Cannot add new row names when append = FALSE
  map1 <- sparseLogMap(character())
  expect_error(mappedRowNames(map1, 1, append = FALSE) <- s1)

  ## Can add existing row names even when append = FALSE
  map1 <- sparseLogMap(s1, 1)
  s2 <- sample(s1, 3)
  mappedRowNames(map1, 2, append = FALSE) <- s2
  expect_true(setequal(mappedRowNames(map1, 2), s2))
})
