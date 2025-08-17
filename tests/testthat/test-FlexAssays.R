
generate_sparse_matrix <- function(nrow, ncol, density = 0.75) {
  mat1 <- Matrix::rsparsematrix(nrow = nrow, ncol = ncol, density = density)
  rownames(mat1) <- sample(LETTERS, nrow(mat1))
  colnames(mat1) <- sample(letters, ncol(mat1))
  return(mat1)
}

test_add_mat <- function(mats, mat, i) {

  mats[[i]] <- mat

  rn <- mappedRowNames(rowMap(mats), i)
  cn <- mappedRowNames(colMap(mats), i)
  expect_equal(mats[[i]], mat[rn, cn])

  expect_in(rownames(mat), rownames(mats))
  expect_in(colnames(mat), colnames(mats))

  mat.list <- assays(mats)
  expect_equal(mats[[i]], mat.list[[i]])

  for (i in seq_along(mat.list)) {
    expect_equal(mat.list[[i]], assay(mats, i))
  }

  mat.list <- assays(mats, withDimnames = FALSE)
  for (i in seq_along(mat.list)) {
    expect_equal(mat.list[[i]], unname(assay(mats, i)))
    expect_equal(mat.list[[i]], assay(mats, i, withDimnames = FALSE))
  }
  mats
}

test_empty_fa <- function(input = list()) {
  mats <- FlexAssays(input)

  mat1 <- generate_sparse_matrix(8, 6)
  mats[[1]] <- mat1
  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(colnames(mats), colnames(mat1))
  expect_equal(dimnames(mats), dimnames(mat1))
  expect_named(mats, NULL)
  expect_length(mats, 1)
  invisible(NULL)
}

test_that("Create a FlexAssays with empty list", {
  test_empty_fa(list())
  test_empty_fa(List())
  test_empty_fa(NULL)
  test_empty_fa(SimpleList())
})

test_that("Create a FlexAssays with only one matrix", {
  mat1 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(mat1)
  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(colnames(mats), colnames(mat1))
  expect_equal(dimnames(mats), dimnames(mat1))
  expect_named(mats, NULL)
  expect_length(mats, 1)
})

test_that("Create a FlexAssays with 'rownames' and 'colnames'", {
  mat1 <- generate_sparse_matrix(8, 6)

  ## Cannot use 'rownames' and 'colnames' when use 'assays'
  expect_error(mats <- FlexAssays(mat1, rownames(mat1), colnames(mat1)))

  ## Create an empty FlexAssays
  mats <- FlexAssays(rownames = rownames(mat1), colnames = colnames(mat1))
  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(colnames(mats), colnames(mat1))
  expect_equal(dimnames(mats), dimnames(mat1))
  expect_length(mats, 0)

  ## Add `mat1`
  mats <- test_add_mat(mats, mat1, "mat1")
  expect_length(mats, 1)
  expect_named(mats, "mat1")

  ## Add subset of `mat1`
  mat1 <- mat1[sample(rownames(mat1), 6), sample(colnames(mat1), 4)]
  mats <- test_add_mat(mats, mat1, 1)
})

test_that("Create a FlexAssays with a list of matrices", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list)

  expect_length(mats, length(mat.list))
  expect_named(mats, NULL)
  expect_equal(rownames(mats), Reduce(union, lapply(mat.list, rownames)))
  expect_equal(colnames(mats), Reduce(union, lapply(mat.list, colnames)))

  mat.list <- list(m1 = mat1, mat2, m3 = mat3)
  mats <- FlexAssays(mat.list)
  expect_named(mats, names(mat.list))
  expect_equal(names(mats)[2], "")

  ## Cannot use 'rownames' and 'colnames' when use 'assays'
  expect_error(mats <- FlexAssays(mat.list, rownames(mat1), colnames(mat1)))
})

test_that("Create a FlexAssays with a list, and rowFlex = 'fixed'", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  rownames(mat1) <- rownames(mat2) <- LETTERS[1:8]

  # Wrong row names (error)
  rownames(mat2)[1:4] <- setdiff(LETTERS, rownames(mat1))[1:4]
  mat.list <- list(mat1, mat2)
  expect_error(mats <- FlexAssays(mat.list, rowFlex = "fixed"))

  # Missed rows (error)
  rownames(mat2) <- sample(rownames(mat1))
  mat.list <- list(mat1, mat2[1:4, ])
  expect_error(mats <- FlexAssays(mat.list, rowFlex = "fixed"))


  mat.list <- list(mat1, mat2)
  mats <- FlexAssays(mat.list, rowFlex = "fixed")
  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(mat2[rownames(mats), colnames(mats[[2]])], mats[[2]])
})

test_that("Create a FlexAssays with a list, and colFlex = 'fixed'", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)

  # Wrong column names (error)
  bad.names <- setdiff(letters, colnames(mat1))
  colnames(mat2) <- c(colnames(mat1)[1:4], bad.names[1:2])
  mat.list <- list(mat1, mat2)
  expect_error(mats <- FlexAssays(mat.list, colFlex = "fixed"))

  # Missed columns (error)
  colnames(mat2) <- sample(colnames(mat1))
  mat.list <- list(mat1, mat2[, 1:4])
  expect_error(mats <- FlexAssays(mat.list, colFlex = "fixed"))

  mat.list <- list(mat1, mat2)
  mats <- FlexAssays(mat.list, colFlex = "fixed")
  expect_equal(colnames(mats), colnames(mat1))
  expect_equal(mat2[rownames(mats[[2]]), colnames(mats)], mats[[2]])
})

test_that("Create a FlexAssays with duplicated dimnames should fail.", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  colnames(mat1)[2] <- colnames(mat1)[1]
  expect_error(FlexAssays(list(mat1, mat2)))

  mat1 <- generate_sparse_matrix(8, 6)
  rownames(mat1)[2] <- rownames(mat1)[1]
  expect_error(FlexAssays(list(mat1, mat2)))
})

test_that("FlexAssays .DollarNames", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1 = mat1, mat2, m3 = mat3)
  mats <- FlexAssays(mat.list)

  ## .DollarNames should get all column names of colData(x)
  dollar.names <- .DollarNames(mats)
  expect_identical(dollar.names, names(assays(mats, withDimnames = FALSE)))

  ## Use `$` to get or set one assay
  expect_identical(mats$mat1, mats[[1]])

  mats$mat2 <- mat2
  rn <- mappedRowNames(rowMap(mats), 2)
  cn <- mappedRowNames(colMap(mats), 2)
  expect_identical(mats$mat2, mat2[rn, cn, drop = FALSE])
})

test_that("'assays' for FlexAssays", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list)

  ## Size, names and S4 class
  expect_s4_class(assays(mats), "SimpleList")
  expect_length(assays(mats), length(mats))
  expect_named(assays(mats), names(mats))

  ## Values are correct
  for (i in seq_along(mat.list)) {
    rn <- mappedRowNames(rowMap(mats), i)
    cn <- mappedRowNames(colMap(mats), i)
    expect_equal(assays(mats)[[i]], mat.list[[i]][rn, cn, drop = FALSE])
  }

  ## `withDimnames = FALSE` works
  for (i in seq_along(mat.list)) {
    rn <- mappedRowNames(rowMap(mats), i)
    cn <- mappedRowNames(colMap(mats), i)
    expect_equal(
      assays(mats, withDimnames = FALSE)[[i]],
      unname(mat.list[[i]][rn, cn, drop = FALSE])
    )
  }

  ## `assay(x, i)` and `assays(x)[[i]]` are equivalences
  for (i in seq_along(mats)) {
    expect_equal(assays(mats)[[i]], assay(mats, i))
    expect_equal(assay(mats, i), mats[[i]])
    expect_equal(
      assays(mats, withDimnames = FALSE)[[i]],
      assay(mats, i, withDimnames = FALSE)
    )
  }
})

test_that("`assays<-` for FlexAssays", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list)

  # Set `assays<-` will automatically try to re-ordered the rows and columns
  mats2 <- mats
  rn <- sample(rownames(mat.list[[1]]))
  cn <- sample(colnames(mat.list[[1]]))
  mat.list[[1]] <- mat.list[[1]][rn, cn]
  assays(mats2) <- mat.list
  expect_equal(assays(mats2), assays(mats))

  # Subset some assays
  mat.list <- assays(mats)
  rn <- rownames(mat.list[[1]])[1:5]
  cn <- colnames(mat.list[[1]])[1:5]
  mat.list[[1]] <- mat.list[[1]][rn, cn]
  assays(mats2) <- mat.list
  expect_equal(mats2[[2]], mats[[2]])
  expect_equal(assay(mats2, 3), assay(mats, 3))
  expect_equal(assay(mats2), assay(mats)[rn, cn])

  # Extend new dim names for some assays
  mats2 <- mats
  mat.list <- assays(mats)
  new.names <- tolower(rownames(mat.list[[1]])[1:3])
  rownames(mat.list[[1]])[1:3] <- new.names
  assays(mats2) <- mat.list
  expect_setequal(rownames(mats2), c(rownames(mats), new.names))
  expect_equal(assay(mats2, 2), assay(mats, 2))
  expect_equal(assay(mats2, 3), assay(mats, 3))
  expect_equal(assay(mats2)[rownames(mat.list[[1]]), ], mat.list[[1]])
})

test_that("`assays<-` for FlexAssays with `rowFlex = 'bounded'`", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list, rowFlex = 'bounded')

  # Set `assays<-` will automatically try to re-ordered the rows and columns
  mats2 <- mats
  rn <- sample(rownames(mat.list[[1]]))
  cn <- sample(colnames(mat.list[[1]]))
  mat.list[[1]] <- mat.list[[1]][rn, cn]
  assays(mats2) <- mat.list
  expect_equal(assays(mats2), assays(mats))

  # Subset some assays
  mat.list <- assays(mats)
  rn <- rownames(mat.list[[1]])[1:5]
  cn <- colnames(mat.list[[1]])[1:5]
  mat.list[[1]] <- mat.list[[1]][rn, cn]
  assays(mats2) <- mat.list
  expect_equal(assay(mats2, 2), assay(mats, 2))
  expect_equal(assay(mats2, 3), assay(mats, 3))
  expect_equal(assay(mats2), assay(mats)[rn, cn])

  # Extend new dim names for some assays (error)
  mats2 <- mats
  mat.list <- assays(mats)
  rownames(mat.list[[1]])[1:3] <- tolower(rownames(mat.list[[1]])[1:3])
  expect_error(assays(mats2) <- mat.list)
})

test_that("`assays<-` for FlexAssays with `rowFlex = 'fixed'`", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)
  rownames(mat1) <- rownames(mat2) <- rownames(mat3) <- LETTERS[1:8]

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list, rowFlex = 'fixed')

  # Set `assays<-` will automatically try to re-ordered the rows and columns
  mats2 <- mats
  rn <- sample(rownames(mat.list[[1]]))
  cn <- sample(colnames(mat.list[[1]]))
  mat.list[[1]] <- mat.list[[1]][rn, cn]
  assays(mats2) <- mat.list
  expect_equal(assays(mats2), assays(mats))

  # Subset some assays (error)
  mat.list <- assays(mats)
  rn <- rownames(mat.list[[1]])[1:5]
  cn <- colnames(mat.list[[1]])[1:5]
  mat.list[[1]] <- mat.list[[1]][rn, cn]
  expect_error(assays(mats2) <- mat.list)

  # Extend new dim names for some assays (error)
  mats2 <- mats
  mat.list <- assays(mats)
  rownames(mat.list[[1]])[1:3] <- tolower(rownames(mat.list[[1]])[1:3])
  expect_error(assays(mats2) <- mat.list)
})

is_identical_dims <- function(fa1, fa2) {
  expect_equal(dimnames(fa1), dimnames(fa2))
  expect_equal(rowFlex(fa1), rowFlex(fa2))
  expect_equal(colFlex(fa1), colFlex(fa2))
  expect_equal(assayClasses(fa1), assayClasses(fa2))
}

test_that("`assasys(x) <- NULL` works", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)
  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list)

  mats2 <- mats
  assays(mats2) <- NULL
  is_identical_dims(mats2, mats)
  expect_length(mats2, 0)

  # rowFlex = "bounded"
  mats <- FlexAssays(mat.list, rowFlex = "bounded")
  mats2 <- mats
  assays(mats2) <- NULL
  is_identical_dims(mats2, mats)
  expect_length(mats2, 0)

  # rowFlex = "fixed"
  rownames(mat1) <- rownames(mat2) <- rownames(mat3) <- LETTERS[1:8]
  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list, rowFlex = "fixed")
  mats2 <- mats
  assays(mats2) <- NULL
  is_identical_dims(mats2, mats)
  expect_length(mats2, 0)
})

test_that("rowFlex = 'free', colFlex = 'free'", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat1)
  mats <- test_add_mat(mats, mat2, 2)
  mats <- test_add_mat(mats, mat3, "m3")
  expect_length(mats, 3)
})

test_that("rowFlex = 'fixed', colFlex = 'free'", {
  mat1 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(mat1, rowFlex = "fixed")

  mat2 <- generate_sparse_matrix(8, 6)
  rownames(mat2) <- sample(rownames(mat1))
  mats <- test_add_mat(mats, mat2, 2)

  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(rownames(mats[[2]]), rownames(mat1))

  ## Wrong row names
  mat3 <- mat2
  bad.names <- setdiff(letters, rownames(mats))[1:3]
  rownames(mat3)[1:3] <- bad.names
  expect_error(mats[[2]] <- mat3)

  ## Missed row names
  mat3 <- mat2
  rownames(mat3) <- sample(rownames(mat1))
  mat3 <- mat3[1:4, ]
  expect_error(mats[[2]] <- mat3)

  ## Out-of-bound rows
  mat3 <- generate_sparse_matrix(10, 6)
  rownames(mat3) <- c(rownames(mat1), "AAA", "BBB")
  expect_error(mats[[3]] <- mat3)
})

test_that("rowFlex = 'free', colFlex = 'fixed'", {
  mat1 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(mat1, colFlex = "fixed")

  mat2 <- generate_sparse_matrix(8, 6)
  colnames(mat2) <- sample(colnames(mat1))

  mats <- test_add_mat(mats, mat2, 2)
  expect_equal(colnames(mats), colnames(mat1))
  expect_equal(colnames(mats[[2]]), colnames(mat1))

  colnames(mat2) <- LETTERS[1:ncol(mat2)]
  expect_error(mats[[2]] <- mat2)

  colnames(mat2) <- sample(colnames(mat1))
  mat2 <- mat2[, 1:4]
  expect_error(mats[[2]] <- mat2)
})

test_that("rowFlex = 'fixed', colFlex = 'fixed'", {
  mat1 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(mat1, colFlex = "fixed", rowFlex = "fixed")

  mat2 <- generate_sparse_matrix(8, 6)
  colnames(mat2) <- sample(colnames(mat1))
  rownames(mat2) <- sample(rownames(mat1))

  mats[[2]] <- mat2
  mat.list <- list(mat1, mat2)
  mats <- test_add_mat(mats, mat2, 2)

  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(colnames(mats), colnames(mat1))
  expect_equal(rownames(mats[[2]]), rownames(mat1))
  expect_equal(colnames(mats[[2]]), colnames(mat1))

  rownames(mat2) <- letters[1:nrow(mat2)]
  expect_error(mats[[2]] <- mat2)

  rownames(mat2) <- sample(rownames(mat1))
  mat2 <- mat2[1:4, ]
  expect_error(mats[[2]] <- mat2)

  colnames(mat2) <- LETTERS[1:ncol(mat2)]
  expect_error(mats[[2]] <- mat2)

  colnames(mat2) <- sample(colnames(mat1))
  mat2 <- mat2[, 1:4]
  expect_error(mats[[2]] <- mat2)
})

test_that("rowFlex = 'bounded', colFlex = 'free'", {
  mat1 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(mat1, rowFlex = "bounded")

  mat2 <- generate_sparse_matrix(8, 6)
  rownames(mat2) <- sample(rownames(mat1))

  mats <- test_add_mat(mats, mat2, 2)
  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(rownames(mats[[2]]), rownames(mat1))

  ## within-bound row names
  rownames(mat2) <- sample(rownames(mat1))
  mat2 <- mat2[1:4, ]

  mats <- test_add_mat(mats, mat2, 2)
  expect_equal(rownames(mats), rownames(mat1))
  expect_in(rownames(mats[[2]]), rownames(mat1))

  ## out-of-bound row names
  new.names <- c(setdiff(LETTERS, rownames(mat1)), rownames(mat1))[1:nrow(mat2)]
  rownames(mat2) <- new.names
  expect_error(mats[[2]] <- mat2)
})

test_that("rowFlex = 'free', colFlex = 'bounded'", {
  mat1 <- generate_sparse_matrix(8, 6)

  # colFlex == "bounded"
  mats <- FlexAssays(mat1, colFlex = "bounded")

  mat2 <- generate_sparse_matrix(8, 6)
  colnames(mat2) <- sample(colnames(mat1))

  mats <- test_add_mat(mats, mat2, 2)
  expect_equal(colnames(mats), colnames(mat1))
  expect_equal(colnames(mats[[2]]), colnames(mat1))

  ## within-bound column names
  colnames(mat2) <- sample(colnames(mat1))
  mat2 <- mat2[, 1:4]

  mats <- test_add_mat(mats, mat2, 2)
  expect_equal(colnames(mats), colnames(mat1))
  expect_in(colnames(mats[[2]]), colnames(mat1))

  ## out-of-bound column names
  new.names <- c(setdiff(letters, colnames(mat1)), colnames(mat1))[1:ncol(mat2)]
  colnames(mat2) <- new.names
  expect_error(mats[[2]] <- mat2)
})

test_that("rowFlex = 'bounded', colFlex = 'bounded'", {
  mat1 <- generate_sparse_matrix(8, 6)

  # both "bounded"
  mats <- FlexAssays(mat1, colFlex = "bounded", rowFlex = "bounded")

  mat2 <- generate_sparse_matrix(8, 6)
  colnames(mat2) <- sample(colnames(mat1))
  rownames(mat2) <- sample(rownames(mat1))

  mats <- test_add_mat(mats, mat2, 2)
  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(colnames(mats), colnames(mat1))
  expect_equal(rownames(mats[[2]]), rownames(mat1))
  expect_equal(colnames(mats[[2]]), colnames(mat1))

  ## within-bound row names and column names
  mat2 <- mat2[1:6, 1:4]
  mats <- test_add_mat(mats, mat2, 2)
  expect_in(rownames(mats[[2]]), rownames(mat1))
  expect_in(colnames(mats[[2]]), colnames(mat1))

  ## out-of-bound row names
  new.names <- c(setdiff(LETTERS, rownames(mat1)), rownames(mat1))[1:nrow(mat2)]
  rownames(mat2) <- new.names
  expect_error(mats[[2]] <- mat2)

  ## out-of-bound column names
  new.names <- c(setdiff(letters, colnames(mat1)), colnames(mat1))[1:ncol(mat2)]
  colnames(mat2) <- new.names
  expect_error(mats[[2]] <- mat2)

})

test_that("rowFlex = 'fixed, colFlex = 'bounded'", {
  mat1 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(mat1, rowFlex = "fixed", colFlex = "bounded")

  mat2 <- generate_sparse_matrix(8, 6)
  rownames(mat2) <- sample(rownames(mat1))
  colnames(mat2) <- sample(colnames(mat1))
  mat2 <- mat2[, 1:4]

  mats <- test_add_mat(mats, mat2, 2)
  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(colnames(mats), colnames(mat1))
  expect_equal(rownames(mats[[2]]), rownames(mat1))
  expect_in(colnames(mats[[2]]), colnames(mat1))

  mat3 <- mat2[1:4, ]
  expect_error(mats[[2]] <- mat3)

  mat3 <- mat2
  new.names <- c(setdiff(LETTERS, rownames(mat1)), rownames(mat1))[1:nrow(mat3)]
  rownames(mat3) <- new.names
  expect_error(mats[[2]] <- mat3)

  mat3 <- mat2
  new.names <- c(setdiff(letters, colnames(mat1)), colnames(mat1))[1:ncol(mat3)]
  colnames(mat3) <- new.names
  expect_error(mats[[2]] <- mat3)
})

test_that("rowFlex = 'bounded, colFlex = 'fixed'", {
  mat1 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(mat1, rowFlex = "bounded", colFlex = "fixed")

  mat2 <- generate_sparse_matrix(8, 6)
  rownames(mat2) <- sample(rownames(mat1))
  colnames(mat2) <- sample(colnames(mat1))
  mat2 <- mat2[1:4, ]

  mats <- test_add_mat(mats, mat2, 2)
  expect_equal(rownames(mats), rownames(mat1))
  expect_equal(colnames(mats), colnames(mat1))
  expect_in(rownames(mats[[2]]), rownames(mat1))
  expect_equal(colnames(mats[[2]]), colnames(mat1))

  mat3 <- mat2[, 1:4]
  expect_error(mats[[2]] <- mat3)

  mat3 <- mat2
  new.names <- c(setdiff(LETTERS, rownames(mat1)), rownames(mat1))[1:nrow(mat3)]
  rownames(mat3) <- new.names
  expect_error(mats[[2]] <- mat3)

  mat3 <- mat2
  new.names <- c(setdiff(letters, colnames(mat1)), colnames(mat1))[1:ncol(mat3)]
  colnames(mat3) <- new.names
  expect_error(mats[[2]] <- mat3)
})

test_that("Remove existing assays", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(list(mat1, mat2, mat3))
  mats[[2]] <- NULL

  expect_length(mats, 2)
  rn <- mappedRowNames(rowMap(mats), 2)
  cn <- mappedRowNames(colMap(mats), 2)
  expect_equal(mats[[2]], mat3[rn, cn, drop = FALSE])
})

test_that("Add new assay without dimension names", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(list(mat1, mat2))

  mat3 <- generate_sparse_matrix(8, 6)
  rownames(mat3) <- NULL
  expect_error(mats[[3]] <- mat3)

  mat3 <- generate_sparse_matrix(8, 6)
  colnames(mat3) <- NULL
  expect_error(mats[[3]] <- mat3)

  mat3 <- generate_sparse_matrix(8, 6)
  mat3 <- unname(mat3)
  expect_error(mats[[3]] <- mat3)
})

test_that("Adding new assay with duplicated dimension names should fail", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)

  mats <- FlexAssays(list(mat1, mat2))

  mat3 <- generate_sparse_matrix(8, 6)
  rownames(mat3)[1] <- rownames(mat3)[2]
  expect_error(mats[[3]] <- mat3)

  mat3 <- generate_sparse_matrix(8, 6)
  colnames(mat3)[1] <- colnames(mat3)[2]
  expect_error(mats[[3]] <- mat3)
})

test_that("'assayClasses' works", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- as.matrix(generate_sparse_matrix(8, 6))

  mat.list <- list(mat1, mat2)
  expect_error(mats <- FlexAssays(mat.list, assayClasses = "CsparseMatrix"))

  mats <- FlexAssays(mat.list, assayClasses = c("matrix", "CsparseMatrix"))
  expect_s4_class(mats[[1]], "CsparseMatrix")
  expect_true(is.matrix(mats[[2]]))

  mat3 <- generate_sparse_matrix(8, 6)
  mat3 <- SVT_SparseArray(mat3)
  expect_error(mats[[3]] <- mat3)
})

test_subset <- function(mat.list, rowFlex = "free", colFlex = "free", ...) {
  mats <- FlexAssays(mat.list, rowFlex = rowFlex, colFlex = colFlex)

  mats2 <- mats[...]

  i <- list(...)$i
  j <- list(...)$j

  for (ii in seq_along(mats2)) {
    if (is.null(i)) {
      rn <- mappedRowNames(rowMap(mats), ii)
    } else {
      rn <- mappedRowNames(rowMap(mats)[i, , drop = FALSE], ii)
    }
    if (is.null(j)) {
      cn <- mappedRowNames(colMap(mats), ii)
    } else {
      cn <- mappedRowNames(colMap(mats)[j, , drop = FALSE], ii)
    }
    expect_equal(mats2[[ii]], mat.list[[ii]][rn, cn, drop = FALSE])

    rn2 <- mappedRowNames(rowMap(mats2), ii)
    cn2 <- mappedRowNames(colMap(mats2), ii)
    expect_equal(mats2[[ii]], mat.list[[ii]][rn2, cn2, drop = FALSE])
  }
  mats2
}

test_that("Index FlexAssays", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list)

  # character indices
  new.rows <- union(sample(rownames(mat1), 5), sample(rownames(mats), 10))
  new.cols <- union(sample(colnames(mat1), 5), sample(colnames(mats), 10))
  mats2 <- test_subset(mat.list, i = new.rows)
  mats2 <- test_subset(mat.list, j = new.cols)
  mats2 <- test_subset(mat.list, i = new.rows, j = new.cols)

  # integer indices
  i <- union(sample(1:5, 5), sample(seq_len(nrow(mats)), 10))
  j <- union(sample(1:4, 4), sample(seq_len(ncol(mats)), 10))
  mats2 <- test_subset(mat.list, i = i)
  mats2 <- test_subset(mat.list, j = j)
  mats2 <- test_subset(mat.list, i = i, j = j)

  # logical indices
  i <- sample(c(TRUE, FALSE), nrow(mats), replace = TRUE)
  j <- sample(c(TRUE, FALSE), ncol(mats), replace = TRUE)
  mats2 <- test_subset(mat.list, i = i)
  mats2 <- test_subset(mat.list, j = j)
  mats2 <- test_subset(mat.list, i = i, j = j)
})

test_that("`drop = TRUE` will remove empty matrices and show warnings", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list)

  new.rows <- setdiff(rownames(mats), rownames(mats[[2]]))
  new.cols <- setdiff(colnames(mats), colnames(mats[[3]]))
  expect_warning(mats2 <- mats[new.rows, new.cols])

  # Subset to empty
  mats2 <- mats[integer(), integer()]
  expect_equal(dim(mats2), c(0, 0))
  expect_length(mats2, 0)
})

test_that("`drop = FALSE` will keep empty matrices", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list)

  new.rows <- setdiff(rownames(mats), rownames(mats[[2]]))
  new.cols <- setdiff(colnames(mats), colnames(mats[[3]]))
  mats2 <- test_subset(mat.list, i = new.rows, drop = FALSE)
  mats2 <- test_subset(mat.list, j = new.cols, drop = FALSE)
  mats2 <- test_subset(mat.list, i = new.rows, j = new.cols, drop = FALSE)
})

test_that("`cleanDimMaps` will remove unused maps", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)

  mat.list <- list(mat1, mat2)
  mats <- FlexAssays(mat.list)

  rownames(mat3) <- setdiff(LETTERS, rownames(mats))[1:nrow(mat3)]
  colnames(mat3) <- setdiff(letters, colnames(mats))[1:ncol(mat3)]
  mats[[3]] <- mat3
  mats[[3]] <- NULL

  expect_in(rownames(mat3), rownames(mats))
  expect_in(colnames(mat3), colnames(mats))

  mats <- cleanDimMaps(mats, check = FALSE)

  expect_length(intersect(rownames(mat3), rownames(mats)), 0)
  expect_length(intersect(colnames(mat3), colnames(mats)), 0)
})

test_that("`cleanDimMaps` with 'bounded'", {
  mat1 <- generate_sparse_matrix(8, 6)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(8, 6)
  rownames(mat3) <- setdiff(
    LETTERS,
    union(rownames(mat1), rownames(mat2))
  )[1:nrow(mat3)]
  colnames(mat3) <- setdiff(
    letters,
    union(colnames(mat1), colnames(mat2))
  )[1:ncol(mat3)]

  mat.list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat.list, rowFlex = "bounded", colFlex = "bounded")

  ## Remove mats[[3]], but still can add it back.
  mats[[3]] <- NULL
  expect_in(rownames(mat3), rownames(mats))
  expect_in(colnames(mat3), colnames(mats))
  mats <- test_add_mat(mats, mat3, 3)

  ## Remove mats[[3]] and remove its dimnames, then cannot add it back.
  mats[[3]] <- NULL
  mats <- cleanDimMaps(mats, check = FALSE)

  expect_length(intersect(rownames(mat3), rownames(mats)), 0)
  expect_length(intersect(colnames(mat3), colnames(mats)), 0)

  expect_error(mats[[3]] <- mat3)
})

test_that("FlexAssays 'names'", {
  mat1 <- generate_sparse_matrix(10, 12)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(4, 7)
  mat_list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat_list)

  ## Set new names
  mats2 <- mats
  new_names <- paste0("AAA_", seq_along(mats))
  names(mats2) <- new_names
  expect_equal_no_attr(names(mats2), new_names)
  expect_equal_no_attr(names(assays(mats2)), new_names)
  expect_equal_no_attr(colnames(mats2@rowMap), new_names)
  expect_equal_no_attr(colnames(mats2@colMap), new_names)

  ## Remove all names
  mats2 <- mats
  names(mats2) <- NULL
  expect_null(names(mats2))
  expect_null(names(assays(mats2)))
})

test_that("FlexAssays 'colnames'", {
  mat1 <- generate_sparse_matrix(10, 12)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(4, 7)
  mat_list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat_list)

  expect_equal_no_attr(
    colnames(mats),
    Reduce(union, lapply(mat_list, colnames))
  )

  # Set new colnames
  mats2 <- mats
  new_names <- paste0("AAA_", colnames(mats))
  colnames(mats2) <- new_names
  expect_equal_no_attr(colnames(mats2), new_names)
  expect_equal_no_attr(rownames(mats2@colMap), new_names)
  expect_in(colnames(mats2[[2]]), new_names)

  # Fail to set empty colnames to non-empty FlexAssays
  mats2 <- mats
  expect_warning(colnames(mats2) <- NULL)
  expect_equal_no_attr(colnames(mats2), colnames(mats))

  # Fail to set duplicated colnames to FlexAssays
  expect_error(colnames(mats2)[1] <- colnames(mats2)[2])
  old.names <- colnames(mats2)
  old.names[1] <- old.names[2]
  expect_error(colnames(mats2) <- old.names)
})

test_that("FlexAssays 'rownames'", {
  mat1 <- generate_sparse_matrix(10, 12)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(4, 7)
  mat_list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat_list)

  expect_equal_no_attr(
    rownames(mats),
    Reduce(union, lapply(mat_list, rownames))
  )

  # Set new rownames
  mats2 <- mats
  new_names <- paste0("AAA_", rownames(mats))
  rownames(mats2) <- new_names
  expect_equal_no_attr(rownames(mats2), new_names)
  expect_equal_no_attr(rownames(mats2@rowMap), new_names)
  expect_in(rownames(mats2[[2]]), new_names)

  # Fail to set empty rownames to non-empty FlexAssays
  mats2 <- mats
  expect_warning(rownames(mats2) <- NULL)
  expect_equal_no_attr(rownames(mats2), rownames(mats))

  # Fail to set duplicated rownames to FlexAssays
  expect_error(rownames(mats2)[1] <- rownames(mats2)[2])
  old.names <- rownames(mats2)
  old.names[1] <- old.names[2]
  expect_error(rownames(mats2) <- old.names)
})

test_that("FlexAssays 'dimnames'", {
  mat1 <- generate_sparse_matrix(10, 12)
  mat2 <- generate_sparse_matrix(8, 6)
  mat3 <- generate_sparse_matrix(4, 7)
  mat_list <- list(mat1, mat2, mat3)
  mats <- FlexAssays(mat_list)

  expect_equal_no_attr(dimnames(mats), list(rownames(mats), colnames(mats)))

  # Set new dimnames
  mats2 <- mats
  new_dimnames <- lapply(dimnames(mats), paste0, "_AAA")
  dimnames(mats2) <- new_dimnames
  expect_equal_no_attr(dimnames(mats2), new_dimnames)

  # Ignore empty dimnames
  mats2 <- mats
  new_names <- paste0(rownames(mats), "_AAA")
  expect_warning(dimnames(mats2) <- list(new_names, NULL))
  expect_equal_no_attr(dimnames(mats2), list(new_names, colnames(mats)))

  mats2 <- mats
  new_names <- paste0(colnames(mats), "_AAA")
  expect_warning(dimnames(mats2) <- list(NULL, new_names))
  expect_equal_no_attr(dimnames(mats2), list(rownames(mats), new_names))

  mats2 <- mats
  expect_warning(dimnames(mats2) <- NULL)
  expect_equal_no_attr(dimnames(mats2), dimnames(mats))
  expect_warning(dimnames(mats2) <- list())
  expect_equal_no_attr(dimnames(mats2), dimnames(mats))

  # Raise error when new dimnames with wrong length
  expect_error(dimnames(mats2) <- list(rownames(mats)))
  expect_error(dimnames(mats2) <- list(colnames(mats)))
  expect_error(dimnames(mats2) <- list(rownames(mats), NULL, colnames(mats)))
})

test_that("Coerce FlexAssays to SimpleList", {
  m1 <- generate_sparse_matrix(8, 6)
  m2 <- generate_sparse_matrix(5, 8)
  m3 <- generate_sparse_matrix(3, 4)

  mats1 <- List(m1, m2, m3)
  mats2 <- FlexAssays(mats1)
  expect_equal(as(mats1, "FlexAssays"), mats2)
  expect_equal(as(as(mats2, "SimpleList"), "FlexAssays"), mats2)
  expect_equal(as(as(mats2, "SimpleList"), "SimpleFlexAssays"), mats2)
})

test_that("Coerce FlexAssays to list", {
  m1 <- generate_sparse_matrix(8, 6)
  m2 <- generate_sparse_matrix(5, 8)
  m3 <- generate_sparse_matrix(3, 4)

  mats1 <- list(m1, m2, m3)
  mats2 <- FlexAssays(mats1)
  expect_equal(as(as(as.list(mats2), "SimpleList"), "FlexAssays"), mats2)
})
