test_that("Count table of mutation types", {

  mut_types <- c("C[T>G]G", "T[C>T]C", "A[T>A]G", "G[C>T]A", "C[C>T]T", "T[C>T]C", "C[T>G]G")

  expected <- data.frame(mut_types = factor(c("A[T>A]G", "C[C>T]T", "C[T>G]G", "G[C>T]A", "T[C>T]C")),
                         Freq = c(1,1,2,1,2))

  expect_equal(count_table(mut_types), expected)
})

