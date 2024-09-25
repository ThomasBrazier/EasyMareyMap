test_that("The marey_map object is conform", {
  m = data.frame(set = "Test",
                 map = "1",
                 mkr = c(1, 2, 3, 4),
                 phys = c(11, 22, 33, 44),
                 gen = c(11, 22, 33, 44),
                 vld = TRUE)
  df = marey_map(m, chromosome = '1')
  expect_equal(class(df), 'marey_map')
})
