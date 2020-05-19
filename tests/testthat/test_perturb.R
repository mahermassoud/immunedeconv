# Exampled melted gene signature matrix
eg.sig <- read_csv("eg_sig_melt.csv")

test_that("truncate signature works", {
  res = truncate_signature(eg.sig, 3)
 
  n_cell_type = eg.sig$cell_type %>% unique() %>% length()
  res = truncate_signature(eg.sig, 3, return_melted=FALSE)
  expect_equal(ncol(res$perturbed), n_cell_type, info="n_cell_type consistent")
})