# Exampled melted gene signature matrix
eg.sig <- read_csv("eg_sig_melt.csv")
eg.nm.sig <- read_csv("eg_sig.csv") %>% as.data.frame() %>% column_to_rownames("gene")

test_that("truncate signature works", {
  n_cell_type = eg.sig$cell_type %>% unique() %>% length()
  res = truncate_signature(eg.sig, 3, return_melted=FALSE)
  expect_equal(ncol(res$perturbed), n_cell_type, info="n_cell_type consistent")
})

test_that("perturb melted input", {
  res = perturb_sig(eg.sig, 0.5)
  n_cell_type = eg.sig$cell_type %>% unique() %>% length()
  expect_equal(ncol(res$perturbed), n_cell_type, info="n_cell_type consistent")
  
  res = perturb_sig(eg.nm.sig, 0.5, input_melted=FALSE)
  expect_equal(ncol(res$perturbed), n_cell_type, info="n_cell_type consistent")
})