#' Functions for perturbing cell-type x gene signatures for simulations
#'
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import tibble
#' 
#' @name perturb
NULL

#' @param signature data.frame where first column is gene, 
#'   second column is cell type, third columns value for that gene,type pair
#'   third column is optional
#' @param n_gene_rm number of genes per cell type to randomly remove
#' @param return_melted boolean whether to return the signature matrix in
#'   melted format or not. If `FALSE`, returns n_gene x n_cell_type data.frame
#'
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import tibble
#' 
#' @name truncate_signature
#' @export
truncate_signature = function(signature, n_gene_rm, return_melted=TRUE) {
  cl.sig <- as_tibble(signature)
  names(cl.sig)[1] <- "gene"
  names(cl.sig)[2] <- "cell_type"
  if(ncol(cl.sig) == 3) {
    names(cl.sig)[3] <- "value"
    cl.sig <- cl.sig %>% filter(value > 0)
  } 

  summ.sig <- cl.sig %>% 
    group_by(cell_type) %>%
    summarise(genes=list(gene), n_gene=n())

  summ.sig <- summ.sig %>%
    mutate(
      n_rm=ifelse(n_gene_rm > as.integer(0.5*n_gene), as.integer(0.5*n_gene), n_gene_rm),
      n_keep=n_gene-n_rm,
      keep_gene=map2(genes, n_keep, function(gs, n) sample(gs, n)),
      rm_gene=map2(genes, keep_gene, setdiff)
    )

  perturb.sig <- summ.sig %>%
    dplyr::select(keep_gene, cell_type) %>%
    unnest(keep_gene) %>%
    rename(gene=1,cell_type=2)
 
  # Merge in values if we need them
  if(ncol(cl.sig) == 3) {
    perturb.sig <- left_join(perturb.sig, cl.sig, by=c("gene","cell_type")) 
  }
  
  if(!return_melted) {
    # Add value column if there isn't one
    if(ncol(cl.sig) < 3) {
      perturb.sig$value <- 1
    }
    perturb.sig <- spread(perturb.sig, key=cell_type, value=value, fill=0) %>%
      as.data.frame() %>%
      column_to_rownames("gene")
  }
   
  out = list(perturbed=perturb.sig, summary=summ.sig)
  return(out)
}