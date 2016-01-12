
################################################################################
# Utilities when processing data
################################################################################

#' @title Expand annotation using phyloseq tables
#' @export
supp_annotation <- function(mvar_object, physeq, n1 = "species", n2 = "site") {
  mvar_object@table[[n1]]@annotation <- data.frame(mvar_object@table[[n1]]@annotation,
                                                   physeq@tax_table)
  mvar_object@table[[n2]]@annotation <- data.frame(mvar_object@table[[n2]]@annotation,
                                                   physeq@sam_data)
  mvar_object
}
