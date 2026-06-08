#' Prepare multi-target chord data
#'
#' @description
#' Creates a link-level data.table where a single sender population is mapped
#' to all potential receiver populations for every treatment.
#'
#' @param circData_dt Cleaned data.table from the input CSV.
#' @param senderPop_v The population to act as the sender.
#' @param receiverPop_v The populations to act as receiver. NULL (default) uses all.
#' @param senderDataCol_v column that has values for sender
#' @param receiverDataCol_v column that has values for receiver
#' @param popCol_v character vector indicating which column of circData_dt has sector info
#' @param arcColorCol_v column indicating what to color arcs by (usually treatment)
#'
#' @return A data.table with [sender, receiver, treatment, senderWidth, receiverWidth]
#' @export
prepareMultiTargetData <- function(circData_dt, senderPop_v, receiverPop_v = NULL, 
                                   senderDataCol_v, receiverDataCol_v,
                                   popCol_v = "Pop", arcColorCol_v = "Treatment") {
  
  ### Subset full data for just provided sender and associated info
  sender_dt <- circData_dt[get(popCol_v) == senderPop_v, mget(c(popCol_v, senderDataCol_v, arcColorCol_v))]
  
  ### Prepare the receiver pool (all populations)
  receiver_dt <- circData_dt[, mget(c(popCol_v, receiverDataCol_v, arcColorCol_v))]
  colnames(receiver_dt)[1] <- "receiverPop"
  
  ### Subset Receiver populations
  if (!is.null(receiverPop_v)) {
    receiver_dt <- receiver_dt[receiverPop %in% receiverPop_v,]
  }
  
  ### Join on Treatment to create the cross-product
  ### This maps the single sender row to all population rows for that treatment
  links_dt <- merge(x = sender_dt, y = receiver_dt, by = "Treatment")
  
  ### Adjust widths
  ### We divide the sender's data by the number of receivers so the 
  ### total sector width reflects the actual value.
  links_dt[, numReceivers := .N, by = get(arcColorCol_v)]
  links_dt[, senderWidthPart := get(senderDataCol_v) / numReceivers]
  
  ### Receiver side doesn't need to be adjusted
  links_dt[, receiverWidthPart := get(receiverDataCol_v)]
  
  links_dt <- setorderv(links_dt, c(arcColorCol_v, "receiverPop"))
  
  return(links_dt)
  
} ### End prepareMultiTargetData