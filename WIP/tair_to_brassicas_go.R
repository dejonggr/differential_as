
tairgo = data.frame(genes = unlist(xx), go_ids = names(unlist(xx))) %>%
  separate(go_ids, c("go_ids", "origin"), sep = "\\.")

