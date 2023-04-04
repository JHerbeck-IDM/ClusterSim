tree_tbl <- read.csv("line_list.csv")


make_phylo = function(tree_tbl, fix_dd_diff = T, fix_nontips = T){
  
  # Note, this parent child mapping assumes that all branches will
  # terminate in a tip, and tip is pre-specifed. An alternate reconstruction
  # of a tip would be to look at names with no descendants, which would be safer
  tmp = tree_tbl %>%
    mutate(tip_fixed = !(name %in% p_name))
  
  if(fix_nontips){
    maxit = 1000
    it = 1
    # iteratively remove internal nodes that don't result in a tip in the data set
    while(sum(with(tmp, tip!=tip_fixed)) > 0 && it < maxit){
      tmp = tmp %>% filter(tip==tip_fixed) %>% mutate(tip_fixed = !(name %in% p_name))
      it = it + 1
    }
  }
  
  tmp = tmp %>%
    arrange(desc(tip), dd, name) %>%
    mutate(idx = as.integer(factor(name, levels = unique(.$name)))) %>%
    mutate(parent = idx[match(p_name, name)],
           child = idx) %>%
    relocate(parent, child) #relocates "parent" and "child" columns to the front
  
  if(fix_dd_diff) tmp = tmp %>% mutate(dd_diff = pmax(dd_diff, 0))
  
  phy = list(edge = tmp %>% filter(!is.na(p_name)) %>% with(cbind(parent, child)),
             edge.length = tmp %>% filter(!is.na(p_name)) %>% pull(dd_diff),
             tip.label = tmp %>% filter(tip_fixed) %>% arrange(idx) %>% pull(name),
             Nnode = nrow(tmp) - sum(tmp$tip_fixed),
             node.label = tmp %>% filter(!tip_fixed) %>% arrange(idx) %>% pull(name))
  class(phy) = "phylo"
  #attr(phy, "order")= "cladewise" # is this necessary? TBD -- ggtree fails without ladderize = F
  attr(phy, "min_date") = min(tmp$dd)
  phy
}

make_treedata = function(tree_tbl){
  phy = make_phylo(tree_tbl)
  as_tibble(phy) %>%
    inner_join(tree_tbl %>%
                 mutate(type = if_else(!tip, as.character(NA), surv_type)) %>%
                 select(label = name, dd, dd_diff, region, type, district),
               by = "label") %>%
    as.treedata()
}


