## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  dpi = 96
)

## ----setup--------------------------------------------------------------------
library(phylorug)

## ----load-sample--------------------------------------------------------------
names(sample_trees)

backbone <- sample_trees[["70p_uce"]]
others   <- sample_trees[names(sample_trees) != "70p_uce"]

## ----check-taxa---------------------------------------------------------------
check_taxa(backbone, others)

## ----Support_type_vector------------------------------------------------------
support_type <- c(
  "70p_ASTRAL_partition_entropy" = "lpp",
  "70p_ASTRAL_uce"               = "lpp",
  "70p_ghost"                    = "ufboot",
  "70p_partition_entropy"        = "ufboot"
)

## ----build-npm----------------------------------------------------------------
npm <- node_presence_matrix(backbone, others, support_col = c(1, 2), support_type = support_type)

## ----quick-presence,  fig.height=5.8, fig.width= 6.5, fig.dpi=300-------------
plot_phylorug(backbone, npm, show_support = FALSE)

## ----presence-file, eval = FALSE----------------------------------------------
# plot_phylorug(backbone, npm, file = "presence.pdf", width = 12, height = 8)

## ----support-default, fig.height=5.8, fig.width= 6.5, fig.dpi=300-------------
plot_phylorug(
  backbone, npm,
  mode         = "support")

## ----support-refined, fig.height=5.8, fig.width= 6.5, fig.dpi=300-------------
# For publication output, add:
# file = "my_support_plot.pdf"
plot_phylorug(
  backbone, npm,
  width             = 6.5,
  height            = 5.8,
  mode              = "support",
  rug_position      = "inside",
  cell_scale        = 0.38,
  cex               = 0.90,
  font              = 3,
  dot_cex           = 1.3,
  show_support      = TRUE,
  support_label_col = "red",
  support_label_cex = 0.4
)

## ----read-culico--------------------------------------------------------------
tree_dir <- system.file("extdata", "culicomorpha", package = "phylorug")
trees    <- read_trees(tree_dir)
names(trees)

## ----root-culico--------------------------------------------------------------
trees <- lapply(trees, function(tr) {
  tr <- ape::root(tr, outgroup = c("Coboldia_fuscipes",
                                   "Phlebotomus_chinensis",
                                   "Clogmia_albipunctata"),
                  resolve.root = TRUE)
  ape::drop.tip(tr, c("Coboldia_fuscipes",
                      "Phlebotomus_chinensis",
                      "Clogmia_albipunctata"))
})

## ----backbone-culico----------------------------------------------------------
backbone <- trees[["Matrix1-kpi_PMSF(H1.guide)"]]
others   <- trees[names(trees) != "Matrix1-kpi_PMSF(H1.guide)"]


## ----check-taxon-consistency-culico-------------------------------------------
check_taxa(backbone, others)

## ----support-type-culico------------------------------------------------------
support_type <- c(
  "Matrix1-kpi_ASTRAL"            = "lpp",
  "Matrix1-smart_ASTRAL"          = "lpp",
  "Matrix1-kpi_partitioning"      = "sh_alrt",
  "Matrix1-kpi_PMSF(H2.guide)"    = "sh_alrt",
  "Matrix1_LG_C20_F_R"            = "sh_alrt",
  "Matrix1-smart_partitioning"    = "sh_alrt",
  "Matrix1-smart_PMSF(H1.guide)"  = "sh_alrt",
  "Matrix1-smart_PMSF(H2.guide)"  = "sh_alrt",
  "Matrix2_partitioning"          = "sh_alrt"
)

## ----npm-culico---------------------------------------------------------------
npm <- node_presence_matrix(backbone, others, support_type = support_type)

## ----plot-culico, fig.width = 8, fig.height = 10, fig.dpi=300-----------------
plot_phylorug(backbone, npm, show_support = FALSE)

## ----plot-culico-support,fig.width = 8, fig.height = 10, fig.dpi=300----------
plot_phylorug(backbone, npm,
              mode               = "support",
              show_support       = TRUE,
              show_support_idx   = 1,
              support_label_cex  = 0.35,
              support_label_col  = "red",
              cell_scale         = 0.3,
              rug_position       = "outside")

## ----read trees-beetles-70p---------------------------------------------------
tree_dir <- system.file("extdata", "beetles_70p", package = "phylorug")
trees    <- read_trees(tree_dir)
names(trees)

## ----root-and-remove-outgroup-beetles-70p-------------------------------------
trees <- lapply(trees, function(tr) {
  tr <- ape::root(tr, outgroup = c("NicorbUCE", "NicvesUCE"),
                  resolve.root = TRUE)
  ape::drop.tip(tr, c("NicorbUCE", "NicvesUCE"))
})

## ----load the lookup table----------------------------------------------------
biogeo_path <- system.file("extdata", "beetles_70p", "biogeo.csv",
                           package = "phylorug")
biogeo <- read.csv(biogeo_path)
head(biogeo)

## ----translate-tips-beetles-70p-----------------------------------------------
trees <- translate_tips(trees, biogeo,
                        from_col = "specimen_code",
                        to_col   = "species_name")

## ----backbone-beetles-70p-----------------------------------------------------
backbone <- trees[["70p_uce"]]
others   <- trees[names(trees) != "70p_uce"]
check_taxa(backbone, others)

## ----support-type-beetles-70p-------------------------------------------------
support_type <- c(
  "70p_partition_entropy"        = "ufboot",
  "70p_ghost"                    = "ufboot",
  "70p_ASTRAL_uce"               = "lpp",
  "70p_ASTRAL_partition_entropy" = "lpp"
)

## ----npm-beetles-70p----------------------------------------------------------
npm <- node_presence_matrix(backbone, others, support_col = c(1, 2),
                            support_type = support_type)

## ----plot phylorug , fig.height= 12, fig.width= 8-----------------------------
plot_phylorug(backbone, npm)

## ----plot-support-beetles-70p, fig.height = 12, fig.width = 10, fig.dpi= 300----
# For publication output, add:
# file = "beetles_support.pdf"
# height = 12 
# width = 10
plot_phylorug(
  backbone,
  npm,
  mode               = "support",
  show_support       = TRUE,
  show_support_idx   = c("1" = "sh_alrt", "2" = "ufboot"),
  support_label_col  = "red",
  support_label_cex  = 0.38,
  cell_scale         = 0.35,
  cex                = 0.8
)

## ----reload-sample------------------------------------------------------------
backbone <- sample_trees[["70p_uce"]]
others   <- sample_trees[names(sample_trees) != "70p_uce"]
support_type <- c(
  "70p_ASTRAL_partition_entropy" = "lpp",
  "70p_ASTRAL_uce"               = "lpp",
  "70p_ghost"                    = "ufboot",
  "70p_partition_entropy"        = "ufboot"
)
npm <- node_presence_matrix(backbone, others, support_col = c(1, 2),
                            support_type = support_type)

## ----dot-on-very-high---------------------------------------------------------
plot_phylorug(backbone, npm, mode = "support", dot_on_very_high = FALSE)

## ----universal_threshold------------------------------------------------------
npm_universal <- node_presence_matrix(backbone, others, support_col = c(1, 2))

plot_phylorug(
  backbone, npm_universal,
  mode = "support"
)

## ----grid-demo----------------------------------------------------------------
# rug with a flattened grid
plot_phylorug(backbone, npm, n_rows = 1, n_cols = 4)

## ----tree-appearance----------------------------------------------------------
plot_phylorug(backbone, npm, cex = 0.6, font = 3, edge.width = 1.5)

## ----include-bb---------------------------------------------------------------
plot_phylorug(backbone, npm, include_backbone = TRUE)

## ----include-bb-support-------------------------------------------------------
plot_phylorug(backbone, npm,
              mode             = "support",
              include_backbone = "ufboot")


## ----custom-thresh, fig.dpi=300-----------------------------------------------
my_thresholds <- list(
  ufboot = c(very_high = 98, high = 90, moderate = 70),
  lpp    = c(very_high = 0.99, high = 0.90, moderate = 0.70)
)

plot_phylorug(
  backbone, npm,
  mode       = "support",
  thresholds = my_thresholds
)

## ----nodes-by-taxa------------------------------------------------------------
plot_phylorug(
  backbone, npm,
  mode = "support",
  nodes = list(
    c("Gyronotus_pumilus_STL10003", "Scarabaeus_westwoodi_STL10034"),
    c("Nanos_dubitatus_STL5001", "Catharsius_sp._STL10033")
  )
)

## ----nodes-by-id--------------------------------------------------------------
plot_phylorug(backbone, npm, nodes = 18)

## ----find-disagreeing-nodes---------------------------------------------------
disagree_ids <- as.integer(rownames(npm$presence))[
  apply(npm$presence, 1, function(row) any(row == 0))
]
disagree_ids

lapply(disagree_ids, function(id) ape::extract.clade(backbone, id)$tip.label)

