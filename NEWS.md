# seqgendiff 0.3.0

This has been a massive rewrite of the seqgendiff package.

- `poisthin()` is now defunct. The two-group model is now implemented in
  the `thin_2group()` function. I'll keep it around since some of my old
  simulation code depends on it.
- The main functions for thinning are now `thin_diff()`, `thin_2group()`, 
  `thin_lib()`, and `thin_gene()`.
- Note that these functions, unlike `poisthin()`, do not have functionality to
  subset count matrices. This is on purpose. I wanted the functionality
  of these thinning functions to be simpler.
- Unlike `poisthin()`, which can only handle the two-group model, `thin_diff()`
  can handle generically any design, while still controlling the level of 
  correlation between the design variables and the surrogate variables.
- Functions now exist to convert the simulation output to common 
  Bioconductor S4 classes. See `ThinDataToSummarizedExperiment()` 
  and `ThinDataToDESeqDataSet()`.

# seqgendiff 0.2.0

* `corassign()` lets you make group assignment that is correlated with hidden factors.
* In `poisthin()`, the `group_assign = "cor"` option uses `corassign()` to make group assignments.
