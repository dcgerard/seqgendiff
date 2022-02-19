# seqgendiff 1.2.3

- Removes `{optmatch}` as a suggested package since it is no longer on CRAN.
- Removes `LazyData: true` from DESCRIPTION since there is no 'data' directory.

# seqgendiff 1.2.2

- Updates all citations to BMC Bioinformatics published manuscript.

# seqgendiff 1.1.1

- Updates to the documentation. I included citations throughout the documentation to direct users to where they can get details on the method.
- The citation file now indicates the preprint on bioRxiv.

# seqgendiff 1.1.0

Fixes a lot of things for CRAN resubmission.

- Adds more information on return values.
- Changes title to less than 65 characters.
- Uses local environments rather than global environment to assess
  messages.
- Removes `seqgendiff::EigenDiff()`. Replaces its usage 
  with `cate::est.factor.num()`. This is fine since it was only 
  used in the now defunct `seqgendiff::poisthin()`.

# seqgendiff 1.0.0

- The biggest change here is that the `{optmatch}` package is now only
  suggested rather than imported. This is because the `{optmatch}` package
  is under a super weird license that I didn't previously know about.
- The user may also now specify the permutation method in the thinner
  functions.
- The Hungarian algorithm, implemented in the `{clue}` package, seems to
  work just as well as `{optmatch}`, and so I added it as an
  option. However, since I used `{optmatch}` in the simulations for the
  paper, I have kept `permute_method = "optmatch"` as the default
  option.

# seqgendiff 0.4.1

- I added `select_counts()`, a function that will subsample the rows (genes)
  and columns (samples) of a RNA-seq count matrix. It is generally
  recommended that you do this subsampling each iteration of a simulation
  study so that your results do not depend on the specific structure of
  your data. The samples are just selected randomly. There are four different
  criteria for selecting the genes.
- I also added citation information. This will be updated after I submit
  the manuscript to bioRxiv.

# seqgendiff 0.4.0

- This version mostly updates the documentation. 
- Beyond this major documentation update, we also added `thin_all()`, 
  a function that uniformly thins all counts.

# seqgendiff 0.3.0

This has been a massive rewrite of the `{seqgendiff}` package.

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

- `corassign()` lets you make group assignment that is correlated with hidden factors.
- In `poisthin()`, the `group_assign = "cor"` option uses `corassign()` to make group assignments.
