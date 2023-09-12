# v0.0.27

* Instead of `include_complement`, the `find_de_neighborhoods` function gains a `add_diff_in_diff` argument. If it is true, the function calculates the difference between the DE results inside the neighborhood vs. outside.
* Change `indices` columns to `neighborhood` and store list of cell name vectors in output of `find_de_neighborhoods`.
* Enforce unique column and row names.

# v0.0.26

* Make the neighborhoods more consistent: (1) include cells which are connected to many cells inside the neighborhood, (2) exclude cells from the neighborhood which are not well connected to the other cells in the neighborhood.
* Add a `control_parameters` argument to `find_de_neighborhoods`.
* Add `BiocNeighbor` as a dependency.

# v0.0.25

* Detect problematic neighborhoods and skip them.
* Replace `test_data_cell_size_factors` by `size_factor_method`, which is more flexibel. Setting `size_factor_method = "ratio"` uses the size factor method described in the original DESeq paper


# v0.0.24

* Fix bug in `find_de_neighborhoods` that meant that accidentally additionally zeros where included in each
neighborhood pseudobulk. The test should have more power now.
* Expose `min_neighborhood_size` argument in `find_de_neighborhoods`.
* Add `test_data_cell_size_factors` argument to `find_de_neighborhoods` which is useful if the function is called
with a subsetted `fit` argument.

# v0.0.23

* Improve alignment functions: simplify algorithm, find linear approximation to Harmony's steps,
include an intercept.
* Avoid calling private methods from `harmony`.
* Convert character columns in `colData` to factors to avoid problems when dividing data into
test and training data.
* Fix bug in `find_de_neighborhoods` where I didn't embrace an argument.
* Remove `BiocNeighbors` dependency.

# v0.0.21

* Minor bug fix in `find_de_neighborhoods`. The function threw an error if `alignment_design != design`. 
* Better error messages if `find_de_neighborhoods` is called without having called `test_de` before.

# v0.0.20

* Change defaults for `find_de_neighborhoods`. Increase the `ridge_penalty` and add a `min_neighborhood_size = 10` argument
to avoid creation of very small neighborhoods.

# v0.0.19

* Add new `test_fraction` argument to `lemur()` function. It automatically defines a hold-out datasets for the fitting step.
These hold-out data is used to infer the differential expression of the neighborhoods in `find_de_neighborhoods`. This change
addresses the double-dipping problem, where it was previously left to the user to provide an independent matrix for the 
`find_de_neighborhoods` function.
* As a consequence of these changes, the structure of `lemur_fit` objects has changed. They gain three new fields called
`fit$test_data`, `fit$training_data`, and `fit$is_test_data`.
* The order and names of the arguments for `find_de_neighborhoods` has changed.

# v0.0.18

* Remove `alignment_method` field from `lemur_fit` objects as it was not used for anything.

# v0.0.17

* Rename argument name for `align_by_template` from `alignment_template` to `template`
* Tweak algorithm for alignment to take cluster sizes into account during optimization

# v0.0.13-v0.0.16

* Change in the alignment model. Previously, the method tried to align cells using
rotations and / or stretching, however, the method could not represent reflections! 
To fix this, I now allow arbitrary linear transformations where $R(x) = (I + sum_k x_k V_k)^{-1}$. The
new alignment is more flexible and easier to infer. The downside is the term inside the parantheses can be 
singular which would lead to an error.
* Skip iteration step: first infer centering and then infer latent space. Previously, I iterated between these steps 
but that either didn't add anything or actually degraded the results.
* Set `center = FALSE` in `find_base_point`. Centering the data before fitting the base point caused
problems and made the data look less integrated in some cases.
* Remove ambient PCA step. This was originally conceived as an performance optimization, however
it had detrimental effects on the inference. Since a few version it was skipped per default, so removing
it should not change the inference.
* Add `linear_coefficient_estimator` to give more flexibility how or if the conditions are centered.
* Reduce the `minimum_cluster_membership` default to `0.001` in `align_harmony` to make it more sensitive.
* Make `test_global` an internal function again until further more extensive testing.
* Remove `base_point` argument from `lemur()`. It wasn't used anyways.

# pre v0.0.13
* Refactor `find_de_neighborhoods`: the function can now combine the results of
different directions, selection criteria, and pseudobulk test (on counts or 
continuous values). To implement this, I changed the names of the arguments and
added parameters.
* Remove many superfluous method generics and only provide the accession via `$`
* Fix documentation warnings
* Rename class from 'lemur_fit_obj' to 'lemur_fit'
* Store 'contrast' in `lemur_fit` after calling `test_de`
* Add option to fit count model in `find_de_neighborhoods` with [edgeR](https://bioconductor.org/packages/edgeR/)
