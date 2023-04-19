# v0.0.15

* Skip iteration step: first infer centering and then infer latent space. Previously, I iterated between these steps 
but that either didn't add anything or actually degraded the results.

# v0.0.13-v0.0.14 (19 of April 2023)

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
