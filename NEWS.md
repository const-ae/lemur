# lemur devel

* Refactor `find_de_neighborhoods`: the function can now combine the results of
different directions, selection criteria, and pseudobulk test (on counts or 
continuous values). To implement this, I changed the names of the arguments and
added parameters.
* Remove many superfluous method generics and only provide the accession via `$`
* Fix documentation warnings
* Rename class from 'lemur_fit_obj' to 'lemur_fit'
* Store 'contrast' in `lemur_fit` after calling `test_de`
* Add option to fit count model in `find_de_neighborhoods` with [edgeR](https://bioconductor.org/packages/edgeR/)
