# CTUtools

Common R functions and reference data used by DKMS CTU.

See `?CTUtools` for list of functions and reference datasets.

Maintained by Bose Falk, <falk@dkms.de>

All functions are documented, use ?functionname to see details of how to use them and examples.

All functions are written to do the classification / operation on one datapoint instead of a whole dataframe at once (i.e. `KIR_first_field(x)` takes as x a single KIR string), so they'll need to be wrapped in `apply()` or `dplyr::rowwise() %>% dplyr::mutate()` inside your scripts.
