`cellAdmix-core` vendors the following third-party header-only dependencies under [include](include):

- `third_party/nanoflann.hpp`
- `third_party/hnswlib/`
- `umappp/`
- `knncolle/`
- `kmeans/`
- `subpar/`
- `aarand/`
- `irlba/`

These headers are used for exact and approximate spatial search, kNN graph
construction, k-means clustering, PCA-related linear algebra helpers, and
UMAP embedding. The vendored copies retain their upstream license notices in
the distributed source files.
