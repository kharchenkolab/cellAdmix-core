// Store-backed fit pipeline that processes molecule rows by cell blocks.

#pragma once

#include <optional>
#include <string>

#include "celladmix/input_store.hpp"
#include "celladmix/pipeline.hpp"
#include "celladmix/run_store.hpp"

namespace celladmix {

struct StorePipelineResult {
  RunManifest manifest;
  SparseNmfResult nmf;
  BasicPipelineTiming timing;
};

StorePipelineResult run_basic_pipeline_store(
    const std::string& store_dir,
    const RunSourceInfo& source,
    const BasicPipelineOptions& options,
    const RunStorageOptions& storage_options,
    const std::string& out_dir,
    const std::optional<std::string>& analysis_crop = std::nullopt,
    bool report_ncv_umap = false,
    const CellCountMatrix* preloaded_counts = nullptr,
    bool verbose = false);

// Loss of a persisted run's sampled training rows under given gene profiles
// (row-major, rank x the run's genes), in the run's own NMF objective. Row
// usages are obtained by fixed-profile projection, so the value compares
// profile sets on the same rows: the run's own profiles against a consensus
// or a profile set from elsewhere. Profiles are expected in raw gene space,
// as stored for weighted least-squares fits; for KL variants with gene
// weights they are carried into the run's transformed feature space.
double training_loss_under_profiles(
    const std::string& run_path_or_dir,
    const std::vector<double>& profiles,
    int rank,
    int num_threads = 1);

}  // namespace celladmix
