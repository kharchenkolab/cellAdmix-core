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

}  // namespace celladmix
