#pragma once

#include "ActsAlignment/Kernel/Alignment.hpp"

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"

namespace detail {

ActsAlignment::AlignedTransformUpdater makeGlobalAlignmentUpdater(
    AlignmentContext& alignCtx);

ActsAlignment::AlignedTransformUpdater makeLocalAlignmentUpdater(
    AlignmentContext& alignCtx);

}  // namespace detail
