#pragma once

#include "ActsAlignment/Kernel/Alignment.hpp"

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"

namespace detail {

ActsAlignment::AlignmentTransformUpdater makeGlobalAlignmentUpdater(
    AlignmentContext& alignCtx);

ActsAlignment::AlignmentTransformUpdater makeLocalAlignmentUpdater(
    AlignmentContext& alignCtx);

}  // namespace detail
