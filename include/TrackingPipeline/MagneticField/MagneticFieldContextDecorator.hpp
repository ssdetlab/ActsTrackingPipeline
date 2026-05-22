#pragma once

#include "Acts/MagneticField/MagneticFieldContext.hpp"

#include <memory>
#include <string>

#include "TrackingPipeline/Infrastructure/IContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldStoreCollection.hpp"

class MagneticFieldContextDecorator : public IContextDecorator {
 public:
  explicit MagneticFieldContextDecorator(
      const std::shared_ptr<MagneticFieldStoreCollection>&
          fieldStoreCollectionPtr)
      : m_fieldStoreCollectionPtr(fieldStoreCollectionPtr) {}

  ProcessCode decorate(AlgorithmContext& ctx) override {
    ctx.magFieldContext = Acts::MagneticFieldContext{
        m_fieldStoreCollectionPtr->magneticFieldStorePtr(ctx.eventNumber)};
    return ProcessCode::SUCCESS;
  }

  const std::string& name() const override { return m_name; };

 private:
  std::shared_ptr<MagneticFieldStoreCollection> m_fieldStoreCollectionPtr;

  std::string m_name = "MagneticFieldContextDecorator";
};
