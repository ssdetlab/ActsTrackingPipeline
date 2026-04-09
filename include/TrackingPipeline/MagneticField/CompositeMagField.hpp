#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>
#include <memory>

class CompositeMagField : public Acts::MagneticFieldProvider {
 public:
  struct FieldComponent {
    std::size_t id;
    Acts::Extent extent;
    std::shared_ptr<Acts::MagneticFieldProvider> fieldProvider;
  };
  using FieldComponents = std::vector<FieldComponent>;

  struct Cache {
    std::unordered_map<std::size_t, Acts::MagneticFieldProvider::Cache>
        componentCaches;
    Cache(const FieldComponents& components,
          const Acts::MagneticFieldContext& mctx) {
      for (const auto& component : components) {
        // std::cout << "COMPONENT ID " << component.id << "\n";
        componentCaches.insert(
            {component.id, component.fieldProvider->makeCache(mctx)});
      }
    }
  };

  CompositeMagField(const FieldComponents& fieldComponents);

  ~CompositeMagField() override;

  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override;

  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& position,
      MagneticFieldProvider::Cache& cache) const override;

  Acts::Result<Acts::Vector3> getFieldGradient(
      const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
      MagneticFieldProvider::Cache& cache) const override;

 private:
  FieldComponents m_fieldComponents;
};
