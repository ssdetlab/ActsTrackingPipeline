#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>
#include <memory>

class CompositeMagField : public Acts::MagneticFieldProvider {
 public:
  struct FieldComponent {
    Acts::Extent extent;
    std::shared_ptr<Acts::MagneticFieldProvider> fieldProvider;
  };
  using FieldComponents = std::unordered_map<std::size_t, FieldComponent>;

  using FieldComponentExtents = std::unordered_map<std::size_t, Acts::Extent>;
  using FieldComponentCaches =
      std::unordered_map<std::size_t, Acts::MagneticFieldProvider::Cache>;

  struct Cache {
    FieldComponentExtents m_componentExtents;
    FieldComponentCaches m_componentCaches;

    explicit Cache(const FieldComponentExtents& componentExtents)
        : m_componentExtents(componentExtents) {}

    Cache(std::size_t id, const FieldComponents& components,
          const Acts::MagneticFieldContext& mctx) {
      for (const auto& [cId, component] : components) {
        m_componentExtents.insert({cId, component.extent});
        m_componentCaches.insert(
            {cId, component.fieldProvider->makeCache(mctx)});
      }
    }
  };

  CompositeMagField(std::size_t id, const FieldComponents& fieldComponents);

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
  std::size_t m_id;
  FieldComponents m_fieldComponents;
};
