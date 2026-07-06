#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>
#include <memory>

#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

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

    Cache(std::size_t id, const FieldComponents& defaultComponents,
          const Acts::MagneticFieldContext& mctx) {
      if (!mctx.hasValue()) {
        for (const auto& [cId, defaultComponent] : defaultComponents) {
          m_componentExtents.insert({cId, defaultComponent.extent});
          m_componentCaches.insert(
              {cId, defaultComponent.fieldProvider->makeCache(mctx)});
        }
      } else {
        const auto& store =
            mctx.get<std::shared_ptr<MagneticFieldStore>&>()->store;
        if (store.contains(id)) {
          auto componentExtents = store.at(id).as<Cache>().m_componentExtents;
          auto componentCaches = store.at(id).as<Cache>().m_componentCaches;

          for (const auto& [cId, defaultComponent] : defaultComponents) {
            if (componentExtents.contains(cId)) {
              m_componentExtents.insert({cId, componentExtents.at(cId)});
              m_componentCaches.insert(
                  {cId, defaultComponent.fieldProvider->makeCache(mctx)});
            } else {
              m_componentExtents.insert({cId, defaultComponent.extent});
              m_componentCaches.insert(
                  {cId, defaultComponent.fieldProvider->makeCache(mctx)});
            }
          }
        } else {
          for (const auto& [cId, defaultComponent] : defaultComponents) {
            m_componentExtents.insert({cId, defaultComponent.extent});
            m_componentCaches.insert(
                {cId, defaultComponent.fieldProvider->makeCache(mctx)});
          }
        }
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
