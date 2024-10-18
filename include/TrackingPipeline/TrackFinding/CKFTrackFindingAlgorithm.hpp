#pragma once

#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"

#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/EventData/TrackContainer.hpp"

/// The map(-like) container accessor
template <typename container_t>
struct SimpleContainerAccessor {
    using Container = container_t;
    using Key = typename container_t::key_type;
    using Value = typename container_t::mapped_type;
    
    /// This iterator adapter is needed to have the deref operator return a single
    /// source link instead of the map pair <GeometryIdentifier,SourceLink>
    struct Iterator {
        using BaseIterator = typename container_t::const_iterator;
    
        using iterator_category = typename BaseIterator::iterator_category;
        using value_type = typename BaseIterator::value_type;
        using difference_type = typename BaseIterator::difference_type;
        using pointer = typename BaseIterator::pointer;
        using reference = typename BaseIterator::reference;
    
        Iterator& operator++() {
        ++m_iterator;
            return *this;
        }
    
        bool operator==(const Iterator& other) const {
            return m_iterator == other.m_iterator;
        }
    
        bool operator!=(const Iterator& other) const { return !(*this == other); }
    
        Acts::SourceLink operator*() const {
            const auto& sl = m_iterator->second;
            return Acts::SourceLink{sl};
        }
    
        BaseIterator m_iterator;
    };
    
    // pointer to the container
    const Container* container = nullptr;
    
    // get the range of elements with requested key
    std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const {
        assert(container != nullptr);
        auto [begin, end] = container->equal_range(surface.geometryId());
        return {Iterator{begin}, Iterator{end}};
    }
};

template <typename propagator_t, typename candidate_container_t>
class CKFTrackFindingAlgorithm : public IAlgorithm {
    public:
        using SimpleSourceLinkContainer =
            std::unordered_multimap<
                Acts::GeometryIdentifier, 
                SimpleSourceLink>;
        using SimpleSourceLinkAccessor = SimpleContainerAccessor<SimpleSourceLinkContainer>;
        using CombinatorialKalmanFilterOptions =
            Acts::CombinatorialKalmanFilterOptions<
                SimpleSourceLinkAccessor::Iterator,
                candidate_container_t>;
        
        /// @brief The nested configuration struct
        struct Config {
            /// CKF finder
            const Acts::CombinatorialKalmanFilter<propagator_t, candidate_container_t>& ckf;
            /// CKF extensions
            Acts::CombinatorialKalmanFilterExtensions<candidate_container_t> extensions;
            /// The input collection
            std::string inputSeeds = "Seed";
            /// The output collection before filtering
            std::string outputTrackCandidates = "TrackCandidates";
            /// Minimum number of source links
            int minSourceLinks = 3;
            /// Maximum number of source links
            int maxSourceLinks = 10;
        };

        /// @brief Constructor
        CKFTrackFindingAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("CKFTrackFindingAlgorithm", level),
            m_cfg(std::move(config)) {
                m_inputSeeds.initialize(m_cfg.inputSeeds);
                m_outputTrackCandidates.initialize(m_cfg.outputTrackCandidates);
        }
        ~CKFTrackFindingAlgorithm() = default;

        CKFTrackFindingAlgorithm(Config config) : m_cfg(config) {}

        /// @brief The execute method
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // Get the input seeds
            // from the context
            auto input = m_inputSeeds(ctx);

            auto options = CombinatorialKalmanFilterOptions(
                ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
                Acts::SourceLinkAccessorDelegate<SimpleSourceLinkAccessor::Iterator>{},
                m_cfg.extensions, Acts::PropagatorPlainOptions(
                    ctx.geoContext, ctx.magFieldContext));

            SimpleSourceLinkAccessor slAccessor;
            options.sourcelinkAccessor.template connect<&SimpleSourceLinkAccessor::range>(
                &slAccessor);

            Seeds trackCandidates;
            for (const auto& seed : input) {
                SimpleSourceLinkContainer ckfSourceLinks;
                for (auto& sl : seed.sourceLinks) {
                    auto ssl = sl.get<SimpleSourceLink>();
                    ckfSourceLinks.insert({ssl.geometryId(), ssl});
                }
    
                slAccessor.container = &ckfSourceLinks;

                Acts::TrackContainer tc{
                    Acts::VectorTrackContainer{},
                    Acts::VectorMultiTrajectory{}};

                // run the CKF for all initial track states
                Acts::CurvilinearTrackParameters ipParameters = seed.ipParameters;
        
                options.propagatorPlainOptions.maxSteps = 10000;
        
                auto res = m_cfg.ckf.findTracks(ipParameters, options, tc);
                if (!res.ok()) {
                    continue;
                }

                for (std::size_t tid = 0u; tid < tc.size(); ++tid) {
                    const auto track = tc.getTrack(tid);
                
                    // check purity of first found track
                    // find the number of hits not originating from the right track
                    std::vector<Acts::SourceLink> sourceLinks;
                    for (const auto trackState : track.trackStatesReversed()) {
                        if (!trackState.hasUncalibratedSourceLink()) {
                            continue;
                        }
                        Acts::SourceLink sl = trackState.getUncalibratedSourceLink();
                        sourceLinks.push_back(sl);
                    }

                    if (sourceLinks.size() < m_cfg.minSourceLinks ||
                        sourceLinks.size() > m_cfg.maxSourceLinks) {
                            continue;
                    }

                    trackCandidates.push_back(Seed{
                        sourceLinks,
                        ipParameters,
                        seed.trackId});
                }
            }

            m_outputTrackCandidates(ctx, std::move(trackCandidates));

            return ProcessCode::SUCCESS;
        }

    private:
        Config m_cfg;

        ReadDataHandle<Seeds> m_inputSeeds
            {this, "InputSeeds"};

        WriteDataHandle<Seeds> m_outputTrackCandidates
            {this, "OutputTrackCandidates"};
};
