#pragma once

#include "ActsLUXEPipeline/IReader.hpp"
#include "ActsLUXEPipeline/ProcessCode.hpp"
#include "ActsLUXEPipeline/AlgorithmContext.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"

#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h" 

/// @brief Set branch addresses for the TChain
///
/// @tparam K type of the keys
/// @tparam T type of the columns
///
/// @param chain input TChain
/// @param keys keys to set the branches for
/// @param columns columns to set the branches for
template <typename K, typename T>
inline void setBranches(
    TChain* chain,
    const K& keys, 
    std::unordered_map<std::string_view, T>& columns) {
    T value = 0;
    for (auto key : keys) {
        columns.insert({key, value});
    }
    for (auto key : keys) {
        chain->SetBranchAddress(key, &columns.at(key));
    }
};

/// @brief Intermediate generalization of the 
/// ROOT file reader to be inhereted from by the
/// readers for the specific tree structures,
/// data types and geometries
///
/// @tparam measurementContainer_t container type 
/// for the measurements to be implemented 
///
/// @note The events are assumed to be ordered
template <typename measurementContainer_t>
class ROOTDataReader : public IReader {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// name of the whiteboard entry
            std::string dataCollection;
            /// Name of the input tree
            std::string treeName = "clusters";
            /// The names of the input files
            std::vector<std::string> filePaths;
            /// The keys we have in the ROOT file
            std::vector<const char*> vector3Keys;
            std::vector<const char*> lorentzKeys;
            std::vector<const char*> vectorIntKeys;
            std::vector<const char*> vectorDoubleKeys;
            std::vector<const char*> intKeys;
        };

        ROOTDataReader(const ROOTDataReader &) = delete;
        ROOTDataReader(const ROOTDataReader &&) = delete;
    
        /// Constructor
        /// @param config The Configuration struct
        ROOTDataReader(const Config &config, Acts::Logging::Level level)
            : IReader(),
            m_cfg(config),
            m_logger(Acts::getDefaultLogger(name(), level)) {
                m_chain = new TChain(m_cfg.treeName.c_str());

                if (m_cfg.filePaths.empty()) {
                    throw std::invalid_argument("Missing input filenames");
                }
                if (m_cfg.treeName.empty()) {
                    throw std::invalid_argument("Missing tree name");
                }
        
                m_outputData.initialize(m_cfg.dataCollection);
        
                // Set the branches
                setBranches(m_chain, m_cfg.vector3Keys, m_vector3Columns);
                setBranches(m_chain, m_cfg.lorentzKeys, m_lorentzColumns);
                setBranches(m_chain, m_cfg.vectorIntKeys, m_vectorIntColumns);
                setBranches(m_chain, m_cfg.vectorDoubleKeys, m_vectorDoubleColumns);
                setBranches(m_chain, m_cfg.intKeys, m_intColumns);

                // Add the files to the chain
                for (auto path : m_cfg.filePaths) {
                    m_chain->Add(path.c_str());
                }

                // Disable all branches and only enable event-id for a first scan of the file
                m_chain->SetBranchStatus("*", false);
                if(!m_chain->GetBranch("eventId")) {
                    throw std::invalid_argument("Missing eventId branch");
                }
                m_chain->SetBranchStatus("eventId", true);
                
                auto nEntries = static_cast<std::size_t>(m_chain->GetEntries());

                // Add the first entry
                m_chain->GetEntry(0);
                m_eventMap.push_back({m_intColumns.at("eventId"), 0ul, 0ul});
    
                // Go through all entries and store the position of the events
                for (auto i = 1ul; i < nEntries; ++i) {
                    m_chain->GetEntry(i);
                    const auto evtId = m_intColumns.at("eventId");
                
                    if (evtId != std::get<0>(m_eventMap.back())) {
                        std::get<2>(m_eventMap.back()) = i;
                        m_eventMap.push_back({evtId, i, i});
                    }
                }
                // Sort by event id
                std::sort(m_eventMap.begin(), m_eventMap.end(),
                    [] (const auto& a, const auto& b) {
                        return std::get<0>(a) < std::get<0>(b);
                    }
                );
        
                std::get<2>(m_eventMap.back()) = nEntries;
        
                // Re-Enable all branches
                m_chain->SetBranchStatus("*", true);
                ACTS_DEBUG("Event range: " << availableEvents().first << " - "
                    << availableEvents().second);
        }
    
        /// Reader name() method
        virtual std::string name() const { return "ROOTDataReader"; }
    
        /// Return the available events range.
        std::pair<std::size_t, std::size_t> 
            availableEvents() const override {
                return {std::get<0>(m_eventMap.front()), 
                    std::get<0>(m_eventMap.back()) + 1};
        }
    
        /// Read out data from the input stream
        ProcessCode read(const AlgorithmContext &context) override {
            auto it = std::find_if(
                m_eventMap.begin(), m_eventMap.end(),
                [&](const auto& a) { return std::get<0>(a) == context.eventNumber; });
        
                if (it == m_eventMap.end()) {
                    // explicitly warn if it happens for the first or last event as that might
                    // indicate a human error
                    if ((context.eventNumber == availableEvents().first) &&
                        (context.eventNumber == availableEvents().second - 1)) {
                            ACTS_WARNING("Reading empty event: " << context.eventNumber);
                    } else {
                        ACTS_DEBUG("Reading empty event: " << context.eventNumber);
                    }
                
                    m_outputData(context, {});
                
                    // Return success flag
                    return ProcessCode::SUCCESS;
                }
                
                // lock the mutex
                std::lock_guard<std::mutex> lock(m_read_mutex);
                
                ACTS_DEBUG("Reading event: " << std::get<0>(*it)
                    << " stored in entries: " << std::get<1>(*it)
                    << " - " << std::get<2>(*it));
        
                // Create the measurements
                measurementContainer_t measurements;
                for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
                    m_chain->GetEntry(entry);
                    prepareMeasurements(context, &measurements);
                }
                
                m_outputData(context, std::move(measurements));
                
                // Return success flag
                return ProcessCode::SUCCESS;
        }
    
        /// Readonly access to the config
        const Config &config() const { return m_cfg; }

    private:
        /// Private access to the logging instance
        const Acts::Logger &logger() const { return *m_logger; }

        /// The config class
        Config m_cfg;

        /// Prepare the measurements
        virtual inline void prepareMeasurements(
            const AlgorithmContext &context, 
            measurementContainer_t* measurements) const = 0;

        WriteDataHandle<measurementContainer_t> m_outputData{this, "OutputData"};
        std::unique_ptr<const Acts::Logger> m_logger;

        /// mutex used to protect multi-threaded reads
        std::mutex m_read_mutex;

        /// Vector of {eventNr, entryMin, entryMax}
        std::vector<std::tuple<uint32_t, std::size_t, std::size_t>> m_eventMap;

        /// The input tree name
        TChain *m_chain = nullptr;

    protected:
        /// The exausitive list of columns
        std::unordered_map<std::string_view, 
            std::vector<std::int32_t>*> m_vectorIntColumns;
        std::unordered_map<std::string_view, 
            std::vector<std::double_t>*> m_vectorDoubleColumns;
        std::unordered_map<std::string_view, 
            std::int32_t> m_intColumns;
        std::unordered_map<std::string_view, 
            std::vector<TVector3>*> m_vector3Columns;
        std::unordered_map<std::string_view, 
            std::vector<TLorentzVector>*> m_lorentzColumns;
};
