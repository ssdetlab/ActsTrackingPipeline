#include "TrackingPipeline/Io/RootMaterialTrackReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>
#include <mutex>
#include <stdexcept>

#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"

RootMaterialTrackReader::RootMaterialTrackReader(const Config& config,
                                                 Acts::Logging::Level level)
    : IReader(),
      m_logger(Acts::getDefaultLogger(name(), level)),
      m_cfg(config) {
  if (m_cfg.filePaths.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  if (m_cfg.filePaths.size() == 1) {
    m_file = new TFile(m_cfg.filePaths.at(0).c_str());
    m_tree = m_file->Get<TTree>(m_cfg.treeName.c_str());
  } else {
    m_chainOwner = new TChain(m_cfg.treeName.c_str());
    // Add the files to the chain
    for (const auto& path : m_cfg.filePaths) {
      m_chainOwner->Add(path.c_str());
    }
    m_tree = dynamic_cast<TTree*>(m_chainOwner);
  }
  // Set event Id branch
  m_tree->SetBranchAddress("event_id", &m_eventId);
  if (m_tree->GetBranch("event_id") == nullptr) {
    throw std::invalid_argument("Missing eventId branch");
  }
  auto nEntries = static_cast<std::size_t>(m_tree->GetEntries());

  // Add the first entry
  m_tree->GetEntry(0);
  m_eventMap.emplace_back(m_eventId, 0, 0);

  for (std::size_t i = 0; i < nEntries; ++i) {
    m_tree->GetEntry(i);
    if (m_eventId != std::get<0>(m_eventMap.back())) {
      std::get<2>(m_eventMap.back()) = i;
      m_eventMap.emplace_back(m_eventId, i, i);
    }
    if (i == nEntries - 1) {
      std::get<2>(m_eventMap.back()) = nEntries;
    }
  }

  // Sort by event id
  std::ranges::sort(m_eventMap, [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  ACTS_DEBUG("Event range: " << availableEvents().first << " - "
                             << availableEvents().second);

  //------------------------------------------------------------------
  // Set the rest of the branches

  // Start global x
  m_tree->SetBranchAddress("v_x", &m_v_x);
  // Start global y
  m_tree->SetBranchAddress("v_y", &m_v_y);
  // Start global z
  m_tree->SetBranchAddress("v_z", &m_v_z);
  // Start global momentum x
  m_tree->SetBranchAddress("v_px", &m_v_px);
  // Start global momentum y
  m_tree->SetBranchAddress("v_py", &m_v_py);
  // Start global momentum z
  m_tree->SetBranchAddress("v_pz", &m_v_pz);
  // Start phi direction
  m_tree->SetBranchAddress("v_phi", &m_v_phi);
  // Start eta direction
  m_tree->SetBranchAddress("v_eta", &m_v_eta);
  // Thickness in X0/L0
  m_tree->SetBranchAddress("t_X0", &m_tX0);
  // Thickness in X0/L0
  m_tree->SetBranchAddress("t_L0", &m_tL0);

  // Step x position
  m_tree->SetBranchAddress("mat_x", &m_step_x);
  // Step y position
  m_tree->SetBranchAddress("mat_y", &m_step_y);
  // Step z position
  m_tree->SetBranchAddress("mat_z", &m_step_z);
  // Step x direction
  m_tree->SetBranchAddress("mat_dx", &m_step_dx);
  // Step y direction
  m_tree->SetBranchAddress("mat_dy", &m_step_dy);
  // Step z direction
  m_tree->SetBranchAddress("mat_dz", &m_step_dz);
  // Step length
  m_tree->SetBranchAddress("mat_step_length", &m_step_length);
  // Step material x0
  m_tree->SetBranchAddress("mat_X0", &m_step_X0);
  // Step material l0
  m_tree->SetBranchAddress("mat_L0", &m_step_L0);
  // Step material A
  m_tree->SetBranchAddress("mat_A", &m_step_A);
  // Step material Z
  m_tree->SetBranchAddress("mat_Z", &m_step_Z);
  // Step material rho
  m_tree->SetBranchAddress("mat_rho", &m_step_rho);

  //------------------------------------------------------------------
  // Initialize the data handles
  m_outputMaterialTracks.initialize(m_cfg.outputMaterialTracks);
}

std::string RootMaterialTrackReader::name() const {
  return "RootMaterialTrackReader";
}

std::pair<std::size_t, std::size_t> RootMaterialTrackReader::availableEvents()
    const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode RootMaterialTrackReader::read(const AlgorithmContext& ctx) {
  auto it = std::ranges::find_if(m_eventMap, [&](const auto& a) {
    return std::get<0>(a) == ctx.eventNumber;
  });

  if (it == m_eventMap.end()) {
    // Explicitly warn if it happens for the first or last event as that might
    // indicate a human error
    if ((ctx.eventNumber == availableEvents().first) &&
        (ctx.eventNumber == availableEvents().second - 1)) {
      ACTS_WARNING("Reading empty event: " << ctx.eventNumber);
    } else {
      ACTS_DEBUG("Reading empty event: " << ctx.eventNumber);
    }

    m_outputMaterialTracks(ctx, {});

    // Return success flag
    return ProcessCode::SUCCESS;
  }

  // Lock the mutex
  std::lock_guard<std::mutex> lock(m_mutex);

  ACTS_DEBUG("Reading event: " << std::get<0>(*it)
                               << " stored in entries: " << std::get<1>(*it)
                               << " - " << std::get<2>(*it));

  // Initialize collection
  std::vector<Acts::RecordedMaterialTrack> mTrackCollection;

  // Loop over the entries for this event
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_tree->GetEntry(entry);

    Acts::RecordedMaterialTrack rmTrack;

    // Fill the position and momentum
    rmTrack.first.first =
        m_cfg.toWorldTransform * Acts::Vector3(m_v_x, m_v_y, m_v_z);
    rmTrack.first.second =
        m_cfg.toWorldTransform * Acts::Vector3(m_v_px, m_v_py, m_v_pz);

    // Fill the individual steps
    std::size_t mSteps = m_step_length->size();
    rmTrack.second.materialInteractions.reserve(mSteps);
    rmTrack.second.materialInX0 = 0.;
    rmTrack.second.materialInL0 = 0.;

    for (std::size_t is = 0; is < mSteps; ++is) {
      double s = m_step_length->at(is);
      if (s == 0) {
        continue;
      }

      double mX0 = m_step_X0->at(is);
      double mL0 = m_step_L0->at(is);

      rmTrack.second.materialInX0 += s / mX0;
      rmTrack.second.materialInL0 += s / mL0;

      // Fill the position & the material
      Acts::MaterialInteraction mInteraction;
      mInteraction.position =
          m_cfg.toWorldTransform *
          Acts::Vector3(m_step_x->at(is), m_step_y->at(is), m_step_z->at(is));
      mInteraction.direction =
          m_cfg.toWorldTransform * Acts::Vector3(m_step_dx->at(is),
                                                 m_step_dy->at(is),
                                                 m_step_dz->at(is));
      mInteraction.materialSlab = Acts::MaterialSlab(
          Acts::Material::fromMassDensity(mX0, mL0, m_step_A->at(is),
                                          m_step_Z->at(is), m_step_rho->at(is)),
          s);

      if (m_cfg.readCachedSurfaceInformation) {
        // Add the surface information to the interaction this allows the
        // mapping to be speed up
        mInteraction.intersectionID =
            Acts::GeometryIdentifier(m_sur_id->at(is));
        mInteraction.intersection =
            m_cfg.toWorldTransform *
            Acts::Vector3(m_sur_x->at(is), m_sur_y->at(is), m_sur_z->at(is));
        mInteraction.pathCorrection = m_sur_pathCorrection->at(is);
      } else {
        mInteraction.intersectionID = Acts::GeometryIdentifier();
        mInteraction.intersection = Acts::Vector3(0, 0, 0);
      }
      rmTrack.second.materialInteractions.push_back(std::move(mInteraction));
    }
    mTrackCollection.push_back(std::move(rmTrack));
  }

  // Write to the collection to the EventStore
  m_outputMaterialTracks(ctx, std::move(mTrackCollection));

  // Return success flag
  return ProcessCode::SUCCESS;
}
