#ifdef __CLING__
#pragma link C++ class vector < TVector3> + ;
#pragma link C++ class vector < TVector2> + ;
#pragma link C++ class vector < TLorentzVector> + ;
#pragma link C++ class vector < TMatrixD> + ;
#pragma link C++ class tuple < int, int> + ;
#pragma link C++ class tuple < double, double> + ;

#pragma link C++ class pixel + ;
#pragma link C++ class chip + ;
#pragma link C++ class stave + ;
#pragma link C++ class event + ;
#pragma link C++ class run_meta_data + ;

#pragma link C++ class chip_event + ;
#pragma link C++ class vector < chip_event> + ;
#pragma link C++ class E320Io::ChipEvent::Cluster + ;
#pragma link C++ class E320Io::ChipEvent + ;
#pragma link C++ class vector < E320Io::ChipEvent> + ;
#pragma link C++ class tlu_event + ;
#pragma link C++ class EpicsFrame + ;
#pragma link C++ class stave_event + ;
#pragma link C++ class E320Io::StaveEvent + ;
#pragma link C++ class detector_event + ;
#pragma link C++ class E320Io::DetectorEvent + ;
#pragma link C++ class detector_event_tlu + ;
#pragma link C++ class ALPIDERegs + ;
#pragma link C++ class chip_run_meta + ;
#pragma link C++ class vector < chip_run_meta> + ;
#pragma link C++ class MOSAICRegs + ;
#pragma link C++ class stave_run_meta + ;
#pragma link C++ class std::vector < stave_run_meta> + ;
#pragma link C++ class run_meta + ;
#endif
