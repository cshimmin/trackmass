#include <iostream>
#include <algorithm>

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"
#include "OutputTree/OutputTree.h"
#include "xAODEventInfo/EventInfo.h"

#include "trackmass/Retrieve.h"

#include "TFile.h"
#include "TH1F.h"

using namespace std;

bool pt_compare(const xAOD::IParticle* j1, const xAOD::IParticle* j2) {
  return j1->pt() > j2->pt();
}

float get_tau21(const xAOD::Jet* jet) {
  float tau1 = jet->auxdata<float>("Tau1_wta");
  float tau2 = jet->auxdata<float>("Tau2_wta");
  return tau1!=0 ? tau2/tau1 : -1;
}

const xAOD::Vertex* get_pv(xAOD::TEvent* evt) {
  const xAOD::VertexContainer* vtxs(nullptr);
  evt->retrieve(vtxs, "PrimaryVertices");
  for (auto vtx : *vtxs) {
    if (vtx->vertexType() == xAOD::VxType::PriVtx) {
      return vtx;
    }
  }
  cout << "Warning! No primary vertex found." << endl;
  return nullptr;
}

bool is_pv_track(const xAOD::TrackParticle* trk, const xAOD::Vertex* pv) {
  if (trk->pt() < 500) return false;
  if (trk->vertex() == pv) return true;

  if (fabs((trk->z0()+trk->vz()-pv->z())*sin(trk->theta()))<3.) return true;

  return false;
}

void process_event(xAOD::TEvent* evt, OutputTree* out, bool verbose) {
  const xAOD::EventInfo* event_info(nullptr);
  evt->retrieve(event_info, "EventInfo");

  const xAOD::Vertex* pv = get_pv(evt);

  float wt = event_info->mcEventWeight();
  out->add_scalar("evt_wt", wt);

  auto fatjets = CopyRetrieve<xAOD::Jet>(evt, "AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets");
  if (verbose) {
    for (auto j : fatjets) {
      cout << "Leading jet pT = " << j->pt() << endl;
      break;
    }
  }
  sort(begin(fatjets), end(fatjets), pt_compare);

  // calculate the TA mass for all the fatjets
  vector<float> fj_mTA;
  vector<float> fj_ntrk;
  for (auto j : fatjets) {
    std::vector<const xAOD::TrackParticle*> ghost_tracks;
    j->getAssociatedObjects<xAOD::TrackParticle>(xAOD::JetAttribute::GhostTrack, ghost_tracks);
    TLorentzVector tsum;
    int ntrk = 0;
    for (auto trk : ghost_tracks) {
      if (! is_pv_track(trk, pv)) continue;
      tsum += trk->p4();
      ntrk++;
    }
    if (ntrk > 0) {
      fj_mTA.push_back( tsum.M() / tsum.Pt() * j->pt() );
    } else {
      fj_mTA.push_back(-1);
    }
    fj_ntrk.push_back(ntrk);
  }

  // calculate the jet tau21
  vector<float> fj_tau21;
  transform(begin(fatjets), end(fatjets), back_inserter(fj_tau21), get_tau21);

  if (verbose) {
    for (float x : fj_tau21) {
      cout << "leading tau21: " << x << endl;
      break;
    }
  }
  
  out->add_jets("fj", fatjets);
  out->add_vector("fj_mTA", fj_mTA);
  out->add_vector("fj_ntrk", fj_ntrk);
  out->add_vector("fj_tau21", fj_tau21);

  // grab the small-radius jets
  auto jets = CopyRetrieve<xAOD::Jet>(evt, "AntiKt4LCTopoJets");
  sort(begin(jets), end(jets), pt_compare);

  out->add_jets("jet", jets);

  if (fatjets.size() > 0) {
    // identify the leading AntiKt10 jet
    const xAOD::Jet* fj1 = fatjets[0];
    const xAOD::Jet* fj2 = fj1;
    int idx1 = 0;
    int idx2 = 0;

    // identify the leading AntiKt10 jet with smaller tau21
    if (fatjets.size() > 1 and fj_tau21[1] < fj_tau21[0]) {
      fj2 = fatjets[1];
      idx2 = 1;
    }

    out->add_jet("res1", fj1);
    out->add_vector("res1_mTA", {fj_mTA[idx1]});
    out->add_vector("res1_ntrk", {fj_ntrk[idx1]});
    out->add_vector("res1_tau21", {fj_tau21[idx1]});
    out->add_jet("res2", fj2);
    out->add_vector("res2_mTA", {fj_mTA[idx2]});
    out->add_vector("res2_ntrk", {fj_ntrk[idx2]});
    out->add_vector("res2_tau21", {fj_tau21[idx2]});

    // find the leading small radius jet that is not overlapping with this large-R jet
    for (auto j : jets) {
      if (fj1->p4().DeltaR(j->p4()) > 1.0) {
        out->add_jet("probe1", j);
        break;
      }
    }
    for (auto j : jets) {
      if (fj2->p4().DeltaR(j->p4()) > 1.0) {
        out->add_jet("probe2", j);
        break;
      }
    }
  }

  out->Fill();
}

int main(int argc, char** argv) {
  string output_filename;
  vector<string> input_files;

  if (argc < 3) {
    cerr << "Please give a filename!" << endl;
    cout << "Usage: " << argv[0] << " output_file input_file,[input_file,...]" << endl;
    return 1;
  }

  output_filename = argv[1];

  string firstfile = argv[2];
  if (firstfile.find(",") != string::npos) {
    // parse comma-separated file list
    size_t pos = 0;
    while ( (pos = firstfile.find(",")) != string::npos) {
      input_files.push_back(firstfile.substr(0, pos));
      firstfile.erase(0, pos + 1);
    }
    if (firstfile.size() > 0) {
      input_files.push_back(firstfile);
    }
  }
  else {
    // read space-delimited file list
    for (int i = 2; i < argc; ++i) {
      input_files.push_back(argv[i]);
    }
  }

  cout << "Listing input files:" << endl;
  for (auto s : input_files) {
    cout << "    " << s << endl;
  }

  cout << "Initializing xAOD..." << endl;
  xAOD::Init();

  TFile* out_file = new TFile(output_filename.c_str(), "recreate");
  out_file->cd();
  OutputTree* out_tree = new OutputTree("tracks");
  TH1F* h_nevt_total = new TH1F("nevt_total", "nevt_total", 1, 0, 1);
  TH1F* h_nevt_total_wt = new TH1F("nevt_total_wt", "nevt_total_wt", 1, 0, 1);

  cout << "Creating TEvent..." << endl;
  xAOD::TEvent *evt = new xAOD::TEvent();

  int itotal = 0;
  for (auto fname : input_files) {
    cout << "Opening file: " << fname << endl;
    TFile* input_file = TFile::Open(fname.c_str());

    cout << "Connecting TEvent to file..." << endl;
    evt->readFrom(input_file);

    cout << "Grabbing MetaData..." << endl;

    const xAOD::CutBookkeeperContainer* cuts(nullptr);
    evt->retrieveMetaInput(cuts, "CutBookkeepers");
    cout << "Found " << cuts->size() << " cutbookkeepers." << endl;
    float nevt_total = -1;
    float nevt_total_wt = -1;
    for (auto cb : *cuts) {
      if (cb->inputStream() != "StreamDAOD_JETM8") continue;
      if (cb->name() != "AllExecutedEvents") continue;
      nevt_total = cb->nAcceptedEvents();
      nevt_total_wt = cb->sumOfEventWeights();
    }
    if ( (nevt_total<0) || (nevt_total_wt<0) ) {
      cerr << "Invalid event counts found!!! Skipping file." << endl;
      continue;
    }
    h_nevt_total->Fill(0., nevt_total);
    h_nevt_total_wt->Fill(0., nevt_total_wt);
    cout << "Total events:    " << nevt_total << endl;
    cout << "Weighted events: " << nevt_total_wt << endl;

    int nevt = evt->getEntries();
    cout << nevt << " entries in file." << endl;
    for (int i = 0; i < nevt; ++i) {
      bool verbose = itotal%100==0;
      itotal++;
      if (verbose) {
        cout << "======= Event " << itotal << " =======" << endl;
      }

      out_tree->clear();
      evt->getEntry(i);
      process_event(evt, out_tree, verbose);

      if (verbose) {
        cout << "========================" << endl << endl;
      }
    }

    input_file->Close();
  }

  out_file->cd();
  out_tree->Write();
  h_nevt_total->Write();
  h_nevt_total_wt->Write();

  out_file->Close();

  return 0;
}
