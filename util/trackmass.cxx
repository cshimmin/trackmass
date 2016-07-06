#include <iostream>

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/PhotonContainer.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/tools/Filter.hh>

#include "TFile.h"

#include "trackmass/OutputTree.h"

using std::cout;
using std::endl;

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

void process_event(xAOD::TEvent* evt, OutputTree* out) {
  //const xAOD::EventInfo* event_info = nullptr;
  //evt->retrieve(event_info, "EventInfo");
  
  
  const xAOD::PhotonContainer* photons(nullptr);
  evt->retrieve(photons, "Photons");

  std::vector<const xAOD::Photon*> good_photons, trig_photons;

  for (auto ph : *photons) {
    if (ph->pt() < 15e3) continue;
    good_photons.push_back(ph);
    out->add_photon("phgood", ph->p4());
    if (ph->pt() < 100e3) continue;
    trig_photons.push_back(ph);
    out->add_photon("ph", ph->p4());
  }
  if (trig_photons.size() < 1) { return; }
  
  const xAOD::Photon* ph0 = trig_photons.at(0);
  const TLorentzVector ph0_v = ph0->p4();

  auto jet_overlaps_photon = [&](const xAOD::Jet* j) { return ph0_v.DeltaR(j->p4()) < 1; };
  
  const xAOD::Vertex* pv = get_pv(evt);

  const xAOD::TrackParticleContainer* tracks(nullptr);
  evt->retrieve(tracks, "InDetTrackParticles");

  std::vector<const xAOD::TrackParticle*> pv_tracks;
  for (auto trk : *tracks) {
    if (!is_pv_track(trk, pv)) continue;
    pv_tracks.push_back(trk);
  }

  if (pv_tracks.size() < 1) return;
  cout << "There are " << pv_tracks.size() << " PV tracks" << endl;

  TLorentzVector trksum;
  std::vector<fastjet::PseudoJet> pseudojets;
  for (auto trk : pv_tracks) {
    TLorentzVector v = trk->p4();
    trksum += v;
    pseudojets.push_back(fastjet::PseudoJet(v.Px(), v.Py(), v.Pz(), v.E()));
  }
  cout << "trksum pt=" << trksum.Pt() << "\tm=" << trksum.M() << endl;
  out->add_jet("trksum", trksum);
  
  static fastjet::JetDefinition jet_def02CA(fastjet::cambridge_algorithm, 0.2);
  static fastjet::JetDefinition jet_def02(fastjet::antikt_algorithm, 0.2);
  static fastjet::JetDefinition jet_def10(fastjet::antikt_algorithm, 1.0);
  static fastjet::JetDefinition jet_def15(fastjet::antikt_algorithm, 1.5);
  static fastjet::JetDefinition jet_def60(fastjet::antikt_algorithm, 6.0);

  fastjet::ClusterSequence cs02CA(pseudojets, jet_def02CA);
  fastjet::ClusterSequence cs02(pseudojets, jet_def02);
  fastjet::ClusterSequence cs10(pseudojets, jet_def10);
  fastjet::ClusterSequence cs15(pseudojets, jet_def15);
  fastjet::ClusterSequence cs60(pseudojets, jet_def60);
  
  std::vector<fastjet::PseudoJet> CAtrkjets02 = fastjet::sorted_by_pt(cs02CA.inclusive_jets());
  std::vector<fastjet::PseudoJet> trkjets02 = fastjet::sorted_by_pt(cs02.inclusive_jets());
  std::vector<fastjet::PseudoJet> trkjets10 = fastjet::sorted_by_pt(cs10.inclusive_jets());
  std::vector<fastjet::PseudoJet> trkjets15 = fastjet::sorted_by_pt(cs15.inclusive_jets());
  std::vector<fastjet::PseudoJet> trkjets60 = fastjet::sorted_by_pt(cs60.inclusive_jets());

  std::vector<fastjet::PseudoJet> CAtrkjets02_OR;
  int nremoved = 0;
  for (auto j : CAtrkjets02) {
    TLorentzVector v;
    v.SetPxPyPzE(j.px(), j.py(), j.pz(), j.e());
    if (ph0_v.DeltaR(v) < 0.2) continue;
    CAtrkjets02_OR.push_back(j);
    nremoved++;
  }
  cout << "Removed " << nremoved << " subjets near the photon." << endl;

  fastjet::ClusterSequence recluster10(CAtrkjets02, jet_def10);
  fastjet::ClusterSequence recluster15(CAtrkjets02, jet_def15);
  fastjet::ClusterSequence recluster10_OR(CAtrkjets02_OR, jet_def10);
  fastjet::ClusterSequence recluster15_OR(CAtrkjets02_OR, jet_def15);
  fastjet::ClusterSequence recluster60_OR(CAtrkjets02_OR, jet_def60);
  std::vector<fastjet::PseudoJet> recjets10 = fastjet::sorted_by_pt(recluster10.inclusive_jets());
  std::vector<fastjet::PseudoJet> recjets15 = fastjet::sorted_by_pt(recluster15.inclusive_jets());
  std::vector<fastjet::PseudoJet> recjets10_OR = fastjet::sorted_by_pt(recluster10_OR.inclusive_jets());
  std::vector<fastjet::PseudoJet> recjets15_OR = fastjet::sorted_by_pt(recluster15_OR.inclusive_jets());
  std::vector<fastjet::PseudoJet> recjets60_OR = fastjet::sorted_by_pt(recluster60_OR.inclusive_jets());

  out->add_jets("rec10", recjets10);
  out->add_jets("rec15", recjets15);

  out->add_jets("rec10OR", recjets10_OR);
  out->add_jets("rec15OR", recjets15_OR);
  out->add_jets("rec60OR", recjets60_OR);

  static fastjet::Filter trim_filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05));

  for (auto j : trkjets10) {
    fastjet::PseudoJet filtered_jet = trim_filter(j);
    out->add_jet("tj10_trim", filtered_jet);
  }
  for (auto j : trkjets60) {
    fastjet::PseudoJet filtered_jet = trim_filter(j);
    out->add_jet("tj60_trim", filtered_jet);
  }
  for (auto j : recjets60_OR) {
    fastjet::PseudoJet filtered_jet = trim_filter(j);
    out->add_jet("rec60OR_trim", filtered_jet);
  }

  cout << "found " << trkjets02.size() << " akt02 jets" << endl;
  cout << "found " << trkjets10.size() << " akt10 jets" << endl;
  cout << "found " << trkjets15.size() << " akt15 jets" << endl;
  cout << "found " << trkjets60.size() << " akt60 jets" << endl;

  if (trkjets10.size() > 0 && trkjets15.size() > 0 && trkjets60.size() > 0) {
    cout << "m10=" << trkjets10.at(0).m() << "\tm15=" << trkjets15.at(0).m() << "\tm60=" << trkjets60.at(0).m() << endl;
  }

  TLorentzVector sum02;
  for (auto j : trkjets02) {
    out->add_jet("tj02", j);
    TLorentzVector v;
    v.SetPxPyPzE(j.px(), j.py(), j.pz(), j.e());
    sum02 += v;
  }
  out->add_jet("sum02", sum02);
  for (auto j : trkjets10) {
    out->add_jet("tj10", j);
  }
  for (auto j : trkjets15) {
    out->add_jet("tj15", j);
  }
  for (auto j : trkjets60) {
    out->add_jet("tj60", j);
  }

  /*
  for (auto j : trkjets15) {
    out->j15_pt.push_back(j.pt());
    out->j15_m.push_back(j.pt());
  }
  for (auto j : trkjets60) {
    out->j60_pt.push_back(j.pt());
    out->j60_m.push_back(j.pt());
  }
  */

  
  const xAOD::JetContainer* fatjets(nullptr);
  evt->retrieve(fatjets, "AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets");

  std::vector<const xAOD::Jet*> fatjetsOR;
  std::remove_copy_if(fatjets->begin(), fatjets->end(), std::back_inserter(fatjetsOR), jet_overlaps_photon);
  out->add_jets("fatjetsOR", fatjetsOR);

  cout << "There are " << fatjets->size() << " fatjets" << endl;
  int nkeep = 0;
  for (auto j : *fatjets) {
    if (ph0_v.DeltaR(j->p4()) < 0.5) { continue; }
    out->add_jet("fj", j->p4());
    nkeep++;

    std::vector<const xAOD::TrackParticle*> ghosttracks;
    j->getAssociatedObjects<xAOD::TrackParticle>(xAOD::JetAttribute::GhostTrack,ghosttracks);
    TLorentzVector tsum;
    int ntrk = 0;
    for (auto trk : ghosttracks) {
      if (! is_pv_track(trk, pv)) continue;
      tsum += trk->p4();
      ntrk++;
    }
    out->add_jet("gaNolim", tsum);
    if (ntrk>1) {
      out->add_jet("ga", tsum);
    }
  }
  cout << "Kept " << nkeep << " fatjets." << endl;

  const xAOD::JetContainer* fatjets_truth(nullptr);
  evt->retrieve(fatjets_truth, "AntiKt10TruthTrimmedPtFrac5SmallR20Jets");
  for (auto j : *fatjets_truth) {
    if (ph0_v.DeltaR(j->p4()) < 0.5) continue;
    out->add_jet("fjtruth", j->p4());
  }

  const xAOD::JetContainer* trackjets10(nullptr);
  evt->retrieve(trackjets10, "AntiKt10PV0TrackJets");
  out->add_jets("trackjets10", trackjets10);


  std::vector<const xAOD::Jet*> trackjets10OR;
  std::remove_copy_if(trackjets10->begin(), trackjets10->end(), std::back_inserter(trackjets10OR), jet_overlaps_photon);

  out->add_jets("trackjets10OR", trackjets10OR);

  const xAOD::JetContainer* trackjets10_trimmed(nullptr);
  evt->retrieve(trackjets10_trimmed, "AntiKt10PV0TrackTrimmedPtFrac5SmallR20Jets");
  out->add_jets("trackjets10_trimmed", trackjets10_trimmed);

  std::vector<const xAOD::Jet*> trackjets10_trimmedOR;
  std::remove_copy_if(trackjets10_trimmed->begin(), trackjets10_trimmed->end(), std::back_inserter(trackjets10_trimmedOR), jet_overlaps_photon);

  out->add_jets("trackjets10_trimmedOR", trackjets10_trimmedOR);

  const xAOD::CaloClusterContainer* clusters(nullptr);
  evt->retrieve(clusters, "CaloCalTopoClusters");
  cout << "There are " << clusters->size() << " clusters." << endl;
  std::vector<const xAOD::CaloCluster*> ta_clusters;
  std::vector<fastjet::PseudoJet> ta_clusters_pj;
  std::map<const xAOD::TrackParticle*, const xAOD::CaloCluster*> track_clusters;
  TLorentzVector cluster_sum, cluster_sum_OR;
  for (auto cl : *clusters) {
    const TLorentzVector cl_v = cl->p4();
    for (auto trk : pv_tracks) {
      if (cl_v.DeltaR(trk->p4()) < 0.1) {
        ta_clusters.push_back(cl);
        ta_clusters_pj.push_back(fastjet::PseudoJet(cl_v.Px(), cl_v.Py(), cl_v.Pz(), cl_v.E()));
        cluster_sum += cl_v;
        if (cl_v.DeltaR(ph0_v) > 0.2) {
          cluster_sum_OR += cl_v;
        }
        if (track_clusters.count(trk)) {
          //cout << "WARNING! this track already had a cluster associated. overwriting." << endl;
        }
        track_clusters[trk] = cl;
        break;
      }
    }
  }
  cout << "Found " << ta_clusters.size() << " clusters near tracks." << endl;

  out->add_jet("clustersum", cluster_sum);
  out->add_jet("clustersumOR", cluster_sum_OR);

  fastjet::ClusterSequence cs_ta10(ta_clusters_pj, jet_def10);
  fastjet::ClusterSequence cs_ta15(ta_clusters_pj, jet_def15);
  std::vector<fastjet::PseudoJet> tacluster10_jets = fastjet::sorted_by_pt(cs_ta10.inclusive_jets());
  std::vector<fastjet::PseudoJet> tacluster15_jets = fastjet::sorted_by_pt(cs_ta15.inclusive_jets());
  std::remove_if(tacluster10_jets.begin(), tacluster10_jets.end(), [&](const fastjet::PseudoJet& j) { TLorentzVector v; v.SetPxPyPzE(j.px(), j.py(), j.pz(), j.e()); return ph0_v.DeltaR(v) < 1.0; });
  std::remove_if(tacluster15_jets.begin(), tacluster15_jets.end(), [&](const fastjet::PseudoJet& j) { TLorentzVector v; v.SetPxPyPzE(j.px(), j.py(), j.pz(), j.e()); return ph0_v.DeltaR(v) < 1.0; });
  out->add_jets("tacluster10", tacluster10_jets);
  out->add_jets("tacluster15", tacluster15_jets);

  std::vector<fastjet::PseudoJet> trkscaled;
  TLorentzVector trkscaled_sum;
  for (auto p : track_clusters) {
    TLorentzVector v;
    v.SetPtEtaPhiM(p.second->pt(), p.first->eta(), p.first->phi(), 0);
    trkscaled.push_back(fastjet::PseudoJet(v.Px(), v.Py(), v.Pz(), v.E()));
    trkscaled_sum += v;
  }
  out->add_jet("trkscaled_sum", trkscaled_sum);

  std::vector<TLorentzVector> boosted_tracks, unboosted_tracks, rboosted_tracks;
  for (auto trk : pv_tracks) {
    boosted_tracks.push_back(trk->p4() - ph0_v);
    rboosted_tracks.push_back(trk->p4() + ph0_v);
    unboosted_tracks.push_back(trk->p4());
  }
  out->add_jets("boosted_tracks", boosted_tracks);
  out->add_jets("rboosted_tracks", rboosted_tracks);
  out->add_jets("unboosted_tracks", unboosted_tracks);

  const xAOD::JetContainer* empflow(nullptr);
  const xAOD::JetContainer* emcpflow(nullptr);
  evt->retrieve(empflow, "AntiKt10EMPFlowTrimmedPtFrac5SmallR20Jets");
  evt->retrieve(emcpflow, "AntiKt10EMCPFlowTrimmedPtFrac5SmallR20Jets");

  //xAOD::JetContainer empflowOR, emcpflowOR;
  std::vector<const xAOD::Jet*> empflowOR, emcpflowOR;
  //for (auto j : *empflow) { empflowOR.push_back(j); }
  //for (auto j : *emcpflow) { emcpflowOR.push_back(j); }


  //std::remove_if(empflowOR.begin(), empflowOR.end(), jet_overlaps_photon);
  std::remove_copy_if(empflow->begin(), empflow->end(), std::back_inserter(empflowOR), jet_overlaps_photon);
  //std::remove_if(emcpflowOR.begin(), emcpflowOR.end(), jet_overlaps_photon);
  std::remove_copy_if(emcpflow->begin(), emcpflow->end(), std::back_inserter(emcpflowOR), jet_overlaps_photon);

  out->add_jets("empf", empflowOR);
  out->add_jets("emcpf", emcpflowOR);

  cout << "emp: " << empflow->size() << "(" << empflowOR.size() << ") emcp: " << emcpflow->size() << "(" << emcpflowOR.size() << ")" << endl;

  cout << "Calling Fill()" << endl;
  out->Fill();
}

int main(int argc, char **argv) {
  int nlim = 500;
  
  std::string input_filename;
  std::string output_filename = "output.root";
  if (argc < 2) {
    std::cerr << "Please give a filename!" << endl;
    cout << "Usage: " << argv[0] << " input_file [output_file [max_events]]" << endl;
    return 1;
  } else {
    input_filename = argv[1];
  }
  if (argc > 2) {
    output_filename = argv[2];
  }
  if (argc > 3) {
    nlim = atoi(argv[3]);
  }

  cout << "Initializing xAOD..." << endl;
  xAOD::Init();

  cout << "Opening file..." << endl;
  TFile *input_file = TFile::Open(input_filename.c_str());

  cout << "Creating TEvent..." << endl;
  xAOD::TEvent *evt = new xAOD::TEvent(input_file);

  int nevt = evt->getEntries();

  TFile* out_file = new TFile(output_filename.c_str(), "recreate");
  OutputTree* out_tree = new OutputTree();

  for (int i = 0; i < nevt; ++i) {
    if (nlim>0 && i>nlim) break;
    cout << "======= Event " << i << " =======" << endl;

    out_tree->clear();
    evt->getEntry(i);
    process_event(evt, out_tree);

    cout << "========================" << endl << endl;
  }

  out_tree->Write();
  //out_file->Write();
  out_file->Close();

  return 0;
}
