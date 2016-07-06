#include "trackmass/OutputTree.h"

#include "TBranch.h"

OutputTree::OutputTree() :
  TTree("tracks", "tracks")
{
}

void OutputTree::add_photon_type(const std::string& name) {
  if (ph_names.count(name) != 0) return;

  TBranch* b_pt = Branch((name+"_pt").c_str(), &ph_vars[name+"_pt"]);
  TBranch* b_eta = Branch((name+"_eta").c_str(), &ph_vars[name+"_eta"]);
  TBranch* b_phi = Branch((name+"_phi").c_str(), &ph_vars[name+"_phi"]);

  int n_entries = GetEntries();
  for (int i = 0; i < n_entries; ++i) {
    b_pt->Fill();
    b_eta->Fill();
    b_phi->Fill();
  }

  ph_names.insert(name);
}

void OutputTree::add_photon(const std::string& name, const TLorentzVector& ph) {
  add_photon_type(name);

  ph_vars[name+"_pt"].push_back(ph.Pt());
  ph_vars[name+"_eta"].push_back(ph.Eta());
  ph_vars[name+"_phi"].push_back(ph.Phi());
}

void OutputTree::add_jet_type(const std::string& name) {
  if (jet_names.count(name) != 0) return;

  TBranch* b_pt = Branch((name+"_pt").c_str(), &jet_vars[name+"_pt"]);
  TBranch* b_eta = Branch((name+"_eta").c_str(), &jet_vars[name+"_eta"]);
  TBranch* b_phi = Branch((name+"_phi").c_str(), &jet_vars[name+"_phi"]);
  TBranch* b_m = Branch((name+"_m").c_str(), &jet_vars[name+"_m"]);

  int n_entries = GetEntries();
  for (int i = 0; i < n_entries; ++i) {
    b_pt->Fill();
    b_eta->Fill();
    b_phi->Fill();
    b_m->Fill();
  }

  jet_names.insert(name);
}

void OutputTree::add_jet(const std::string& name, const TLorentzVector& j) {
  add_jet_type(name);

  jet_vars[name+"_pt"].push_back(j.Pt());
  jet_vars[name+"_eta"].push_back(j.Eta());
  jet_vars[name+"_phi"].push_back(j.Phi());
  jet_vars[name+"_m"].push_back(j.M());
}

void OutputTree::add_jet(const std::string& name, const fastjet::PseudoJet& j) {
  add_jet_type(name);

  jet_vars[name+"_pt"].push_back(j.pt());
  jet_vars[name+"_eta"].push_back(j.eta());
  jet_vars[name+"_phi"].push_back(j.phi());
  jet_vars[name+"_m"].push_back(j.m());
}

void OutputTree::add_jets(const std::string& name, const std::vector<TLorentzVector>& jets) {
  for (auto j : jets) {
    add_jet(name, j);
  }
}

void OutputTree::add_jets(const std::string& name, const std::vector<fastjet::PseudoJet>& jets) {
  for (auto j : jets) {
    add_jet(name, j);
  }
}

void OutputTree::add_jets(const std::string& name, const xAOD::JetContainer* jets) {
  for (auto j : *jets) {
    add_jet(name, j->p4());
  }
}

void OutputTree::add_jets(const std::string& name, const std::vector<const xAOD::Jet*>& jets) {
  for (auto j : jets) {
    add_jet(name, j->p4());
  }
}

void OutputTree::clear() {
  for (auto &p : ph_vars) {
    p.second.clear();
  }
  for (auto &p : jet_vars) {
    p.second.clear();
  }
}
