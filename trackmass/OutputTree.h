#include "TTree.h"
#include "TLorentzVector.h"

#include "fastjet/PseudoJet.hh"

#include "xAODJet/JetContainer.h"

class OutputTree : public TTree {
public:
  OutputTree();
  
  void add_photon_type(const std::string& name);

  void add_photon(const std::string& name, const TLorentzVector& ph);

  void add_jet_type(const std::string& name);
  
  void add_jet(const std::string& name, const TLorentzVector& j);
  void add_jet(const std::string& name, const fastjet::PseudoJet& j);
  void add_jets(const std::string& name, const std::vector<TLorentzVector>& jets);
  void add_jets(const std::string& name, const std::vector<fastjet::PseudoJet>& jets);
  void add_jets(const std::string& name, const xAOD::JetContainer* jets);
  void add_jets(const std::string& name, const std::vector<const xAOD::Jet*>& jets);

  void clear();
  
private:
  std::map<std::string, std::vector<float> > jet_vars;
  std::set<std::string> jet_names;

  std::map<std::string, std::vector<float> > ph_vars;
  std::set<std::string> ph_names;
};

