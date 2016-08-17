#include <iostream>

#include <string>
#include <stdexcept>
#include <algorithm>

#include "xAODBase/IParticle.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleContainer.h"

#include "OutputTree/OutputTree.h"
#include "trackmass/Retrieve.h"

#include "TFile.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

void process_event(xAOD::TEvent* evt1, xAOD::TEvent* evt3, OutputTree* out) {
  // Get the MC Event information
  const xAOD::EventInfo* info1(nullptr);
  const xAOD::EventInfo* info3(nullptr);

  evt1->retrieve(info1, "EventInfo");
  evt3->retrieve(info3, "EventInfo");

  // Check that the event numbers are consistent between files.
  unsigned int evt1_number = info1->eventNumber();
  unsigned int evt3_number = info3->eventNumber();
  if (evt1_number != evt3_number) {
    throw std::runtime_error("Event numbers do not match!");
  }

  // Get a list of partons from TRUTH1
  auto truth_particles = Retrieve<xAOD::TruthParticle>(evt1, "TruthParticles");

  std::vector<const xAOD::TruthParticle* > partons;
  std::vector<const xAOD::TruthParticle* > zprimes;
  for (auto p : *truth_particles) {
    // check for the intermediate Z' state
    if (p->status() == 22 && p->pdgId() == 101) {
      //cout << "Found Zprime   with status: " << p->status() << "  ID: " << p->pdgId() << "\tpT: " << p->pt() << "\tm: " << p->m() << endl;
      zprimes.push_back(p);
      continue;
    }

    // look at only final-state particles from the hard process
    if (p->status() != 23) continue;

    // keep only quarks and gluons
    if (abs(p->pdgId()) > 5 && p->pdgId() != 21) continue;

    //cout << "Found parton   with status: " << p->status() << "  ID: " << p->pdgId() << "\tpT: " << p->pt() << "\tm: " << p->m() << endl;

    partons.push_back(p);
  }

  // make a function that comapres particle pT
  auto pt_comp = [](const xAOD::IParticle* p1, const xAOD::IParticle* p2) -> bool {
      return p1->pt() > p2->pt();
  };

  // sort particles by their pT
  std::sort(partons.begin(), partons.end(), pt_comp);
  std::sort(zprimes.begin(), zprimes.end(), pt_comp);
      
  // write the truth particles to the ntuple branches.
  out->add_truths("zp", zprimes);
  out->add_truths("parton", partons);


  // grab the jets into a mutable vector container
  auto jets = CopyRetrieve<xAOD::Jet>(evt1, "AntiKt4TruthJets");
  auto fatjets = CopyRetrieve<xAOD::Jet>(evt3, "TrimmedAntiKt10TruthJets");

  // sort them by pT
  std::sort(jets.begin(), jets.end(), pt_comp);
  std::sort(fatjets.begin(), fatjets.end(), pt_comp);

  // write to tree
  out->add_jets("jet", jets);
  out->add_jets("fjet", fatjets);

  std::vector<float> fj_dR;
  if (zprimes.size() > 0) {
    auto z4 = zprimes[0]->p4();
    for (auto j : fatjets) {
      fj_dR.push_back(z4.DeltaR(j->p4()));
    }
  }
  out->add_vector("fjet_dR", fj_dR);

  // commit this event to the output ntuple.
  out->Fill();
}

void usage(int /*argc*/, char **argv) {
  cout << "Usage: " << argv[0] << " truth1_file truth3_file output_file" << endl;
}

int main(int argc, char **argv) {
  string truth1_filename;
  string truth3_filename;
  string output_filename;

  if (argc < 4) {
    cerr << "Not enough arguments supplied!" << endl;
    usage(argc, argv);
    return 1;
  }
  
  truth1_filename = argv[1];
  truth3_filename = argv[2];
  output_filename = argv[3];

  cout << "Initializing xAOD..." << endl;
  xAOD::Init();

  cout << "Opening TRUTH1 file:" << endl;
  cout << truth1_filename << endl;
  TFile *truth1_file = TFile::Open(truth1_filename.c_str());

  cout << "Opening TRUTH3 file:" << endl;
  cout << truth3_filename << endl;
  TFile *truth3_file = TFile::Open(truth3_filename.c_str());

  cout << "Creating TEvents..." << endl;
  xAOD::TEvent *evt1 = new xAOD::TEvent(truth1_file);
  xAOD::TEvent *evt3 = new xAOD::TEvent(truth3_file);

  int nevt_truth1 = evt1->getEntries();
  int nevt_truth3 = evt3->getEntries();

  cout << "There are " << nevt_truth1 << " : " << nevt_truth3 << " events in the files." << endl;

  if (nevt_truth1 != nevt_truth3) {
    cerr << "Error! TRUTH1 and TRUTH3 have different number of events. Abort." << endl;
    return 1;
  }

  TFile* out_file = new TFile(output_filename.c_str(), "recreate");
  OutputTree* out_tree = new OutputTree("truth");

  for (int i = 0; i < nevt_truth1; ++i) {
    //if (i > 100) break;
    
    if (i%500==0) {
      cout << "Processing event " << i << endl;
    }

    out_tree->clear();
    evt1->getEntry(i);
    evt3->getEntry(i);
    process_event(evt1, evt3, out_tree);

  }

  out_tree->Write();
  out_file->Close();

  return 0;
}
