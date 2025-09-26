#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "PU14/EventMixer.hh"
#include "PU14/CmdLine.hh"
#include "PU14/PU14.hh"
#include "PU14/EventSource.hh"

#include "include/jetCollection.hh"
#include "include/softDropGroomer.hh"
#include "include/softDropCounter.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/csSubtractor.hh"
#include "include/csSubtractorFullEvent.hh"
#include "include/csSubFullEventIterative.hh"
#include "include/Angularity.hh"
#include "include/AliceFastSim.hh"
#include "include/thermalAlice.hh"

using namespace std;
using namespace fastjet;

// This class runs time reclustering with background
//./runTimeClus_MC -hard input/vectorTree_LHC25a2b.root -hardtype ROOT -hardtreename o2Tracks -hardvarname particle_truth -reco input/vectorTree_LHC25a2b.root -recotype ROOT -recotreename o2Tracks -recovarname particle_reco -nev 100

int main (int argc, char ** argv) {

  auto start_time = std::chrono::steady_clock::now();
  AliceFastSim fastSim = AliceFastSim();//Bas

  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  std::cout << "will run on " << nEvent << " events" << std::endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  TFile *fout = new TFile("JetToyHIResultTimeClus_MC.root","RECREATE");
  
  treeWriter trwSig("jetTreeSig");
 
  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  fastjet::JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 100.;//3.0;
  fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( mixer.next_event() && iev < nEvent )
  {  
    if(mixer.particles().size()==0){
      std::cout<<"no event info"<<std::endl;
      continue;
    }
    // increment event number    
    iev++;
    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    std::cout<<"iev: "<<iev<<std::endl;
    std::vector<fastjet::PseudoJet> particlesMergedAll = mixer.particles();
    std::vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
   // eventWeight.push_back(mixer.pu_weight());

    std::cout<<"MergedAll.size(): "<<particlesMergedAll.size()<<std::endl;

    // run as --hard
    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesTruth = sig_selector(particlesMergedAll);// all input particles (PYTHIA/MC)

    std::cout<<"particlesTruth.size(): "<<particlesTruth.size()<<std::endl;
    // run as --reco
    fastjet::Selector detector_selector = SelectorVertexNumber(99);
    vector<PseudoJet> particlesDetector = detector_selector(particlesMergedAll);

    std::cout<<"particlesDetector.size(): "<<particlesDetector.size()<<std::endl;
    //---------------------------------------------------------------------------
    //   fastsim
    //---------------------------------------------------------------------------
    //fastSim.setInputEvent(particlesSig);
    //vector<PseudoJet> particlesTruth = fastSim.AliceAcceptance();// all particles from PYTHIA that fit into alice ranges (eta, pt, only charged)
    //vector<PseudoJet> particlesDetector = fastSim.AliceDetector();// particles smeared by detector 

    //std::cout << "#particles: " << particlesSig.size() << " of which charged: " << particlesSigCh.size() << std::endl;
 
    double hbarc = 0.19732697;
    double GeVtofm = 1./hbarc; //~5.068;
    int id = 0;

    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea truth(particlesTruth, jet_def, area_def);
    std::cout<<"after clustersequence"<<std::endl;
    jetCollection jetCollectionTruth(sorted_by_pt(jet_selector(truth.inclusive_jets(15))));//was 10

    fastjet::ClusterSequenceArea det(particlesDetector, jet_def, area_def);
    jetCollection jetCollectionDet(sorted_by_pt(jet_selector(det.inclusive_jets(15))));//was 10

    jetCollection jetCollectionDetectorMatched(jetCollectionDet);

    //match csSig jets to csDet
    jetMatcher jmDetector(0.6*R);// check with o2
    jmDetector.setBaseJets(jetCollectionTruth);
    jmDetector.setTagJets(jetCollectionDetectorMatched);
    jmDetector.matchJets();
    jmDetector.reorderedToBase(jetCollectionDetectorMatched);// reorder detector to sig jets 

    std::cout<<"after matching"<<std::endl;
    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for truth MC jets
    //---------------------------------------------------------------------------
    //std::cout << "Reclustering " << std::endl;
    softDropCounter sdcTruth(0.0,0.0,R,0.0);
    sdcTruth.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcTruth.run(jetCollectionTruth);

    softDropCounter sdcTauTruth(0.0,0.0,R,0.0);
    sdcTauTruth.setRecursiveAlgo(3);
    sdcTauTruth.run(jetCollectionTruth);

    jetCollectionTruth.addVector("truthJetRecur_jetpt",     sdcTruth.getPts());
    jetCollectionTruth.addVector("truthJetRecur_jeteta",    sdcTruth.getEtas());
    jetCollectionTruth.addVector("truthJetRecur_kt",        sdcTruth.getKts());
    jetCollectionTruth.addVector("truthJetRecur_z",         sdcTruth.getZgs());
    jetCollectionTruth.addVector("truthJetRecur_dr12",      sdcTruth.getDRs());
    jetCollectionTruth.addVector("truthJetRecur_tf",        sdcTruth.getTfs());//high E & soft & collinear
    jetCollectionTruth.addVector("truthJetRecur_tfe",       sdcTruth.getTfes());//high E (1-cos)
    jetCollectionTruth.addVector("truthJetRecur_tfe2",      sdcTruth.getTfes2());//1->3 split

    jetCollectionTruth.addVector("truthJetRecurTau_jetpt",     sdcTauTruth.getPts());
    jetCollectionTruth.addVector("truthJetRecurTau_jeteta",    sdcTauTruth.getEtas());
    jetCollectionTruth.addVector("truthJetRecurTau_kt",        sdcTauTruth.getKts());
    jetCollectionTruth.addVector("truthJetRecurTau_z",         sdcTauTruth.getZgs());
    jetCollectionTruth.addVector("truthJetRecurTau_dr12",      sdcTauTruth.getDRs());
    jetCollectionTruth.addVector("truthJetRecurTau_tf",        sdcTauTruth.getTfs());//high E & soft & collinear
    jetCollectionTruth.addVector("truthJetRecurTau_tfe",       sdcTauTruth.getTfes());//high E (1-cos)
    jetCollectionTruth.addVector("truthJetRecurTau_tfe2",      sdcTauTruth.getTfes2());//1->3 split
 
    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for truth MC jets - with zcut
    //---------------------------------------------------------------------------

    softDropCounter sdcTruthzcut(0.1,0.0,R,0.0);
    sdcTruthzcut.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcTruthzcut.run(jetCollectionTruth);

    softDropCounter sdcTauTruthzcut(0.1,0.0,R,0.0);
    sdcTauTruthzcut.setRecursiveAlgo(3);
    sdcTauTruthzcut.run(jetCollectionTruth);
    
    jetCollectionTruth.addVector("truthJetRecurZcut_jetpt",     sdcTruthzcut.getPts());
    jetCollectionTruth.addVector("truthJetRecurZcut_jeteta",    sdcTruthzcut.getEtas());
    jetCollectionTruth.addVector("truthJetRecurZcut_kt",        sdcTruthzcut.getKts());
    jetCollectionTruth.addVector("truthJetRecurZcut_z",         sdcTruthzcut.getZgs());
    jetCollectionTruth.addVector("truthJetRecurZcut_dr12",      sdcTruthzcut.getDRs());
    jetCollectionTruth.addVector("truthJetRecurZcut_tf",        sdcTruthzcut.getTfs());//high E & soft & collinear
    jetCollectionTruth.addVector("truthJetRecurZcut_tfe",       sdcTruthzcut.getTfes());//high E (1-cos)
    jetCollectionTruth.addVector("truthJetRecurZcut_tfe2",      sdcTruthzcut.getTfes2());//1->3 split

    jetCollectionTruth.addVector("truthJetRecurTauZcut_jetpt",     sdcTauTruthzcut.getPts());
    jetCollectionTruth.addVector("truthJetRecurTauZcut_jeteta",    sdcTauTruthzcut.getEtas());
    jetCollectionTruth.addVector("truthJetRecurTauZcut_kt",        sdcTauTruthzcut.getKts());
    jetCollectionTruth.addVector("truthJetRecurTauZcut_z",         sdcTauTruthzcut.getZgs());
    jetCollectionTruth.addVector("truthJetRecurTauZcut_dr12",      sdcTauTruthzcut.getDRs());
    jetCollectionTruth.addVector("truthJetRecurTauZcut_tf",        sdcTauTruthzcut.getTfs());//high E & soft & collinear
    jetCollectionTruth.addVector("truthJetRecurTauZcut_tfe",       sdcTauTruthzcut.getTfes());// high E (1-cos)
    jetCollectionTruth.addVector("truthJetRecurTauZcut_tfe2",      sdcTauTruthzcut.getTfes2());//1->3 split

    //----------------------------------------------------------------------------------
    //   Recursive Soft Drop for detector jets, matched and reordered according to truth
    //----------------------------------------------------------------------------------
    //std::cout << "Reclustering " << std::endl;
    softDropCounter sdcDet(0.0,0.0,R,0.0);
    sdcDet.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcDet.run(jetCollectionDetectorMatched);

    softDropCounter sdcTauDet(0.0,0.0,R,0.0);
    sdcTauDet.setRecursiveAlgo(3);
    sdcTauDet.run(jetCollectionDetectorMatched);

    jetCollectionDetectorMatched.addVector("detectorJetRecur_jetpt",     sdcDet.getPts());
    jetCollectionDetectorMatched.addVector("detectorJetRecur_jeteta",    sdcDet.getEtas());
    jetCollectionDetectorMatched.addVector("detectorJetRecur_kt",        sdcDet.getKts());
    jetCollectionDetectorMatched.addVector("detectorJetRecur_z",         sdcDet.getZgs());
    jetCollectionDetectorMatched.addVector("detectorJetRecur_dr12",      sdcDet.getDRs());
    jetCollectionDetectorMatched.addVector("detectorJetRecur_tf",        sdcDet.getTfs());//high E & soft & collinear
    jetCollectionDetectorMatched.addVector("detectorJetRecur_tfe",       sdcDet.getTfes());//high E (1-cos)
    jetCollectionDetectorMatched.addVector("detectorJetRecur_tfe2",      sdcDet.getTfes2());//1->3 split

    jetCollectionDetectorMatched.addVector("detectorJetRecurTau_jetpt",     sdcTauDet.getPts());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTau_jeteta",    sdcTauDet.getEtas());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTau_kt",        sdcTauDet.getKts());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTau_z",         sdcTauDet.getZgs());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTau_dr12",      sdcTauDet.getDRs());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTau_tf",        sdcTauDet.getTfs());//high E & soft & collinear
    jetCollectionDetectorMatched.addVector("detectorJetRecurTau_tfe",       sdcTauDet.getTfes());//high E (1-cos)
    jetCollectionDetectorMatched.addVector("detectorJetRecurTau_tfe2",      sdcTauDet.getTfes2());//1->3 split

    //----------------------------------------------------------------------------------------------
    //   Recursive Soft Drop for detector jets, matched and reordered according to truth - with zcut
    //----------------------------------------------------------------------------------------------

    softDropCounter sdcDetzcut(0.1,0.0,R,0.0);
    sdcDetzcut.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcDetzcut.run(jetCollectionDetectorMatched);

    softDropCounter sdcTauDetzcut(0.1,0.0,R,0.0);
    sdcTauDetzcut.setRecursiveAlgo(3);
    sdcTauDetzcut.run(jetCollectionDetectorMatched);

    jetCollectionDetectorMatched.addVector("detectorJetRecurZcut_jetpt",     sdcDetzcut.getPts());
    jetCollectionDetectorMatched.addVector("detectorJetRecurZcut_jeteta",    sdcDetzcut.getEtas());
    jetCollectionDetectorMatched.addVector("detectorJetRecurZcut_kt",        sdcDetzcut.getKts());
    jetCollectionDetectorMatched.addVector("detectorJetRecurZcut_z",         sdcDetzcut.getZgs());
    jetCollectionDetectorMatched.addVector("detectorJetRecurZcut_dr12",      sdcDetzcut.getDRs());
    jetCollectionDetectorMatched.addVector("detectorJetRecurZcut_tf",        sdcDetzcut.getTfs());//high E & soft & collinear
    jetCollectionDetectorMatched.addVector("detectorJetRecurZcut_tfe",       sdcDetzcut.getTfes());//high E (1-cos)
    jetCollectionDetectorMatched.addVector("detectorJetRecurZcut_tfe2",      sdcDetzcut.getTfes2());//1->3 split

    jetCollectionDetectorMatched.addVector("detectorJetRecurTauZcut_jetpt",     sdcTauDetzcut.getPts());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTauZcut_jeteta",    sdcTauDetzcut.getEtas());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTauZcut_kt",        sdcTauDetzcut.getKts());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTauZcut_z",         sdcTauDetzcut.getZgs());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTauZcut_dr12",      sdcTauDetzcut.getDRs());
    jetCollectionDetectorMatched.addVector("detectorJetRecurTauZcut_tf",        sdcTauDetzcut.getTfs());//high E & soft & collinear
    jetCollectionDetectorMatched.addVector("detectorJetRecurTauZcut_tfe",       sdcTauDetzcut.getTfes());// high E (1-cos)
    jetCollectionDetectorMatched.addVector("detectorJetRecurTauZcut_tfe2",      sdcTauDetzcut.getTfes2());//1->3 split


    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trwSig.addCollection("eventWeight",   eventWeight);
    trwSig.addCollection("truthJet",        jetCollectionTruth);//truth
    trwSig.addCollection("matchedDetJet", jetCollectionDetectorMatched);
    std::cout<<"after addingCOllections"<<std::endl;
    trwSig.fillTree();  //signal jets
    std::cout<<"after filling the tree"<<std::endl;
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  fout->cd();
  TTree *trOut = trwSig.getTree();
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;//was 1000
  //std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
