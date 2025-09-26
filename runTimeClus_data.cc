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
//./runTimeClus_data -hard vectorTree.root -hardtype ROOT -hardtreename o2Tracks -hardvarname particle_data -nev 10

int main (int argc, char ** argv) {

  auto start_time = std::chrono::steady_clock::now();
  AliceFastSim fastSim = AliceFastSim();//Bas
  thermalAlice thrmEvent;//Bas

  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  std::cout << "will run on " << nEvent << " events" << std::endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  TFile *fout = new TFile("JetToyHIResultTimeClus_Data.root","RECREATE");
  
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

    std::vector<fastjet::PseudoJet> particlesMergedAll = mixer.particles();
    std::vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    // select final state particles from hard event only
    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    //vector<PseudoJet> particlesMerged = particlesBkg;//Marta
    vector<PseudoJet> particlesMerged = particlesSig;//Bas
    //particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );//Marta

    double hbarc = 0.19732697;
    double GeVtofm = 1./hbarc; //~5.068;
    int id = 0;

    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    //std::cout << "jet clustering" << std::endl;
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea sig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(sig.inclusive_jets(15.))));

    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for signal jets
    //---------------------------------------------------------------------------
    //std::cout << "Reclustering " << std::endl;
    softDropCounter sdcSig(0.0,0.0,R,0.0);
    sdcSig.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSig.run(jetCollectionSig);

    softDropCounter sdcTau(0.0,0.0,R,0.0);
    sdcTau.setRecursiveAlgo(3);
    sdcTau.run(jetCollectionSig);

    jetCollectionSig.addVector("sigJetRecur_jetpt",     sdcSig.getPts());
    jetCollectionSig.addVector("sigJetRecur_jeteta",    sdcSig.getEtas());
    jetCollectionSig.addVector("sigJetRecur_kt",        sdcSig.getKts());
    jetCollectionSig.addVector("sigJetRecur_z",         sdcSig.getZgs());
    jetCollectionSig.addVector("sigJetRecur_dr12",      sdcSig.getDRs());
    jetCollectionSig.addVector("sigJetRecur_tf",        sdcSig.getTfs());//high E & soft & collinear
    jetCollectionSig.addVector("sigJetRecur_tfe",       sdcSig.getTfes());//high E (1-cos)
    jetCollectionSig.addVector("sigJetRecur_tfe2",      sdcSig.getTfes2());//1->3 split

    jetCollectionSig.addVector("sigJetRecurTau_jetpt",     sdcTau.getPts());
    jetCollectionSig.addVector("sigJetRecurTau_jeteta",    sdcTau.getEtas());
    jetCollectionSig.addVector("sigJetRecurTau_kt",        sdcTau.getKts());
    jetCollectionSig.addVector("sigJetRecurTau_z",         sdcTau.getZgs());
    jetCollectionSig.addVector("sigJetRecurTau_dr12",      sdcTau.getDRs());
    jetCollectionSig.addVector("sigJetRecurTau_tf",        sdcTau.getTfs());//high E & soft & collinear
    jetCollectionSig.addVector("sigJetRecurTau_tfe",       sdcTau.getTfes());//high E (1-cos)
    jetCollectionSig.addVector("sigJetRecurTau_tfe2",      sdcTau.getTfes2());//1->3 split

    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for signal jets - with zcut
    //---------------------------------------------------------------------------

    softDropCounter sdcSigzcut(0.1,0.0,R,0.0);
    sdcSigzcut.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSigzcut.run(jetCollectionSig);

    softDropCounter sdcSigTauzcut(0.1,0.0,R,0.0);
    sdcSigTauzcut.setRecursiveAlgo(3);
    sdcSigTauzcut.run(jetCollectionSig);
    
    jetCollectionSig.addVector("sigJetRecurZcut_jetpt",     sdcSigzcut.getPts());
    jetCollectionSig.addVector("sigJetRecurZcut_jeteta",    sdcSigzcut.getEtas());
    jetCollectionSig.addVector("sigJetRecurZcut_kt",        sdcSigzcut.getKts());
    jetCollectionSig.addVector("sigJetRecurZcut_z",         sdcSigzcut.getZgs());
    jetCollectionSig.addVector("sigJetRecurZcut_dr12",      sdcSigzcut.getDRs());
    jetCollectionSig.addVector("sigJetRecurZcut_tf",        sdcSigzcut.getTfs());//high E & soft & collinear
    jetCollectionSig.addVector("sigJetRecurZcut_tfe",       sdcSigzcut.getTfes());//high E (1-cos)
    jetCollectionSig.addVector("sigJetRecurZcut_tfe2",      sdcSigzcut.getTfes2());//1->3 split
    jetCollectionSig.addVector("sigJetRecurZcut_nSD",       sdcSigzcut.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecurZcut_zSD",       sdcSigzcut.calculateNSD(1.0));

    jetCollectionSig.addVector("sigJetRecurTauZcut_jetpt",     sdcSigTauzcut.getPts());
    jetCollectionSig.addVector("sigJetRecurTauZcut_jeteta",    sdcSigTauzcut.getEtas());
    jetCollectionSig.addVector("sigJetRecurTauZcut_kt",        sdcSigTauzcut.getKts());
    jetCollectionSig.addVector("sigJetRecurTauZcut_z",         sdcSigTauzcut.getZgs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_dr12",      sdcSigTauzcut.getDRs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_tf",        sdcSigTauzcut.getTfs());//high E & soft & collinear
    jetCollectionSig.addVector("sigJetRecurTauZcut_tfe",       sdcSigTauzcut.getTfes());// high E (1-cos)
    jetCollectionSig.addVector("sigJetRecurTauZcut_tfe2",      sdcSigTauzcut.getTfes2());//1->3 split
   
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.fillTree();  //signal jets
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
