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

#include "include/jetCollection.hh"
#include "include/softDropGroomer.hh"
#include "include/softDropCounter.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/csSubtractor.hh"
#include "include/csSubtractorFullEvent.hh"
#include "include/csSubFullEventIterative.hh"
#include "include/Angularity.hh"
#include "include/dyGroomer.hh"

using namespace std;
using namespace fastjet;

// This class runs time and C/A reclustering for jewel
// histo->Fill(var, weight); (root)
// ./runTimeClusJewel -hard /data/alice/jlomker/JEWEL/eventfiles/2.2.0.hepmc -hardtype HepMC2 -nev 100

int main (int argc, char ** argv) {

  auto start_time = std::chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  std::cout << "will run on " << nEvent << " events" << std::endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  TFile *fout = new TFile("JetToyHIResultTimeClusJewel.root","RECREATE");
  treeWriter trwSig("jetTreeSig");
 
  //Jet definition
  double R                   = 0.2;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  fastjet::JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 3;//3.0;
  fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  // Adding width and pTD from lambda variable for cross checks
  Angularity width(1.,1.,R);
  Angularity pTD(0.,2.,R);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;//what is this doing ?
  while ( mixer.next_event() && iev < nEvent )
  {
    // increment event number    
    iev++;
    // std::cout << "begin " << std::endl;

    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    std::vector<fastjet::PseudoJet> particlesMergedAll = mixer.particles();

    std::vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    // extract thermal and dummies (parton -> dummy)
    fastjet::Selector dummy_selector = SelectorVertexNumber(-1);
    vector<PseudoJet> particlesDummy = dummy_selector(particlesMergedAll);

    // select final state particles from hard event only
    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    // select final state particles from background event only
    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );

    //charged particles
    fastjet::Selector charged_selector = SelectorIsCharged();
    vector<PseudoJet> particlesSigCh = charged_selector(particlesSig);

    //std::cout << "#particles: " << particlesSig.size() << " of which charged: " << particlesSigCh.size() << std::endl;

    //remove ghosts from jewel dummies
    for(int i = 0; i < (int)particlesDummy.size(); i++){
      if(particlesDummy[i].perp() < 1e-5 && fabs(particlesDummy[i].pz()) > 2000){ // where r these values from ?
        particlesDummy.erase(particlesDummy.begin() + i);
        i = i - 1;
      }
    }

    
    //---------------------------------------------------------------------------
    //   subtract medium response in full event -- This or csSubtractor ?
    //---------------------------------------------------------------------------
   /*
    fastjet::contrib::ConstituentSubtractor subtractor;
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR); // distance in eta-phi plane
    subtractor.set_max_distance(0.5); // free parameter for the maximal allowed particle i and ghost k
    subtractor.set_alpha(0.); // free parameter for the distance measure (the exponent of particle pt). Note that in older versions of the packe alpha was multiplied by two bit in newe versions this is not the case anymore
    //subtractor.set_scale_fourmomentum(); // what is this, why is this not in use ?
    subtractor.set_remove_all_zero_pt_particles(true);

    std::vector<fastjet::PseudoJet> subtracted_particles = subtractor.do_subtraction(particlesSig, particlesDummy);


	*/
    //-----------------------------------------------------------------------------
    //   look at first splitting of hard parton -- correct/add if 1st split info available 
    //-----------------------------------------------------------------------------


    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // std::cout << "jet clustering" << std::endl;
    // run the clustering, extract the CS signal jets
    //fastjet::ClusterSequenceArea csSig(subtracted_particles, jet_def, area_def);//switch to subtracted particles
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(15.))));

    // run the clustering, extract the signal charged jets
    fastjet::ClusterSequenceArea csSigCh(particlesSigCh, jet_def, area_def);
    jetCollection jetCollectionSigCh(sorted_by_pt(jet_selector(csSigCh.inclusive_jets(15.))));

    //---------------------------------------------------------------------------
    //   Recursive Soft Drop for signal jets
    //---------------------------------------------------------------------------
    //  std::cout << "Reclustering" << std::endl;
    softDropCounter sdcSig(0.0,0.0,R,0.0);
    sdcSig.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSig.run(jetCollectionSig);

    softDropCounter sdcTau(0.0,0.0,R,0.0);
    sdcTau.setRecursiveAlgo(3);// tau algo.
    sdcTau.run(jetCollectionSig);


    jetCollectionSig.addVector("sigJetRecur_jetpt",     sdcSig.getPts());
    jetCollectionSig.addVector("sigJetRecur_jeteta",    sdcSig.getEtas());
    jetCollectionSig.addVector("sigJetRecur_kt",        sdcSig.getKts());
    jetCollectionSig.addVector("sigJetRecur_z",         sdcSig.getZgs());
    jetCollectionSig.addVector("sigJetRecur_dr12",      sdcSig.getDRs());
    jetCollectionSig.addVector("sigJetRecur_erad",      sdcSig.getErads());
    jetCollectionSig.addVector("sigJetRecur_logdr12",   sdcSig.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecur_logztheta", sdcSig.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecur_nSD",       sdcSig.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecur_zSD",       sdcSig.calculateNSD(1.0));
    jetCollectionSig.addVector("sigJetRecur_tf",        sdcSig.getTfs());//johanna cos()
    jetCollectionSig.addVector("sigJetRecur_tfe",       sdcSig.getTfes());//morgan (pt's)
    jetCollectionSig.addVector("sigJetRecur_tfe2",      sdcSig.getTfes2());//1->3 split

    jetCollectionSig.addVector("sigJetRecurTau_jetpt",     sdcTau.getPts());
    jetCollectionSig.addVector("sigJetRecurTau_jeteta",    sdcTau.getEtas());
    jetCollectionSig.addVector("sigJetRecurTau_kt",        sdcTau.getKts());
    jetCollectionSig.addVector("sigJetRecurTau_z",         sdcTau.getZgs());
    jetCollectionSig.addVector("sigJetRecurTau_dr12",      sdcTau.getDRs());
    jetCollectionSig.addVector("sigJetRecurTau_erad",      sdcTau.getErads());
    jetCollectionSig.addVector("sigJetRecurTau_logdr12",   sdcTau.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecurTau_logztheta", sdcTau.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecurTau_nSD",       sdcTau.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecurTau_zSD",       sdcTau.calculateNSD(1.0));
    jetCollectionSig.addVector("sigJetRecurTau_tf",        sdcTau.getTfs());//high E & soft & collinear
    jetCollectionSig.addVector("sigJetRecurTau_tfe",       sdcTau.getTfes());//high E (1-cos)
    jetCollectionSig.addVector("sigJetRecurTau_tfe2",      sdcTau.getTfes2());//1->3 split

    // calculate some angularities
    vector<double> widthSig; widthSig.reserve(jetCollectionSig.getJet().size());
    vector<double> pTDSig;   pTDSig.reserve(jetCollectionSig.getJet().size());
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      widthSig.push_back(width.result(jet));
      pTDSig.push_back(pTD.result(jet));
    }
    jetCollectionSig.addVector("widthSig", widthSig);
    jetCollectionSig.addVector("pTDSig", pTDSig);

    //find closest parton for each jet
    //and matching
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
    jetCollectionSig.addVector("sigJetRecurZcut_erad",      sdcSigzcut.getErads());
    jetCollectionSig.addVector("sigJetRecurZcut_logdr12",   sdcSigzcut.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecurZcut_logztheta", sdcSigzcut.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecurZcut_tf",        sdcSigzcut.getTfs());//high E & soft & collinear
    jetCollectionSig.addVector("sigJetRecurZcut_tfe",       sdcSigzcut.getTfes());//high E (1-cos)
    jetCollectionSig.addVector("sigJetRecurZcut_tfe2",      sdcSigzcut.getTfes2());//1->3 split
    jetCollectionSig.addVector("sigJetRecurZcut_nSD",       sdcSigzcut.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecurZcut_zSD",       sdcSigzcut.calculateNSD(1.0));
    jetCollectionSig.addVector("sigJetRecurZcut_droppedBFS",sdcSigzcut.getDBFSs());
    jetCollectionSig.addVector("sigJetRecurZcut_droppedTfBFS",sdcSigzcut.getDTfBFSs());
    jetCollectionSig.addVector("sigJetRecurZcut_droppedTfeBFS",sdcSigzcut.getDTfeBFSs());
    jetCollectionSig.addVector("sigJetRecurZcut_droppedPts",sdcSigzcut.getPtBFSs());
    jetCollectionSig.addVector("sigJetRecurZcut_droppedKts",sdcSigzcut.getKtBFSs());
    jetCollectionSig.addVector("sigJetRecurZcut_droppedLog1drBFS",sdcSigzcut.getLog1DrBFSs());
    jetCollectionSig.addVector("sigJetRecurZcut_tau21",     sdcSigzcut.getTau21s());
    jetCollectionSig.addVector("sigJetRecurZcut_tau32",     sdcSigzcut.getTau32s());


    jetCollectionSig.addVector("sigJetRecurTauZcut_jetpt",     sdcSigTauzcut.getPts());
    jetCollectionSig.addVector("sigJetRecurTauZcut_jeteta",    sdcSigTauzcut.getEtas());
    jetCollectionSig.addVector("sigJetRecurTauZcut_kt",        sdcSigTauzcut.getKts());
    jetCollectionSig.addVector("sigJetRecurTauZcut_z",         sdcSigTauzcut.getZgs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_dr12",      sdcSigTauzcut.getDRs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_erad",      sdcSigTauzcut.getErads());
    jetCollectionSig.addVector("sigJetRecurTauZcut_logdr12",   sdcSigTauzcut.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_logztheta", sdcSigTauzcut.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_tf",        sdcSigTauzcut.getTfs());//high E & soft & collinear
    jetCollectionSig.addVector("sigJetRecurTauZcut_tfe",       sdcSigTauzcut.getTfes());// high E (1-cos)
    jetCollectionSig.addVector("sigJetRecurTauZcut_tfe2",      sdcSigTauzcut.getTfes2());//1->3 split
    jetCollectionSig.addVector("sigJetRecurTauZcut_nSD",       sdcSigTauzcut.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecurTauZcut_zSD",       sdcSigTauzcut.calculateNSD(1.0));
    jetCollectionSig.addVector("sigJetRecurTauZcut_droppedBFS",sdcSigTauzcut.getDBFSs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_droppedTfBFS",sdcSigTauzcut.getDTfBFSs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_droppedTfeBFS",sdcSigTauzcut.getDTfeBFSs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_droppedPts",sdcSigTauzcut.getPtBFSs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_droppedKts",sdcSigTauzcut.getKtBFSs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_droppedLog1drBFS",sdcSigTauzcut.getLog1DrBFSs());

    jetCollectionSig.addVector("sigJetRecurTauZcut_tau21",     sdcSigTauzcut.getTau21s());
    jetCollectionSig.addVector("sigJetRecurTauZcut_tau32",     sdcSigTauzcut.getTau32s());

    /*
    // match Ch jets to signal jets
    jetMatcher jmCh(R);
    jmCh.setBaseJets(jetCollectionSigCh);
    jmCh.setTagJets(jetCollectionSig);
    jmCh.matchJets();
    jmCh.reorderedToTag(jetCollectionSigCh);
    */
    //---------------------------------------------------------------------------
    //   Full event constituent subtraction -- this one
    //---------------------------------------------------------------------------

    csSubFullEventIterative csSubEmbedded( {0.0} , {0.1}, 0.005,ghostRapMax); // alpha, rParam, ghA, ghRapMax
    csSubEmbedded.setInputParticles(particlesMerged);
    csSubEmbedded.setMaxEta(1.0);
    fastjet::ClusterSequenceArea csEmbedded(csSubEmbedded.Subtract(), jet_def, area_def);
    jetCollection jetCollectionCSFull(sorted_by_pt(jet_selector(csEmbedded.inclusive_jets(15.))));

    std::vector<double> rho;
    std::vector<double> rhom;
    rho.push_back(csSubEmbedded.getRho());
    rhom.push_back(csSubEmbedded.getRhoM());

    //match CSFull jets to signal jets
    jetMatcher jmCSFull(R);
    jmCSFull.setBaseJets(jetCollectionCSFull);
    jmCSFull.setTagJets(jetCollectionSig);
    jmCSFull.matchJets();

    jmCSFull.reorderedToTag(jetCollectionCSFull);


    //---------------------------------------------------------------------------
    //  recluster subtracted (full event CS) jets
    //---------------------------------------------------------------------------
    softDropCounter sdcCSFull(0.0,0.0,R,0.0);
    sdcCSFull.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcCSFull.run(jetCollectionCSFull);

    jetCollectionCSFull.addVector("csFullJetRecur_jetpt",     sdcCSFull.getPts());
    jetCollectionCSFull.addVector("csFullJetRecur_jeteta",    sdcCSFull.getEtas());
    jetCollectionCSFull.addVector("csFullJetRecur_kt",        sdcCSFull.getKts());
    jetCollectionCSFull.addVector("csFullJetRecur_z",         sdcCSFull.getZgs());
    jetCollectionCSFull.addVector("csFullJetRecur_dr12",      sdcCSFull.getDRs());
    jetCollectionCSFull.addVector("csFullJetRecur_erad",      sdcCSFull.getErads());
    jetCollectionCSFull.addVector("csFullJetRecur_logdr12",   sdcCSFull.getLog1DRs());
    jetCollectionCSFull.addVector("csFullJetRecur_logztheta", sdcCSFull.getLogzDRs());
    jetCollectionCSFull.addVector("csFullJetRecur_tf",        sdcCSFull.getTfs());
    jetCollectionCSFull.addVector("csFullJetRecur_tfe",       sdcCSFull.getTfes());
    jetCollectionCSFull.addVector("csFullJetRecur_tfe2",      sdcCSFull.getTfes2());//highE 1->3 split
    jetCollectionCSFull.addVector("csFullJetRecur_nSD",       sdcCSFull.calculateNSD(0.0));
    jetCollectionCSFull.addVector("csFullJetRecur_zSD",       sdcCSFull.calculateNSD(1.0));

    softDropCounter sdcCSFullTau(0.0,0.0,R,0.0);
    sdcCSFullTau.setRecursiveAlgo(3);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcCSFullTau.run(jetCollectionCSFull);

    jetCollectionCSFull.addVector("csFullJetRecurTau_jetpt",     sdcCSFullTau.getPts());
    jetCollectionCSFull.addVector("csFullJetRecurTau_jeteta",    sdcCSFullTau.getEtas());
    jetCollectionCSFull.addVector("csFullJetRecurTau_kt",        sdcCSFullTau.getKts());
    jetCollectionCSFull.addVector("csFullJetRecurTau_z",         sdcCSFullTau.getZgs());
    jetCollectionCSFull.addVector("csFullJetRecurTau_dr12",      sdcCSFullTau.getDRs());
    jetCollectionCSFull.addVector("csFullJetRecurTau_erad",      sdcCSFullTau.getErads());
    jetCollectionCSFull.addVector("csFullJetRecurTau_logdr12",   sdcCSFullTau.getLog1DRs());
    jetCollectionCSFull.addVector("csFullJetRecurTau_logztheta", sdcCSFullTau.getLogzDRs());
    jetCollectionCSFull.addVector("csFullJetRecurTau_tf",        sdcCSFullTau.getTfs());
    jetCollectionCSFull.addVector("csFullJetRecurTau_tfe",       sdcCSFullTau.getTfes());
    jetCollectionCSFull.addVector("csFullJetRecurTau_tfe2",      sdcCSFullTau.getTfes2());//highE 1->3 split
    jetCollectionCSFull.addVector("csFullJetRecurTau_nSD",       sdcCSFullTau.calculateNSD(0.0));
    jetCollectionCSFull.addVector("csFullJetRecurTau_zSD",       sdcCSFullTau.calculateNSD(1.0));


    softDropCounter sdcCSFullzcut(0.1,0.0,R,0.0);
    sdcCSFullzcut.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcCSFullzcut.run(jetCollectionCSFull);


    jetCollectionCSFull.addVector("csFullJetRecurZcut_jetpt",     sdcCSFullzcut.getPts());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_jeteta",    sdcCSFullzcut.getEtas());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_kt",        sdcCSFullzcut.getKts());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_z",         sdcCSFullzcut.getZgs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_dr12",      sdcCSFullzcut.getDRs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_erad",      sdcCSFullzcut.getErads());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_logdr12",   sdcCSFullzcut.getLog1DRs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_logztheta", sdcCSFullzcut.getLogzDRs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_tf",        sdcCSFullzcut.getTfs());//high E & soft & colliner
    jetCollectionCSFull.addVector("csFullJetRecurZcut_tfe",       sdcCSFullzcut.getTfes());//high E (1-cos)
    jetCollectionCSFull.addVector("csFullJetRecurZcut_tfe2",      sdcCSFullzcut.getTfes2());//highE 1->3 split
    jetCollectionCSFull.addVector("csFullJetRecurZcut_nSD",       sdcCSFullzcut.calculateNSD(0.0));
    jetCollectionCSFull.addVector("csFullJetRecurZcut_zSD",       sdcCSFullzcut.calculateNSD(1.0));
    jetCollectionCSFull.addVector("csFullJetRecurZcut_droppedBFS",sdcCSFullzcut.getDBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_droppedTfBFS",sdcCSFullzcut.getDTfBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_droppedTfeBFS",sdcCSFullzcut.getDTfeBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_droppedPts",sdcCSFullzcut.getPtBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_droppedKts",sdcCSFullzcut.getKtBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_droppedLog1drBFS",sdcCSFullzcut.getLog1DrBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_tau21",     sdcCSFullzcut.getTau21s());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_tau32",     sdcCSFullzcut.getTau32s());

    // add the tau here
    softDropCounter sdcCSFullTauzcut(0.1,0.0,R,0.0);
    sdcCSFullTauzcut.setRecursiveAlgo(3);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcCSFullTauzcut.run(jetCollectionCSFull);

    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_jetpt",     sdcCSFullTauzcut.getPts());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_jeteta",    sdcCSFullTauzcut.getEtas());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_kt",        sdcCSFullTauzcut.getKts());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_z",         sdcCSFullTauzcut.getZgs());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_dr12",      sdcCSFullTauzcut.getDRs());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_erad",      sdcCSFullTauzcut.getErads());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_logdr12",   sdcCSFullTauzcut.getLog1DRs());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_logztheta", sdcCSFullTauzcut.getLogzDRs());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_tf",        sdcCSFullTauzcut.getTfs());//highE & soft & collinear
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_tfe",       sdcCSFullTauzcut.getTfes());//high E (1-cos)
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_tfe2",      sdcCSFullTauzcut.getTfes2());//highE 1->3 split
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_nSD",       sdcCSFullTauzcut.calculateNSD(0.0));
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_zSD",       sdcCSFullTauzcut.calculateNSD(1.0));
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_droppedBFS",sdcCSFullTauzcut.getDBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_droppedTfBFS",sdcCSFullTauzcut.getDTfBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_droppedTfeBFS",sdcCSFullTauzcut.getDTfeBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_droppedPts",sdcCSFullTauzcut.getPtBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_droppedKts",sdcCSFullTauzcut.getKtBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurTauZcut_droppedLog1drBFS",sdcCSFullTauzcut.getLog1DrBFSs());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_tau21",     sdcCSFullzcut.getTau21s());
    jetCollectionCSFull.addVector("csFullJetRecurZcut_tau32",     sdcCSFullzcut.getTau32s());


    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trwSig.addCollection("eventWeight",   eventWeight);

    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.addCollection("sigJetCh",      jetCollectionSigCh);

    trwSig.addCollection("csFullJet",         jetCollectionCSFull);

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
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
