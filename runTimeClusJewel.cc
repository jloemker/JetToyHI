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
#include "include/Angularity.hh"
#include "include/dyGroomer.hh"

using namespace std;
using namespace fastjet;

// This class runs time and C/A reclustering for jewel
// histo->Fill(var, weight); (root)
// ./runTimeClusBkg -hard samples/PythiaEventsTune14PtHat120_10k.pu14 -pileup samples/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10

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
  treeWriter trwSig("jetTreeSig");
 
  //Jet definition
  double R                   = 0.4;
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
    //   subtract medium response in full event
    //---------------------------------------------------------------------------
    fastjet::contrib::ConstituentSubtractor subtractor;
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR); // distance in eta-phi plane
    subtractor.set_max_distance(0.5); // free parameter for the maximal allowed particle i and ghost k
    subtractor.set_alpha(0.); // free parameter for the distance measure (the exponent of particle pt). Note that in older versions of the packe alpha was multiplied by two bit in newe versions this is not the case anymore
    //subtractor.set_scale_fourmomentum(); // what is this, why is this not in use ?
    subtractor.set_remove_all_zero_pt_particles(true);

    std::vector<fastjet::PseudoJet> subtracted_particles = subtractor.do_subtraction(particlesSig, particlesDummy);

    //-----------------------------------------------------------------------------
    //   look at first splitting of hard parton 
    //-----------------------------------------------------------------------------

    /*
    std::vector<double> drsplit;
    std::vector<double> tfsplit;
    double hbarc = 0.19732697;
    double GeVtofm = 1./hbarc; //~5.068;
    int id = 0;
    for(int ip = 0; ip<partons.size(); ++ip) {
      std::cout << "1st split hard parton " << std::endl;

      PseudoJet p = partons[ip];
      PseudoJet d1 = partonsFirstSplit[id];
      PseudoJet d2 = partonsFirstSplit[id+1];
      double dr = d1.delta_R(d2);
      drsplit.push_back(dr);
      double z1 = max(d1.e(),d2.e())/p.e();
      double z2 = min(d1.e(),d2.e())/p.e();
      tfsplit.push_back(1./(2.*z1*z2*p.e()*GeVtofm*(1-fastjet::cos_theta(d1,d2))));
      std::cout << "end of calculation " << std::endl;
 
      id+=2;
    }
    */


    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // std::cout << "jet clustering" << std::endl;
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea csSig(subtracted_particles, jet_def, area_def);//switch to subtracted particles
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    // run the clustering, extract the signal charged jets
    fastjet::ClusterSequenceArea csSigCh(particlesSigCh, jet_def, area_def);
    jetCollection jetCollectionSigCh(sorted_by_pt(jet_selector(csSigCh.inclusive_jets(10.))));
    //jetCollection jetCollectionSigJewel(GetCorrectedJets(jetCollectionSig.getJet(), particlesDummy));

    // run clustering on subtracted event
    fastjet::ClusterSequenceArea csSigCS(subtracted_particles, jet_def, area_def);//our csSig is already there. 
    jetCollection jetCollectionSigCS(sorted_by_pt(jet_selector(csSigCS.inclusive_jets(10.))));

    //--------------------------------------------------------------------------
    //   dynamical Grooming
    //--------------------------------------------------------------------------
    dyGroomer dygESig(1);// i somehow set a here - which i beliee is the a from the kappa calculation ~= 1/tau_form(but the dtheta has exp. a instead of 2 now 1)
    dygESig.setRecursiveAlgo(0);// CA jetdef
    //jetCollection jetCollectionSigDYE(dygESig.doGrooming(jetCollectionSig));
    dygESig.doGrooming(jetCollectionSig);
    jetCollectionSig.addVector("sigJetkappaDYE",    dygESig.getKappas());
    jetCollectionSig.addVector("sigJetDYE_dr12",       dygESig.getDR12());
    jetCollectionSig.addVector("sigJetDYE_z",       dygESig.getZgs());
    jetCollectionSig.addVector("sigJetDYE_kt",       dygESig.getKts());
    jetCollectionSig.addVector("sigJetDYE_tf",       dygESig.getTfs());
    jetCollectionSig.addVector("sigJetDYE_tfe",       dygESig.getTfes());
    jetCollectionSig.addVector("sigJetDYE_tau21",       dygESig.getTau21());
    jetCollectionSig.addVector("sigJetDYE_tau32",       dygESig.getTau32());
    jetCollectionSig.addVector("sigJetDYE_drBranches",       dygESig.getNDroppedSubjets());

    dyGroomer dygESigTau(1);
    dygESigTau.setRecursiveAlgo(3);// Tau jetdef
    dygESigTau.doGrooming(jetCollectionSig);
    jetCollectionSig.addVector("sigJetTaukappaDYE",    dygESig.getKappas());
    jetCollectionSig.addVector("sigJetTauDYE_dr12",       dygESig.getDR12());
    jetCollectionSig.addVector("sigJetTauDYE_z",       dygESig.getZgs());
    jetCollectionSig.addVector("sigJetTauDYE_kt",       dygESig.getKts());
    jetCollectionSig.addVector("sigJetTauDYE_tf",       dygESig.getTfs());
    jetCollectionSig.addVector("sigJetTauDYE_tfe",       dygESig.getTfes());
    jetCollectionSig.addVector("sigJetTauDYE_tau21",       dygESig.getTau21());
    jetCollectionSig.addVector("sigJetTauDYE_tau32",       dygESig.getTau32());
    jetCollectionSig.addVector("sigJetTauDYE_drBranches",       dygESig.getNDroppedSubjets());

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
    jetCollectionSig.addVector("sigJetRecur_z",         sdcSig.getZgs());
    jetCollectionSig.addVector("sigJetRecur_dr12",      sdcSig.getDRs());
    jetCollectionSig.addVector("sigJetRecur_erad",      sdcSig.getErads());
    jetCollectionSig.addVector("sigJetRecur_logdr12",   sdcSig.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecur_logztheta", sdcSig.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecur_tf",        sdcSig.getTfs());
    jetCollectionSig.addVector("sigJetRecur_tfe",       sdcSig.getTfes());
    jetCollectionSig.addVector("sigJetRecur_nSD",       sdcSig.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecur_zSD",       sdcSig.calculateNSD(1.0));

    jetCollectionSig.addVector("sigJetRecurTau_jetpt",     sdcTau.getPts());
    jetCollectionSig.addVector("sigJetRecurTau_z",         sdcTau.getZgs());
    jetCollectionSig.addVector("sigJetRecurTau_dr12",      sdcTau.getDRs());
    jetCollectionSig.addVector("sigJetRecurTau_erad",      sdcTau.getErads());
    jetCollectionSig.addVector("sigJetRecurTau_logdr12",   sdcTau.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecurTau_logztheta", sdcTau.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecurTau_tf",        sdcTau.getTfs());
    jetCollectionSig.addVector("sigJetRecurTau_tfe",       sdcTau.getTfes());
    jetCollectionSig.addVector("sigJetRecurTau_nSD",       sdcTau.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecurTau_zSD",       sdcTau.calculateNSD(1.0));

    // calculate some angularities
    vector<double> widthSig; widthSig.reserve(jetCollectionSig.getJet().size());
    vector<double> pTDSig;   pTDSig.reserve(jetCollectionSig.getJet().size());
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      widthSig.push_back(width.result(jet));
      pTDSig.push_back(pTD.result(jet));
    }
    jetCollectionSig.addVector("widthSig", widthSig);
    jetCollectionSig.addVector("pTDSig", pTDSig);

    /*
    //find closest parton for each jet
    std::vector<int> partonmatch;
    std::vector<double> partonmatchdr;
    std::vector<fastjet::PseudoJet> sigJets =  jetCollectionSig.getJet();
    for(fastjet::PseudoJet p : sigJets) {
      int ipmin = -1;
      double drmin = 999.;
      for(int ip = 0; ip<partons.size(); ++ip) {
        double dr = p.delta_R(partons[ip]);
        if(dr<drmin) {
          drmin = dr;
          ipmin = ip;
        }
      }
      partonmatch.push_back(ipmin);
      partonmatchdr.push_back(drmin);
    }
    jetCollectionSig.addVector("sigJetRecur_partonMatchID", partonmatch);
    jetCollectionSig.addVector("sigJetRecur_partonMatchDr", partonmatchdr);
    */
    // std::cout << "SD " << std::endl;

    softDropCounter sdcSigzcut(0.1,0.0,R,0.0);
    sdcSigzcut.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSigzcut.run(jetCollectionSig);

    softDropCounter sdcSigTauzcut(0.1,0.0,R,0.0);
    sdcSigTauzcut.setRecursiveAlgo(3);
    sdcSigTauzcut.run(jetCollectionSig);
    
    jetCollectionSig.addVector("sigJetRecurZcut_jetpt",     sdcSigzcut.getPts());
    jetCollectionSig.addVector("sigJetRecurZcut_z",         sdcSigzcut.getZgs());
    jetCollectionSig.addVector("sigJetRecurZcut_dr12",      sdcSigzcut.getDRs());
    jetCollectionSig.addVector("sigJetRecurZcut_erad",      sdcSigzcut.getErads());
    jetCollectionSig.addVector("sigJetRecurZcut_logdr12",   sdcSigzcut.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecurZcut_logztheta", sdcSigzcut.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecurZcut_tf",        sdcSigzcut.getTfs());
    jetCollectionSig.addVector("sigJetRecurZcut_tfe",       sdcSigzcut.getTfes());
    jetCollectionSig.addVector("sigJetRecurZcut_nSD",       sdcSigzcut.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecurZcut_zSD",       sdcSigzcut.calculateNSD(1.0));
    jetCollectionSig.addVector("sigJetRecurZcut_tau21",     sdcSigzcut.getTau21s());
    jetCollectionSig.addVector("sigJetRecurZcut_tau32",     sdcSigzcut.getTau32s());
    
    jetCollectionSig.addVector("sigJetRecurTauZcut_jetpt",     sdcSigTauzcut.getPts());
    jetCollectionSig.addVector("sigJetRecurTauZcut_z",         sdcSigTauzcut.getZgs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_dr12",      sdcSigTauzcut.getDRs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_erad",      sdcSigTauzcut.getErads());
    jetCollectionSig.addVector("sigJetRecurTauZcut_logdr12",   sdcSigTauzcut.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_logztheta", sdcSigTauzcut.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_tf",        sdcSigTauzcut.getTfs());
    jetCollectionSig.addVector("sigJetRecurTauZcut_tfe",       sdcSigTauzcut.getTfes());
    jetCollectionSig.addVector("sigJetRecurTauZcut_nSD",       sdcSigTauzcut.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecurTauZcut_zSD",       sdcSigTauzcut.calculateNSD(1.0));
    jetCollectionSig.addVector("sigJetRecurTauZcut_tau21",     sdcSigTauzcut.getTau21s());
    jetCollectionSig.addVector("sigJetRecurTauZcut_tau32",     sdcSigTauzcut.getTau32s());

    //---------------------------------------------------------------------------
    //   jet clustering of charged-particle signal jets
    //---------------------------------------------------------------------------
    softDropCounter sdcSigzcutCh(0.1,0.0,R,0.0);
    sdcSigzcutCh.setRecursiveAlgo(0);//0 = CA 1 = AKT 2 = KT  3=gen_kt t-form ordered
    sdcSigzcutCh.run(jetCollectionSigCh);

    softDropCounter sdcSigTauzcutCh(0.1,0.0,R,0.0);//SD z = 0.1
    sdcSigTauzcutCh.setRecursiveAlgo(3);// tau 
    sdcSigTauzcutCh.run(jetCollectionSigCh);
    
    jetCollectionSigCh.addVector("sigJetChRecurZcut_jetpt",     sdcSigzcutCh.getPts());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_z",         sdcSigzcutCh.getZgs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_dr12",      sdcSigzcutCh.getDRs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_erad",      sdcSigzcutCh.getErads());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_logdr12",   sdcSigzcutCh.getLog1DRs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_logztheta", sdcSigzcutCh.getLogzDRs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_tf",        sdcSigzcutCh.getTfs());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_tfe",       sdcSigzcutCh.getTfes());
    jetCollectionSigCh.addVector("sigJetChRecurZcut_nSD",       sdcSigzcutCh.calculateNSD(0.0));
    jetCollectionSigCh.addVector("sigJetChRecurZcut_zSD",       sdcSigzcutCh.calculateNSD(1.0));
 
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_jetpt",     sdcSigTauzcutCh.getPts());               
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_z",         sdcSigTauzcutCh.getZgs());               
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_dr12",      sdcSigTauzcutCh.getDRs());               
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_erad",      sdcSigTauzcutCh.getErads());             
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_logdr12",   sdcSigTauzcutCh.getLog1DRs());           
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_logztheta", sdcSigTauzcutCh.getLogzDRs());           
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_tf",        sdcSigTauzcutCh.getTfs());               
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_tfe",       sdcSigTauzcutCh.getTfes());              
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_nSD",       sdcSigTauzcutCh.calculateNSD(0.0));      
    jetCollectionSigCh.addVector("sigJetRecurTauZcutCh_zSD",       sdcSigTauzcutCh.calculateNSD(1.0)); 


    /*
    //find closest parton for each charged jet
    std::vector<int> partonmatchCh;
    std::vector<double> partonmatchdrCh;
    std::vector<fastjet::PseudoJet> sigJetsCh =  jetCollectionSigCh.getJet();
    for(fastjet::PseudoJet p : sigJetsCh) {
      int ipmin = -1;
      double drmin = 999.;
      for(int ip = 0; ip<partons.size(); ++ip) {
        double dr = p.delta_R(partons[ip]);
        if(dr<drmin) {
          drmin = dr;
          ipmin = ip;
        }
      }
      partonmatchCh.push_back(ipmin);
      partonmatchdrCh.push_back(drmin);
    }
    jetCollectionSigCh.addVector("sigJetChRecur_partonMatchID", partonmatchCh);
    jetCollectionSigCh.addVector("sigJetChRecur_partonMatchDr", partonmatchdrCh);
    */

    // match CS jets to signal jets
    jetMatcher jmCS(R);
    jmCS.setBaseJets(jetCollectionSigCS);
    jmCS.setTagJets(jetCollectionSig);
    jmCS.matchJets();

    // some grooming of a jetcollection that i dont need

    jmCS.reorderedToTag(jetCollectionSigCS);

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trwSig.addCollection("eventWeight",   eventWeight);

    //trwSig.addPartonCollection("partons",       partons);
    //trwSig.addPartonCollection("partonsFirstSplit",       partonsFirstSplit);
    //trwSig.addDoubleCollection("drSplit", drsplit);
    //trwSig.addDoubleCollection("tfSplit", tfsplit);
    
    trwSig.addCollection("sigJet",        jetCollectionSig);
    trwSig.addCollection("sigJetCh",      jetCollectionSigCh);

    //trwSig.addCollection("sigJetCS",      jetCollectionSigCS);
    
    trwSig.fillTree();  //signal jets
  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TFile *fout = new TFile("JetToyHIResultTimeClusJewel.root","RECREATE");
  trwSig.getTree()->Write();
  
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;//was 1000
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
