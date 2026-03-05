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

// ./runTimeClusBkg -hard samples/PythiaEventsTune14PtHat120_10k.pu14 -pileup samples/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10
bool checkNanValues(int id, vector<PseudoJet> partonsFirstSplit){
  if(partonsFirstSplit[id].m()<0.0000005 || partonsFirstSplit[id+1].m()<0.0000005 || partonsFirstSplit[id].pt()>500000000 || partonsFirstSplit[id+1].pt()>500000000 || partonsFirstSplit[id].pt()<0.05 || partonsFirstSplit[id+1].pt()<0.05 || !partonsFirstSplit[id].pt() || !partonsFirstSplit[id+1].pt() || !partonsFirstSplit[id].m() || !partonsFirstSplit[id+1].m() || partonsFirstSplit[id].rap() < 0.00001 || partonsFirstSplit[id+1].rap()<0.00001 || partonsFirstSplit[id].m() =='inf' || partonsFirstSplit[id+1].m() =='inf' || partonsFirstSplit[id].m() =='-nan' || partonsFirstSplit[id+1].m() =='-nan'){
   return false;
   }else{
   return true;
   }
}

int main (int argc, char ** argv) {

  auto start_time = std::chrono::steady_clock::now();
  AliceFastSim fastSim = AliceFastSim();//Bas
  thermalAlice thrmEvent;//Bas

  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  int bin = cmdline.value<int>("-bin",20)
  //bool verbose = cmdline.present("-verbose");

  std::cout << "will run on " << nEvent << " events" << std::endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  TFile *fout = new TFile("JetToyHIResultTimeClusBkg"+bin+".root","RECREATE");
  
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

  // Adding width and pTD from lambda variable for cross checks
  Angularity width(1.,1.,R);
  Angularity pTD(0.,2.,R);

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
    std::vector<fastjet::PseudoJet> particlesPileup = thrmEvent.createThermalEventAlice();//Bas
    thrmEvent.createThermalEventAlice();//Bas

    std::vector<fastjet::PseudoJet> particlesMergedAll = mixer.particles();
    std::vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    // extract hard partons that initiated the jets
    fastjet::Selector parton_selector = SelectorVertexNumber(-1);
    vector<PseudoJet> partons = parton_selector(particlesMergedAll);

    // extract hard partons from first splitting
    fastjet::Selector parton_selector_split = SelectorVertexNumber(-2);
    vector<PseudoJet> partonsFirstSplit = parton_selector_split(particlesMergedAll);

    // select final state particles from hard event only
    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    // select final state particles from background event only
    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    //vector<PseudoJet> particlesMerged = particlesBkg;//Marta
    vector<PseudoJet> particlesMerged = particlesSig;//Bas
    //particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );//Marta
    particlesMerged.insert( particlesMerged.end(), particlesPileup.begin(), particlesPileup.end() );
    //charged particles
    fastjet::Selector charged_selector = SelectorIsCharged();
    vector<PseudoJet> particlesSigCh = charged_selector(particlesSig);

    //std::cout << "#particles: " << particlesSig.size() << " of which charged: " << particlesSigCh.size() << std::endl;

 
    //---------------------------------------------------------------------------
    //   look at first splitting of hard partons
    //---------------------------------------------------------------------------
    std::vector<double> drsplit;
    std::vector<double> tfesplit;
    std::vector<double> zgsplit;
    std::vector<double> ktsplit;
    std::vector<double> partonWd1Wod2pt;
    std::vector<double> partonWod1Wd2pt;
    std::vector<double> partonWd1d2pt;
    std::vector<double> partonWod1d2pt;
    double hbarc = 0.19732697;
    double GeVtofm = 1./hbarc; //~5.068;
    int id = 0;

    //trwSig.addDoubleCollection("partonWoConst", partonWoConst);// to check for parton pt of those that don't split 
    //for(int ip = 0; ip<std::min(particlesMergedAll.size(),6); ++ip) {
    //std::cout<<"N event: "<<iev<<std::endl;
    int particlesSize = particlesMergedAll.size();
    for(int ip = 0; ip< std::min(particlesSize,6); ++ip) {
      bool d1True=false;
      bool d2True=false;
      PseudoJet p = particlesMergedAll[ip];
      PseudoJet d1;
      PseudoJet d2;
      int vtx = p.user_info<PU14>().vertex_number();
      if(vtx != -1) continue;
      
      if(vtx == -1) {
        if( (ip+1) < particlesMergedAll.size() ){
          d1 = particlesMergedAll[ip+1];
          if(d1.user_info<PU14>().vertex_number()==-2){ // std::cout << "yay daughter 1" << std::endl;
	    d1True=true;
	  }else{
	    std::cout<<" no 1st daughter, but enough entries "<<std::endl;
	  }
          if( (ip+2) < particlesMergedAll.size() ){
	    d2 = particlesMergedAll[ip+2];
	    if(d2.user_info<PU14>().vertex_number()==-2){ // std::cout << "yay daughter 2" << std::endl;
	      d2True=true;
	    }else{
	      std::cout<<" no 2nd daughter, but enough entries "<<std::endl;
	    }
	  }else{
	    std::cout<<"no entries for 2nd daughter"<<std::endl;
	  }  
	}else{
          std::cout<<"no entries 1st daughter"<<std::endl;
	}

        if( (d1True==true) && (d2True==true)){ // std::cout << "yay both daughters" << std::endl;
          partonWd1d2pt.push_back(p.pt());
          double dr = std::sqrt(d1.squared_distance(d2));
          double kt = min(d1.pt(), d2.pt())*dr;
          drsplit.push_back(dr);
          ktsplit.push_back(kt);
          double zg = min(d1.pt(), d2.pt()) / (d1.pt() + d2.pt());
	  zgsplit.push_back(zg);
          tfesplit.push_back(1./(2.*zg*(1-zg)*p.perp()*GeVtofm*(1-cos(dr/R))));//correcting this
        }else{//to make sure that the matching still works out !
	 /* drsplit.push_back(-99);
	  ktsplit.push_back(-99);
	  zgsplit.push_back(-99);
	  tfesplit.push_back(-99);
	*/
	  std::cout<<"missing daughter(s)";
	}
        if( (d1True==true) && (d2True==false)){
          partonWd1Wod2pt.push_back(p.pt());
        }
        if( (d1True==false) && (d2True==true)){
          partonWod1Wd2pt.push_back(p.pt());
        }
        if( (d1True==false) && (d2True==false)){
          partonWod1d2pt.push_back(p.pt());
        }
      }
    }
    //std::cout<<"event "<<iev<<std::endl;

    /*
    // old version
    for(int ip = 0; ip<partons.size(); ++ip) {
      std::cout<<"N event: "<<iev<<std::endl;
      std::cout << "1st split hard parton "<<" ip: "<<ip<<" partons[ip]: "<<partons[ip]<<"\n";
      std::cout<<"partonsFirstSplit[id]: "<<partonsFirstSplit[id]<<"\n";
      std::cout<<"partonsFirstSplit[id+1]: "<<partonsFirstSplit[id+1] << std::endl;
      PseudoJet p = partons[ip];
      PseudoJet d1 = partonsFirstSplit[id];
      PseudoJet d2 = partonsFirstSplit[id+1];
      double dr = std::sqrt(d1.squared_distance(d2));
      double kt = min(d1.pt(), d2.pt())*dr; 
      drsplit.push_back(dr);
      ktsplit.push_back(kt);
      double z1 = max(d1.e(),d2.e())/p.e();
      double z2 = min(d1.e(),d2.e())/p.e();
      double zg = min(d1.pt(), d2.pt()) / (d1.pt() + d2.pt());
      tfesplit.push_back(1./(2.*zg*(1-zg)*p.perp()*GeVtofm*(1-cos(dr/R))));
      id+=2;
      
    }
    */
    
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    //std::cout << "jet clustering" << std::endl;
    // run the clustering, extract the signal jets
    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(15.))));

    // run the clustering, extract the signal charged jets
    fastjet::ClusterSequenceArea csSigCh(particlesSigCh, jet_def, area_def);
    jetCollection jetCollectionSigCh(sorted_by_pt(jet_selector(csSigCh.inclusive_jets(15.))));
    
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
    jetCollectionSig.addVector("sigJetRecur_erad",      sdcSig.getErads());
    jetCollectionSig.addVector("sigJetRecur_logdr12",   sdcSig.getLog1DRs());
    jetCollectionSig.addVector("sigJetRecur_logztheta", sdcSig.getLogzDRs());
    jetCollectionSig.addVector("sigJetRecur_nSD",       sdcSig.calculateNSD(0.0));
    jetCollectionSig.addVector("sigJetRecur_zSD",       sdcSig.calculateNSD(1.0));
    jetCollectionSig.addVector("sigJetRecur_tf",        sdcSig.getTfs());//high E & soft & collinear
    jetCollectionSig.addVector("sigJetRecur_tfe",       sdcSig.getTfes());//high E (1-cos)
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
    std::vector<int> partonmatch;
    std::vector<double> partonmatchdr;
    std::vector<fastjet::PseudoJet> sigJets =  jetCollectionSig.getJet();
    for(fastjet::PseudoJet j : sigJets) {
      int ipmin = -1;
      double drmin = 999.;
      int particlesSize = particlesMergedAll.size();
      for(int ip = 0; ip< std::min(particlesSize,6); ++ip) {
        bool d1True=false;
        bool d2True=false;
        PseudoJet d1;
        PseudoJet d2;
        PseudoJet p = particlesMergedAll[ip];
        int vtx = p.user_info<PU14>().vertex_number();
	if(vtx != -1) continue;
        if(vtx == -1) {
          if( (ip+1) < particlesMergedAll.size() ){
            d1 = particlesMergedAll[ip+1];
            if(d1.user_info<PU14>().vertex_number()==-2){ // std::cout << "yay daughter 1" << std::endl;
              d1True=true;
	    }else{
              std::cout<<" no 1st daughter for match, but enough entries "<<std::endl;
            }
            if( (ip+2) < particlesMergedAll.size() ){
              d2 = particlesMergedAll[ip+2];
              if(d2.user_info<PU14>().vertex_number()==-2){ // std::cout << "yay daughter 2" << std::endl;
                d2True=true;
              }else{
                std::cout<<" no 2nd daughter for match, but enough entries "<<std::endl;
              }
            }else{
              std::cout<<"no entries for 2nd daughter match"<<std::endl;
            }
          }else{
            std::cout<<"no entries 1st daughter match"<<std::endl;
          }
          if( (d1True==true) && (d2True==true)){
            double dr = std::sqrt(j.squared_distance(partons[ip]));// distance jet to parton 
	  //for(int ip = 0; ip<partons.size(); ++ip) {
          //double dr = p.delta_R(partons[ip]);//I could belive that smth is wrng here ..  
	    if(dr<drmin) {//though it looks all reasonable, the fact that the partons ipmin is always 0 or 1 is suspicious to me...
              drmin = dr;
              ipmin = ip;
            }
          }
        }
      }// ip loop
      //std::cout<<"drmin: "<<drmin<<"      ipmin: "<<ipmin<<std::endl;
      partonmatch.push_back(ipmin);
      partonmatchdr.push_back(drmin);
    }
    jetCollectionSig.addVector("sigJetRecur_partonMatchID", partonmatch);
    jetCollectionSig.addVector("sigJetRecur_partonMatchDr", partonmatchdr);

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
   
    //---------------------------------------------------------------------------
    //  Bkg subtraction - core dump if i dont execute this part !
    //---------------------------------------------------------------------------
   // csSubtractor csSub(R, 0., 0.1, 0.005,1);//Values from Bas - here it breaks.
   // double rJet = 0.4, double alpha = 1., double rParam = -1., double ghostArea = 0.005, double ghostRapMax = 3.0, double jetRapMax = 3.0
    csSubtractor csSub(R, 0., -1, 0.005,ghostRapMax,2.5);//Values from Marta jettoy - here it runs
    csSub.setInputParticles(particlesMerged);
    jetCollection jetCollectionCS(csSub.doSubtraction());

    //Background densities used by constituent subtraction
    
    //std::vector<double> rho;
    //std::vector<double> rhom;
    //rho.push_back(csSub.getRho());
    //rhom.push_back(csSub.getRhoM());
    
    //match CS jets to signal jets
    jetMatcher jmCS(R);
    jmCS.setBaseJets(jetCollectionCS);
    jmCS.setTagJets(jetCollectionSig);
    jmCS.matchJets();

    jmCS.reorderedToTag(jetCollectionCS);

    //---------------------------------------------------------------------------
    //   Full event constituent subtraction -- this one
    //---------------------------------------------------------------------------

    //csSubtractorFullEvent csSubFull( 0., 0.1, 0.005, 1);  // alpha, rParam, ghA, ghRapMax -> values from Bas = here it doesn't break ?!
    //csSubtractorFullEvent csSubFull( 0., 0.25, 0.005, 2.5);  // alpha, rParam, ghA, ghRapMax -> double check with bas
    // csSubFull.setRho(csSub.getRho());//  Or should this be inherited from the jet-by-jet CS above ?!
    // csSubFull.setRhom(csSub.getRhoM());
    // csSubFull.setMaxEta(1.0);
    // csSubFull.setInputParticles(particlesMerged);
    //New CS subtraction - iteratively:

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

    //find closest parton for each jet
    std::vector<int> partonmatchCSFull;
    std::vector<double> partonmatchdrCSFull;
    std::vector<fastjet::PseudoJet> csFullJets =  jetCollectionCSFull.getJet();
    for(fastjet::PseudoJet j : csFullJets) {
      int ipmin = -1;
      double drmin = 999.;
      int particlesSize = particlesMergedAll.size();
      for(int ip = 0; ip< std::min(particlesSize,6); ++ip) {
        bool d1True=false;
        bool d2True=false;
        PseudoJet d1;
        PseudoJet d2;
        PseudoJet p = particlesMergedAll[ip];
        int vtx = p.user_info<PU14>().vertex_number();
        if(vtx != -1) continue;
        if(vtx == -1) {
          if( (ip+1) < particlesMergedAll.size() ){
            d1 = particlesMergedAll[ip+1];
            if(d1.user_info<PU14>().vertex_number()==-2){
              d1True=true;
	    }else{
              std::cout<<" no 1st daughter for match csMatch, but enough entries "<<std::endl;
            }
            if( (ip+2) < particlesMergedAll.size() ){
              d2 = particlesMergedAll[ip+2];
              if(d2.user_info<PU14>().vertex_number()==-2){ 
                d2True=true;
              }else{
                std::cout<<" no 2nd daughter for match csMatch, but enough entries "<<std::endl;
              }
            }else{
              std::cout<<"no entries for 2nd daughter csMatch"<<std::endl;
            }
          }else{
            std::cout<<"no entries 1st daughter csMatch"<<std::endl;
          }
          if( (d1True==true) && (d2True==true)){
            double dr = j.delta_R(partons[ip]);
            if(dr<drmin) {
              drmin = dr;
              ipmin = ip;
            }
         }
      }
    }// ip loop
    partonmatchCSFull.push_back(ipmin);
    partonmatchdrCSFull.push_back(drmin);
    }

    jetCollectionCSFull.addVector("csFullJetRecur_partonMatchID", partonmatchCSFull);
    jetCollectionCSFull.addVector("csFullJetRecur_partonMatchDr", partonmatchdrCSFull);

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
    trwSig.addCollection("csRho",         rho);
    trwSig.addCollection("csRhom",        rhom);

    trwSig.addPartonCollection("partons",       partons);
    trwSig.addPartonCollection("partonsFirstSplit",       partonsFirstSplit);
    trwSig.addDoubleCollection("drsplit", drsplit);
    trwSig.addDoubleCollection("tfesplit", tfesplit);
    trwSig.addDoubleCollection("zgsplit", zgsplit);
    trwSig.addDoubleCollection("ktsplit", ktsplit);

    trwSig.addDoubleCollection("partonWd1Wod2pt", partonWd1Wod2pt);
    trwSig.addDoubleCollection("partonWod1Wd2pt", partonWod1Wd2pt);
    trwSig.addDoubleCollection("partonWd1d2pt", partonWd1d2pt);
    trwSig.addDoubleCollection("partonWod1d2pt", partonWod1d2pt);

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
  //std::cout << "runFromFile: " << time_in_seconds << std::endl;
}
