#define jetTreeSig_cxx
#include "jetTreeSig.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>


int returnBinNumber(double Value,int Bins,int BinLow,int BinUp){
  float binWidth = abs(BinUp-BinLow);;
  binWidth = binWidth/Bins;
 // std::cout<<"binWidth: "<<binWidth;
 // std::cout<<"BinUp: "<<BinUp<<" binLow: "<<BinLow<<std::endl;
 // std::cout<<"value: "<<Value<<" NBins: "<<Bins<<std::endl;
  int binNr = Value/binWidth;
 // std::cout<<"Bin NR: "<<binNr<<std::endl;
  return binNr;
}
void jetTreeSig::Loop(Bool_t reduced, TString generator, TString duration, TString fin, double R)
{
//   In a ROOT session, you can do:
//      root> .L jetTreeSig.C
//      root> jetTreeSig t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//    Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //Long64_t nentries = fChain->GetEntries();
   int nentries = fChain->GetEntriesFast();
   std::cout << "nentries: " << nentries << std::endl;
   TH1 *hEntries = new TH1F("hEntries","hEntries",10,0.,1.);
   
   //Jet Energy Scale 
   TH1F *h1JESCS = new TH1F("h1JESCS","hJESCS",100,0.,2.);
   TH2F *h2PtJESCS = new TH2F("h2PtJESCS","h2PtJESCS",200.,0.,400.,100,0.,2.);
   //Crosscheck - add weight ("eventWeight",   eventWeight);
   TH1F *h1JPtW = new TH1F("h1JPtW","h1JPtW",40,0.,400.);
   TH1F *h1JPt = new TH1F("h1JPt","h1JPt",40,0.,400.);
   //Crosscheck - parton pt's for 1, 2 or none daughters - partons that split/don't split
   TH1F *h1PartonNoDaughterPt = new TH1F("h1PartonNoDaughterPt","h1PartonNoDaughterPt",50,100.,400.);
   TH1F *h1Parton2DaughtersPt = new TH1F("h1Parton2DaughtersPt","h1Parton2DaughtersPt",50,100.,400.);
   TH1F *h1PartonDaughter1Pt = new TH1F("h1PartonDaughter1Pt","h1PartonDaughter1Pt",50,100.,400.);
   TH1F *h1PartonDaughter2Pt = new TH1F("h1PartonDaughter2Pt","h1PartonDaughter2Pt",50,100.,400.);

   //Differences in Pt Sig-Tau,Ca 
   TH1F *h1SigTauCs = new TH1F("h1SigTauCs","h1SigTauCs",400,-400.,400.);//sig-taucs
   TH1F *h1SigCaCs = new TH1F("h1SigCaCs","h1SigCaCs",400,-400,400.);//sig-cacs
   TH1F *h1CaTauCs = new TH1F("h1CaTauCs","h1CaTauCs",400,-400.,400.);//ca cs - tau cs
   TH1F *h1Ca = new TH1F("h1Ca","h1Ca",400,-400.,400.);//sig-cs
   TH1F *h1Tau = new TH1F("h1Tau","h1Tau",400,-400.,400.);//sig-tau

   TH1F *h1JsigTauPtW = new TH1F("h1JsigTauPtW","h1JsigTauPtW",40,0.,400.);
   TH1F *h1JsigTauPt = new TH1F("h1JsigTauPt","h1JsigTauPt",40,0.,400.);

   TH1F *h1umsigTauPt = new TH1F("h1umsigTauPt","h1umsigTauPt",40,0.,400.);//unmatched
   TH1F *h1matchedsigTauPt = new TH1F("h1matchedsigTauPt","h1matchedsigTauPt",40,0.,400.);//matched
   TH1F *h1umsigCAPt = new TH1F("h1umsigCAPt","h1umsigCAPt",40,0.,400.);
   TH1F *h1matchedsigCAPt = new TH1F("h1matchedsigCAPt","h1matchedsigCAPt",40,0.,400.);
   TH1F *h1umsigTauEta = new TH1F("h1umsigTauEta","h1umsigTauEta",40,-2,2);
   TH1F *h1matchedsigTauEta = new TH1F("h1matchedsigTauEta","h1matchedsigTauEta",40,-2,2);
   TH1F *h1umsigCAEta = new TH1F("h1umsigCAEta","h1umsigCAEta",40,-2,2);
   TH1F *h1matchedsigCAEta = new TH1F("h1matchedsigCAEta","h1matchedsigCAEta",40,-2,2);

   TH1F *h1JsigCAPtW = new TH1F("h1JsigCAPtW","h1JsigCAPtW",40,0.,400.);
   TH1F *h1JsigCAPt = new TH1F("h1JsigCAPt","h1JsigCAPt",40,0.,400.);

   //N formation times
   TH1F *h1Ntf2ZCA = new TH1F("h1Ntf2ZCA","h1Ntf2ZCA",100,0.,10);//1->3
   TH1F *h1Ntf2ZTau = new TH1F("h1Ntf2ZTau","h1Ntf2ZTau",100,0.,10);
   TH1F *h1Ntf2CA = new TH1F("h1Ntf2CA","h1Ntf2CA",100,0.,10);
   TH1F *h1Ntf2Tau = new TH1F("h1Ntf2Tau","h1Ntf2Tau",100,0.,10);
   TH1F *h1NtfZCA = new TH1F("h1NtfZCA","h1NtfZCA",100,0.,10);//1->2
   TH1F *h1NtfZTau = new TH1F("h1NtfZTau","h1NtfZTau",100,0.,10);
   TH1F *h1NtfCA = new TH1F("h1NtfCA","h1NtfCA",100,0.,10);
   TH1F *h1NtfTau = new TH1F("h1NtfTau","h1NtfTau",100,0.,10);
   //cs
   TH1F *h1Ncstf2ZCA = new TH1F("h1Ncstf2ZCA","h1Ncstf2ZCA",100,0.,10);//1->3
   TH1F *h1Ncstf2ZTau = new TH1F("h1Ncstf2ZTau","h1Ncstf2ZTau",100,0.,10);
   TH1F *h1Ncstf2CA = new TH1F("h1Ncstf2CA","h1Ncstf2CA",100,0.,10);
   TH1F *h1Ncstf2Tau = new TH1F("h1Ncstf2Tau","h1Ncstf2Tau",100,0.,10);
   TH1F *h1NcstfZCA = new TH1F("h1NcstfZCA","h1NcstfZCA",100,0.,10);//1->2
   TH1F *h1NcstfZTau = new TH1F("h1NcstfZTau","h1NcstfZTau",100,0.,10);
   TH1F *h1NcstfCA = new TH1F("h1NcstfCA","h1NcstfCA",100,0.,10);
   TH1F *h1NcstfTau = new TH1F("h1NcstfTau","h1NcstfTau",100,0.,10);
   //Formation time
   TH2F *h2TfeRespCA = new TH2F("h2TfeRespCA","h2TfeRespCA",40,0.,9., 40,0.,5.);//1->2 tfe cs/sig
   TH2F *h2TfeRespTau = new TH2F("h2TfeRespTau","h2TfeRespTau",40,0.,9.,40,0.,5.);
   TH2F *h2Tfe2RespCA = new TH2F("h2Tfe2RespCA","h2Tfe2RespCA",40,0.,9., 40,0.,5.);//1->3 tfe
   TH2F *h2Tfe2RespTau = new TH2F("h2Tfe2RespTau","h2Tfe2RespTau",40,0.,9., 40,0.,5.);
   //with zcut
   TH2F *h2TfeRespCAZcut = new TH2F("h2TfeRespCAZcut","h2TfeRespCAZcut",40,0.,9., 40,0.,5.);//1->2 tfe cs/sig
   TH2F *h2TfeRespTauZcut = new TH2F("h2TfeRespTauZcut","h2TfeRespTauZcut",40,0.,9., 40,0.,5.);
   TH2F *h2Tfe2RespCAZcut = new TH2F("h2Tfe2RespCAZcut","h2Tfe2RespCAZcut",40,0.,9., 40,0.,5.);//1->3 tfe
   TH2F *h2Tfe2RespTauZcut = new TH2F("h2Tfe2RespTauZcut","h2Tfe2RespTauZcut",40,0.,9., 40,0.,5.);

   //theta vs zg
   TH2F *h2ThetaZgCAZcut = new TH2F("h2ThetaZgCAZcut","h2ThetaZgCAZcut",50,0.,1, 50,0.,1.);//
   TH2F *h2ThetaZgTauZcut = new TH2F("h2ThetaZgTauZcut","h2ThetaZgTauZcut",50,0.,1, 50,0.,1.);//
   TH2F *h2ThetaZgCA = new TH2F("h2ThetaZgCA","h2ThetaZgCA",50,0.,1, 50,0.,1.);//
   TH2F *h2ThetaZgTau = new TH2F("h2ThetaZgTau","h2ThetaZgTau",50,0.,1, 50,0.,1.);//

   TH1F *h1LogTfSig = new TH1F("h1LogTfSig","h1LogTfSig",50,-10.,10.);
   TH1F *h1LogTfSigTau = new TH1F("h1LogTfSigTau","h1LogTfSigTau",50,-10.,10.);
   TH1F *h1LogTfSigZ = new TH1F("h1LogTfSigZ","h1LogTfSigZ",50,-10.,10.);
   TH1F *h1LogTfSigZTau = new TH1F("h1LogTfSigZTau","h1LogTfSigZTau",50,-10.,10.);
   TH1F *h1LogTfSigToPartonZ = new TH1F("h1LogTfSigToPartonZ","h1LogTfSigToPartonZ",50,-10.,10.);
   TH1F *h1LogTfSigToPartonZTau = new TH1F("h1LogTfSigToPartonZTau","h1LogTfSigToPartonZTau",50,-10.,10.);

   TH2F *h2SigTfeTauTfe = new TH2F("h2SigTfeTauTfe","h2SigTfeTauTfe",100,0.,10.,100,0.,10.);
   TH2F *h2SigTfeTauTfeZcut = new TH2F("h2SigTfeTauTfeZcut","h2SigTfeTauTfeZcut",100,0.,10.,100,0.,10.);
   TH2F *h2SigTfeTauTfeCs = new TH2F("h2SigTfeTauTfeCS","h2SigTfeTauTfeCs",100,0.,10.,100,0.,10.);
   TH2F *h2SigTfeTauTfeZcutCs = new TH2F("h2SigTfeTauTfeZcutCs","h2SigTfeTauTfeZcutCs",100,0.,10.,100,0.,10.);
   
   
   TH2F *h2SigRgTauRgZcut = new TH2F("h2SigRgTauRgZcut","h2SigRgTauRgZcut",100,0.,1.,100,0.,1.);// Rg corr
   TH2F *h2SigZgTauZgZcut = new TH2F("h2SigZgTauZgZcut","h2SigZgTauZgZcut",100,0.,0.6,100,0.,0.6);// Zg corr
   //Pythia information
   TH3F *h3tfs = new TH3F("h3tfs","h3tfs",11,0.,10.,100,0.,20.,100,0,20);// len, max, sum
   TH3F *h3tfes = new TH3F("h3tfes","h3tfes",11,0.,10.,100,0.,20.,100,0,20);// len, max, sum
   TH3F *h3drs = new TH3F("h3drs","h3drs",11,0.,10.,100,0.,20.,100,0,20);// len, max, sum
   TH3F *h3tftsfesdrs = new TH3F("h3tftsfesdrs","h3tftsfesdrs",100,0.,20.,100,0.,20.,100,0,10);// tf, tfe, dr
   TH3F *h3jetsplittfs = new TH3F("h3jetsplittfs","h3jetsplittfs",100, 0, 99, 10,0.,9.,100,0.,20.);// jet_i, split_j, tf
   TH3F *h3jetsplittfes = new TH3F("h3jetsplittfes","h3jetsplittfes",100, 0, 99, 10,0.,9.,100,0.,20.);// jet_i, split_j, tfe
   TH3F *h3jetsplitdrs = new TH3F("h3jetsplitdrs","h3jetsplitdrs",100, 0, 99, 10,0.,9.,100,0.,20.);// jet_i, split_j, dr

   //Number of dropped branches before 1st split yields the dr of the dropped emissions)
   TH2F *h2DroppedDrTauZcut = new TH2F("h2DroppedDrTauZcut","h2DroppedDrTauZcut",21,0.,20.,100,0.,1.);
   TH2F *h2DroppedPtTauZcut = new TH2F("h2DroppedPtTauZcut","h2DroppedPtTauZcut",21,0.,20.,400,0.,400.);
   TH2F *h2DroppedKtTauZcut = new TH2F("h2DroppedKtTauZcut","h2DroppedKtTauZcut",21,0.,20.,400,0.,20.);
   TH2F *h2DroppedLPTauZcut = new TH2F("h2DroppedLPTauZcut","h2DroppedLPTauZcut",100,0.,5.,100,-10.,5.);// log dr, log kt  100, 0,5, 100, -10, 2)
   TH2F *h2DroppedDrZcut = new TH2F("h2DroppedDrZcut","h2DroppedDrZcut",21,0.,20.,100,0.,1.);
   TH2F *h2DroppedPtZcut = new TH2F("h2DroppedPtZcut","h2DroppedPtZcut",21,0.,20.,400,0.,400.);
   TH2F *h2DroppedKtZcut = new TH2F("h2DroppedKtZcut","h2DroppedKtZcut",21,0.,20.,400,0.,20.);
   TH2F *h2DroppedLPZcut = new TH2F("h2DroppedLPZcut","h2DroppedLPZcut",100,0.,5.,100,-10.,5.);

   //Compare 1->2 with 1->3 body formation time
   TH2F *h2CATfe2TfeSig = new TH2F("h2CATfe2TfeSig", "h2CATfe2TfeSig", 100, 0,20, 100, 0,20);
   TH2F *h2TauTfe2TfeSig = new TH2F("h2TauTfe2TfeSig", "h2TauTfe2TfeSig", 100, 0,20, 100, 0,20);
   TH2F *h2CATfe2Tfe = new TH2F("h2CATfe2Tfe", "h2CATfe2Tfe", 100, 0,20, 100, 0,20);
   TH2F *h2TauTfe2Tfe = new TH2F("h2TauTfe2Tfe", "h2TauTfe2Tfe", 100, 0,20, 100, 0,20);
   TH2F *h2CATfe2TfeCs = new TH2F("h2CATfe2TfeCs", "h2CATfe2TfeCs", 100, 0,20, 100, 0,20);
   TH2F *h2TauTfe2TfeCs = new TH2F("h2TauTfe2TfeCs", "h2TauTfe2TfeCs", 100, 0,20, 100, 0,20);
   TH2F *h2CATfe2TfeCsUg = new TH2F("h2CATfe2TfeCsUg", "h2CATfe2TfeCsUg", 100, 0,20, 100, 0,20);
   TH2F *h2TauTfe2TfeCsUg = new TH2F("h2TauTfe2TfeCsUg", "h2TauTfe2TfeCsUg", 100, 0,20, 100, 0,20);

   TH2F *h2SumCATfe2Tfe = new TH2F("h2SumCATfe2Tfe", "h2SumCATfe2Tfe", 100, 0,20, 100, 0,20);
   TH2F *h2SumTauTfe2Tfe = new TH2F("h2SumTauTfe2Tfe", "h2SumTauTfe2Tfe", 100, 0,20, 100, 0,20);
   TH2F *h2SumCATfe2TfeCs = new TH2F("h2SumCATfe2TfeCs", "h2SumCATfe2TfeCs", 100, 0,20, 100, 0,20);
   TH2F *h2SumTauTfe2TfeCs = new TH2F("h2SumTauTfe2TfeCs", "h2SumTauTfe2TfeCs", 100, 0,20, 100, 0,20);
   //Pt effects in 2D
   TH2F *h2SigTauCs = new TH2F("h2SigTauCs", "h2SigTauCs", 100, 0,200, 100, 0,200);
   TH2F *h2SigCaCs = new TH2F("h2SigCaCs", "h2SigCaCs", 100, 0,200, 100, 0,200);
   TH2F *h2CaTauCs = new TH2F("h2CaTauCs", "h2CaTauCs", 100, 0,200, 100, 0,200);
   TH2F *h2Ca = new TH2F("h2Ca", "h2Ca", 100, 0,200, 100, 0,200);
   TH2F *h2Tau = new TH2F("h2Tau", "h2Tau", 100, 0,200, 100, 0,200);

   //Lund Planes Novel Tools Fig. 3
   TH2F *h2CALundPlaneNT = new TH2F("h2CALundPlaneNT", "h2CALundPlaneNT", 100, 0, 5, 100, -10,2);
   TH2F *h2TauLundPlaneNT = new TH2F("h2TauLundPlaneNT", "h2TauLundPlaneNT", 100, 0,5, 100, -10, 2);
   //Lund Planes single Jet 
   //log (kt) = log(pt*dr/R) vs log(tfe)   
   TH3F *h3CALundPlaneSigTfeSingleIncl = new TH3F("h3CALundPlaneSigTfeSingleIncl", "h3CALundPlaneSigTfeSingleIncl", 100, -6,12, 100, -10, 5, 21, 0, 20);
   TH3F *h3TauLundPlaneSigTfeSingleIncl = new TH3F("h3TauLundPlaneSigTfeSingleIncl", "h3TauLundPlaneSigTfeSingleIncl", 100, -6,12, 100, -10, 5, 21, 0, 20);
   TH3F *h3CALundPlaneSigTfeSingleInclZ = new TH3F("h3CALundPlaneSigTfeSingleInclZ", "h3CALundPlaneSigTfeSingleInclZ", 100, -6,12, 100, -10, 5, 21, 0, 20);
   TH3F *h3TauLundPlaneSigTfeSingleInclZ = new TH3F("h3TauLundPlaneSigTfeSingleInclZ", "h3TauLundPlaneSigTfeSingleInclZ", 100, -6,12, 100, -10, 5, 21, 0, 20);
   TH2F *h2CALundPlaneSigTfeSingle = new TH2F("h2CALundPlaneSigTfeSingle", "h2CALundPlaneSigTfeSingle", 100, -6,12, 100, -10, 5);
   TH2F *h2TauLundPlaneSigTfeSingle = new TH2F("h2TauLundPlaneSigTfeSingle", "h2TauLundPlaneSigTfeSingle", 100, -6,12, 100, -10, 5);
   //log (kt) = log(pt*dr/R) vs log(R/dr)
   TH3F *h3CALundPlaneSigSingleIncl = new TH3F("h3CALundPlaneSigSingleIncl", "h3CALundPlaneSigSingleIncl", 100, -6,12, 100, -10, 5, 21, 0, 20);
   TH3F *h3TauLundPlaneSigSingleIncl = new TH3F("h3TauLundPlaneSigSingleIncl", "h3TauLundPlaneSigSingleIncl", 100, -6,12, 100, -10, 5, 21, 0, 20);
   TH3F *h3CALundPlaneSigSingleInclZ = new TH3F("h3CALundPlaneSigSingleInclZ", "h3CALundPlaneSigSingleInclZ", 100, -6,12, 100, -10, 5, 21, 0, 20);
   TH3F *h3TauLundPlaneSigSingleInclZ = new TH3F("h3TauLundPlaneSigSingleInclZ", "h3TauLundPlaneSigSingleInclZ", 100, -6,12, 100, -10, 5, 21, 0, 20);
   TH2F *h2CALundPlaneSigSingle = new TH2F("h2CALundPlaneSigSingle", "h2CALundPlaneSigSingle", 100, -6,12, 100, -10, 5);
   TH2F *h2TauLundPlaneSigSingle = new TH2F("h2TauLundPlaneSigSingle", "h2TauLundPlaneSigSingle", 100, -6,12, 100, -10, 5);
   // Ungroomed Lund Planes
   // true parton info log (kt) = log(pt*dr) vs log(R/dr)
   TH2F *h2PartonLogTfKt = new TH2F("h2PartonLogTfKt", "h2PartonLogTfKt", 100, -6,12, 100, -10, 5);
   TH2F *h2PartonLogDrKt = new TH2F("h2PartonLogDrKt", "h2PartonLogDrKt", 100, -6,10, 100, -10, 5);
   TH2F *h2PartonLundPlaneNT = new TH2F("h2PartonLundPlaneNT", "h2PartonLundPlaneNT", 100, -6,12, 100, -10, 5);
   
   //with zcut
   TH2F *h2PartonLogTfKtZcut = new TH2F("h2PartonLogTfKtZcut", "h2PartonLogTfKtZcut", 100, -6,12, 100, -10, 5);
   TH2F *h2PartonLogDrKtZcut = new TH2F("h2PartonLogDrKtZcut", "h2PartonLogDrKtZcut", 100, -6,10, 100, -10, 5);
   //log (kt) = log(pt*dr) vs log(R/dr)
   TH2F *h2CALundPlaneSig = new TH2F("h2CALundPlaneSig", "h2CALundPlaneSig", 100, -6,10, 100, -10, 5);
   TH2F *h2TauLundPlaneSig = new TH2F("h2TauLundPlaneSig", "h2TauLundPlaneSig", 100, -6,10, 100, -10, 5);
   TH2F *h2CALundPlaneCs = new TH2F("h2CALundPlaneCs", "h2CALundPlaneCs", 100, -6, 10, 100, -10, 5);
   TH2F *h2TauLundPlaneCs = new TH2F("h2TauLundPlaneCs", "h2TauLundPlaneCs", 100, -6,10, 100, -10, 5);
   //log (kt) = log(pt*dr) vs log(tfe)   
   TH2F *h2CALundPlaneSigTfe = new TH2F("h2CALundPlaneSigTfe", "h2CALundPlaneSigTfe", 100, -6,12, 100, -10, 5);
   TH2F *h2TauLundPlaneSigTfe = new TH2F("h2TauLundPlaneSigTfe", "h2TauLundPlaneSigTfe", 100, -6,12, 100, -10, 5);
   TH2F *h2CALundPlaneCsTfe = new TH2F("h2CALundPlaneCsTfe", "h2CALundPlaneCsTfe", 100, -6,12, 100, -10, 5);
   TH2F *h2TauLundPlaneCsTfe = new TH2F("h2TauLundPlaneCsTfe", "h2TauLundPlaneCsTfe", 100, -6,12, 100, -10, 5);
   //Lund Planes log (kt) = log(pt*dr) vs log(R/dr)
   TH2F *h2LogLundPlane = new TH2F("h2LogLundPlane", "h2LogLundPlane", 100, 0,10, 100, -10, 5);
   TH2F *h2LogLundPlaneTf = new TH2F("h2LogLundPlaneTf", "h2LogLundPlaneTf", 100, -6,12, 100, -10, 5);
   TH2F *h2LogLundPlaneTfe = new TH2F("h2LogLundPlaneTfe", "h2LogLundPlaneTfe", 100, -6,12, 100, -10, 5);
   TH2F *h2LogLundPlaneTau = new TH2F("h2LogLundPlaneTau", "h2LogLundPlaneTau", 100, 0,10, 100, -10, 5);
   TH2F *h2LogLundPlaneTauTf = new TH2F("h2LogLundPlaneTauTf", "h2LogLundPlaneTauTf", 100, -6,12, 100, -10, 5);
   TH2F *h2LogLundPlaneTauTfe = new TH2F("h2LogLundPlaneTauTfe", "h2LogLundPlaneTauTfe", 100, -6,12, 100, -10, 5);
   TH2F *h2LogLundPlaneCs = new TH2F("h2LogLundPlaneCs", "h2LogLundPlaneCs", 100, 0,10, 100, -10, 5);
   TH2F *h2LogLundPlaneTfCs = new TH2F("h2LogLundPlaneTfCs", "h2LogLundPlaneTfCs", 100, -6,12, 100, -10, 5);
   TH2F *h2LogLundPlaneTfeCs = new TH2F("h2LogLundPlaneTfeCs", "h2LogLundPlaneTfeCs", 100, -6,12, 100, -10, 5);
   TH2F *h2LogLundPlaneTauCs = new TH2F("h2LogLundPlaneTauCs", "h2LogLundPlaneTauCs", 100, 0,10, 100, -10, 5);
   TH2F *h2LogLundPlaneTauTfCs = new TH2F("h2LogLundPlaneTauTfCs", "h2LogLundPlaneTauTfCs", 100, -6,12, 100, -10, 5);
   TH2F *h2LogLundPlaneTauTfeCs = new TH2F("h2LogLundPlaneTauTfeCs", "h2LogLundPlaneTauTfeCs", 100, -6,12, 100, -10, 5);
   TH2F *h2LogLundPlaneBFS = new TH2F("h2LogLundPlaneBFS", "h2LogLundPlaneBFS", 100, 0,10, 100, -10, 5);
   TH2F *h2LogLundPlaneTauBFS = new TH2F("h2LogLundPlaneTauBFS", "h2LogLundPlaneTauBFS", 100, 0,10, 100, -10, 5);

   TH2F *h2Log1stLundPlane = new TH2F("h2Log1stLundPlane", "h2Log1stLundPlane", 100, 0,10, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTf = new TH2F("h2Log1stLundPlaneTf", "h2Log1stLundPlaneTf", 100, -6,12, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTfe = new TH2F("h2Log1stLundPlaneTfe", "h2Log1stLundPlaneTfe", 100, -6,12, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTau = new TH2F("h2Log1stLundPlaneTau", "h2Log1stLundPlaneTau", 100, 0,10, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTauTf = new TH2F("h2Log1stLundPlaneTauTf", "h2Log1stLundPlaneTauTf", 100, -6,12, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTauTfe = new TH2F("h2Log1stLundPlaneTauTfe", "h2Log1stLundPlaneTauTfe", 100, -6,12, 100, -10, 5);
   TH2F *h2Log1stLundPlaneCs = new TH2F("h2Log1stLundPlaneCs", "h2Log1stLundPlaneCs", 100, 0,10, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTfCs = new TH2F("h2Log1stLundPlaneTfCs", "h2Log1stLundPlaneTfCs", 100, -6,12, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTfeCs = new TH2F("h2Log1stLundPlaneTfeCs", "h2Log1stLundPlaneTfeCs", 100, -6,12, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTauCs = new TH2F("h2Log1stLundPlaneTauCs", "h2Log1stLundPlaneTauCs", 100, 0,10, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTauTfCs = new TH2F("h2Log1stLundPlaneTauTfCs", "h2Log1stLundPlaneTauTfCs", 100, -6,12, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTauTfeCs = new TH2F("h2Log1stLundPlaneTauTfeCs", "h2Log1stLundPlaneTauTfeCs", 100, -6,12, 100, -10, 5);
   TH2F *h2Log1stLundPlaneBFS = new TH2F("h2Log1stLundPlaneBFS", "h2Log1stLundPlaneBFS", 100, 0,10, 100, -10, 5);
   TH2F *h2Log1stLundPlaneTauBFS = new TH2F("h2Log1stLundPlaneTauBFS", "h2Log1stLundPlaneTauBFS", 100, 0,10, 100, -10, 5);
   //Number of calculated formation times in jet
   TH1F *h1SigNCaTf = new TH1F("h1SigNCaTf","h1SigNCaTf",20,0.,20.);
   TH1F *h1SigNCaTfZcut = new TH1F("h1SigNCaTfZcut","h1SigNCaTfZcut",20,0.,20.);
   TH1F *h1SigNTauTf = new TH1F("h1SigNTauTf","h1SigNTauTf",20,0.,20.);
   TH1F *h1SigNTauTfZcut = new TH1F("h1SigNTauTfZcut","h1SigNTauTfZcut",20,0.,20.);
   //Dropped branches
   TH1F *h1DroppedBrTau = new TH1F("h1DroppedBrTau","h1DroppedBrTau",20,0.,20);
   TH1F *h1DroppedBr = new TH1F("h1DroppedBr","h1DroppedBr",20,0.,20);
   //Nsubjettiness (dy)
   TH2F *h2Nsubj21 = new TH2F("h2Nsubj21","h2Nsubj21",200,0.,1.,200,0.,1.);
   TH2F *h2Nsubj32 = new TH2F("h2Nsubj32","h2Nsubj32",200,0.,1.,200,0.,1.);
   TH2F *h2TauNsubj21 = new TH2F("h2TauNsubj21","h2TauNsubj21",200,0.,10.,200,0.,1.);
   //Nsubjettiness (ca/tau)
   TH2F *h2CaTauNsubj21 = new TH2F("h2CaTauNsubj21","h2CaTauNsubj21",20,0.,1.,20,0.,1.);
   TH2F *h2CaTauNsubj32 = new TH2F("h2CaTauNsubj32","h2CaTauNsubj32",20,0.,1.,20,0.,1.);
   //Formation time scatter plots
   TH2F *h2PartonSigTf = new TH2F("h2PartonSigTf","h2PartonSigTf",100,0.,5.,100,0.,5.);
   TH2F *h2PartonSigTfTau = new TH2F("h2PartonSigTfTau","h2PartonSigTfTau",100,0.,5.,100,0.,5.);
   TH2F *h2PartonSigTfZcut = new TH2F("h2PartonSigTfZcut","h2PartonSigTfZcut",100,0.,5.,100,0.,5.);
   TH2F *h2PartonSigTfZcutTau = new TH2F("h2PartonSigTfZcutTau","h2PartonSigTfZcutTau",100,0.,5.,100,0.,5.);

   TH2F *h2PartonSigTfeZcut = new TH2F("h2PartonSigTfeZcut","h2PartonSigTfeZcut",100,0.,5.,100,0.,5.);
   TH2F *h2PartonSigTfeZcutTau = new TH2F("h2PartonSigTfeZcutTau","h2PartonSigTfeZcutTau",100,0.,5.,100,0.,5.);

   TH2F *h2SigTfTfe = new TH2F("h2SigTfTfe","h2SigTfTfe",100,0.,5.,100,0.,5.);//tf vs tfe
   TH2F *h2SigTfTfeTau = new TH2F("h2SigTfTfeTau","h2SigTfTfeTau",100,0.,5.,100,0.,5.);//tf vs tfe
   TH2F *h2SigTfTfeZcut = new TH2F("h2SigTfTfeZcut","h2SigTfTfeZcut",100,0.,5.,100,0.,5.);//tf vs tfe
   TH2F *h2SigTfTfeZcutTau = new TH2F("h2SigTfTfeZcutTau","h2SigTfTfeZcutTau",100,0.,5.,100,0.,5.);//tf vs tfe

   TH2F *h2CsPartonTfe = new TH2F("h2CsPartonTfe","h2CsPartonTfe",100,0.,5.,100,0.,5.);// parton vs cs ca- zcut
   TH2F *h2CsPartonTfeTau = new TH2F("h2CsPartonTfeTau","h2CsPartonTfeTau",100,0.,5.,100,0.,5.);//parton vs cs tau -zcut 
   TH2F *h2CsSigTfeZcut = new TH2F("h2CsSigTfeZcut","h2CsSigTfeZcut",100,0.,5.,100,0.,5.);// cs vs sig ca
   TH2F *h2CsSigTfeZcutTau = new TH2F("h2CsSigTfeZcutTau","h2CsSigTfeZcutTau",100,0.,5.,100,0.,5.);// cs vs sig tau
   

   TList *fSplitList = new TList();
   fSplitList->SetName("fSplitList");
   fSplitList->SetOwner(kTRUE);
   const int Nsplit = 2;

   TH1F *h1TfSigTau[Nsplit]; 
   TH1F *h1TfSigCA[Nsplit];  
   TH1F *h1ZgSigCA[Nsplit]; 
   TH1F *h1ZgSigTau[Nsplit];
   TH1F *h1RgCASig[Nsplit];
   TH1F *h1RgTauSig[Nsplit];
   TH1F *h1EradSigCA[Nsplit];
   TH1F *h1EradSigTau[Nsplit];
   TH1F *h1LogDr12SigCA[Nsplit];
   TH1F *h1LogDr12SigTau[Nsplit];
   TH1F *h1LogZthetaSigCA[Nsplit];
   TH1F *h1LogZthetaSigTau[Nsplit];

   for (int nsplit = 0; nsplit < 2; nsplit++){
    h1TfSigTau[nsplit] = NULL; 
    h1TfSigCA[nsplit] = NULL; 
    h1ZgSigCA[nsplit] = NULL;
    h1ZgSigTau[nsplit] = NULL;
    h1RgCASig[nsplit] = NULL;
    h1RgTauSig[nsplit] = NULL;
    h1EradSigCA[nsplit] = NULL;
    h1EradSigTau[nsplit] = NULL;
    h1LogDr12SigCA[nsplit] = NULL;
    h1LogDr12SigTau[nsplit] = NULL;
    h1LogZthetaSigCA[nsplit] = NULL;
    h1LogZthetaSigTau[nsplit] = NULL;
   }


   TList *fBins = new TList();
   fBins->SetName("fBins");
   fBins->SetOwner(kTRUE);
   const int NBins = 10;

   //tau1->3/tau1->2 in bins of tau1->2
   TH1F *h1CATfe2TfeSig[NBins];
   TH1F *h1TauTfe2TfeSig[NBins];
   TH1F *h1CATfe2Tfe[NBins];
   TH1F *h1TauTfe2Tfe[NBins];
   TH1F *h1CATfe2TfeCs[NBins];
   TH1F *h1TauTfe2TfeCs[NBins];
   TH1F *h1CATfe2TfeCsUg[NBins];
   TH1F *h1TauTfe2TfeCsUg[NBins];
   //tau1->2/tau1->3 in bins of tau1->3
   TH1F *h1CATfeTfe2Sig[NBins];
   TH1F *h1TauTfeTfe2Sig[NBins];
   TH1F *h1CATfeTfe2[NBins];
   TH1F *h1TauTfeTfe2[NBins];
   TH1F *h1CATfeTfe2Cs[NBins];
   TH1F *h1TauTfeTfe2Cs[NBins];
   TH1F *h1CATfeTfe2CsUg[NBins];
   TH1F *h1TauTfeTfe2CsUg[NBins];
   //(tau1->3-tau1->2)/tau1->2 in bins of tau1->2
   TH1F *hdiff1CATfe2TfeSig[NBins];
   TH1F *hdiff1TauTfe2TfeSig[NBins];
   TH1F *hdiff1CATfe2Tfe[NBins];
   TH1F *hdiff1TauTfe2Tfe[NBins];
   TH1F *hdiff1CATfe2TfeCs[NBins];
   TH1F *hdiff1TauTfe2TfeCs[NBins];
   TH1F *hdiff1CATfe2TfeCsUg[NBins];
   TH1F *hdiff1TauTfe2TfeCsUg[NBins];
   //(tau1->2-tau1->3)/tau1->3 in bins of tau1->3
   TH1F *hdiff1CATfeTfe2Sig[NBins];
   TH1F *hdiff1TauTfeTfe2Sig[NBins];
   TH1F *hdiff1CATfeTfe2[NBins];
   TH1F *hdiff1TauTfeTfe2[NBins];
   TH1F *hdiff1CATfeTfe2Cs[NBins];
   TH1F *hdiff1TauTfeTfe2Cs[NBins];
   TH1F *hdiff1CATfeTfe2CsUg[NBins];
   TH1F *hdiff1TauTfeTfe2CsUg[NBins];

   for (int NBin = 0; NBin < NBins; NBin++){
    h1CATfe2TfeSig[NBin] = NULL;
    h1TauTfe2TfeSig[NBin] = NULL;
    h1CATfe2Tfe[NBin] = NULL;
    h1TauTfe2Tfe[NBin] = NULL;
    h1CATfe2TfeCs[NBin] = NULL;
    h1TauTfe2TfeCs[NBin] = NULL;
    h1CATfe2TfeCsUg[NBin] = NULL;
    h1TauTfe2TfeCsUg[NBin] = NULL;
    h1CATfeTfe2Sig[NBin] = NULL;
    h1TauTfeTfe2Sig[NBin] = NULL;
    h1CATfeTfe2[NBin] = NULL;
    h1TauTfeTfe2[NBin] = NULL;
    h1CATfeTfe2Cs[NBin] = NULL;
    h1TauTfeTfe2Cs[NBin] = NULL;
    h1CATfeTfe2CsUg[NBin] = NULL;
    h1TauTfeTfe2CsUg[NBin] = NULL;

    hdiff1CATfe2TfeSig[NBin] = NULL;
    hdiff1TauTfe2TfeSig[NBin] = NULL;
    hdiff1CATfe2Tfe[NBin] = NULL;
    hdiff1TauTfe2Tfe[NBin] = NULL;
    hdiff1CATfe2TfeCs[NBin] = NULL;
    hdiff1TauTfe2TfeCs[NBin] = NULL;
    hdiff1CATfe2TfeCsUg[NBin] = NULL;
    hdiff1TauTfe2TfeCsUg[NBin] = NULL;
    hdiff1CATfeTfe2Sig[NBin] = NULL;
    hdiff1TauTfeTfe2Sig[NBin] = NULL;
    hdiff1CATfeTfe2[NBin] = NULL;
    hdiff1TauTfeTfe2[NBin] = NULL;
    hdiff1CATfeTfe2Cs[NBin] = NULL;
    hdiff1TauTfeTfe2Cs[NBin] = NULL;
    hdiff1CATfeTfe2CsUg[NBin] = NULL;
    hdiff1TauTfeTfe2CsUg[NBin] = NULL;
   }

   
   TList *fEmmission = new TList();
   fEmmission->SetName("fEmmission");
   fEmmission->SetOwner(kTRUE);
   const int Nemission = 9;

   TH3F *h3CALundPlaneSigTfeSingle[Nemission]; 
   TH3F *h3TauLundPlaneSigTfeSingle[Nemission]; 
   TH3F *h3CALundPlaneSigSingle[Nemission]; 
   TH3F *h3TauLundPlaneSigSingle[Nemission]; 

   TH3F *h3CALundPlaneSigTfeSingleZ[Nemission]; 
   TH3F *h3TauLundPlaneSigTfeSingleZ[Nemission]; 
   TH3F *h3CALundPlaneSigSingleZ[Nemission]; 
   TH3F *h3TauLundPlaneSigSingleZ[Nemission]; 

   for (int emission = 0; emission < Nemission; emission++){
    h3CALundPlaneSigTfeSingle[emission] = NULL; 
    h3TauLundPlaneSigTfeSingle[emission] = NULL; 
    h3CALundPlaneSigSingle[emission] = NULL; 
    h3TauLundPlaneSigSingle[emission] = NULL; 

    h3CALundPlaneSigTfeSingleZ[emission] = NULL; 
    h3TauLundPlaneSigTfeSingleZ[emission] = NULL; 
    h3CALundPlaneSigSingleZ[emission] = NULL; 
    h3TauLundPlaneSigSingleZ[emission] = NULL; 
   }

   for (int Ne = 0; Ne < Nemission; Ne++){
    h3CALundPlaneSigTfeSingle[Ne] = new TH3F(Form("h3CALundPlaneSigTfeSingle[%d]",Ne), Form("h3CALundPlaneSigTfeSingle[%d]",Ne), 100, -6,12, 100, -10, 5, 21, 0, 20);
    fEmmission->Add(h3CALundPlaneSigTfeSingle[Ne]);
    h3TauLundPlaneSigTfeSingle[Ne] = new TH3F(Form("h3TauLundPlaneSigTfeSingle[%d]",Ne), Form("h3TauLundPlaneSigTfeSingle[%d]",Ne), 100, -6,12, 100, -10, 5, 21, 0, 20);
    fEmmission->Add(h3TauLundPlaneSigTfeSingle[Ne]);
    h3CALundPlaneSigSingle[Ne] = new TH3F(Form("h3CALundPlaneSigSingle[%d]",Ne), Form("h3CALundPlaneSigSingle[%d]",Ne), 100, -6,12, 100, -10, 5, 21, 0, 20);
    fEmmission->Add(h3CALundPlaneSigSingle[Ne]);
    h3TauLundPlaneSigSingle[Ne] = new TH3F(Form("h3TauLundPlaneSigSingle[%d]",Ne), Form("h3TauLundPlaneSigSingle[%d]",Ne), 100, -6,12, 100, -10, 5, 21, 0, 20);
    fEmmission->Add(h3TauLundPlaneSigSingle[Ne]);

    h3CALundPlaneSigTfeSingleZ[Ne] = new TH3F(Form("h3CALundPlaneSigTfeSingleZ[%d]",Ne), Form("h3CALundPlaneSigTfeSingleZ[%d]",Ne), 100, -6,12, 100, -10, 5, 21, 0, 20);
    fEmmission->Add(h3CALundPlaneSigTfeSingleZ[Ne]);
    h3TauLundPlaneSigTfeSingleZ[Ne] = new TH3F(Form("h3TauLundPlaneSigTfeSingleZ[%d]",Ne), Form("h3TauLundPlaneSigTfeSingleZ[%d]",Ne), 100, -6,12, 100, -10, 5, 21, 0, 20);
    fEmmission->Add(h3TauLundPlaneSigTfeSingleZ[Ne]);
    h3CALundPlaneSigSingleZ[Ne] = new TH3F(Form("h3CALundPlaneSigSingleZ[%d]",Ne), Form("h3CALundPlaneSigSingleZ[%d]",Ne), 100, -6,12, 100, -10, 5, 21, 0, 20);
    fEmmission->Add(h3CALundPlaneSigSingleZ[Ne]);
    h3TauLundPlaneSigSingleZ[Ne] = new TH3F(Form("h3TauLundPlaneSigSingleZ[%d]",Ne), Form("h3TauLundPlaneSigSingleZ[%d]",Ne), 100, -6,12, 100, -10, 5, 21, 0, 20);
    fEmmission->Add(h3TauLundPlaneSigSingleZ[Ne]);
   }


   for (int NBin = 0; NBin < NBins; NBin++){// from 0 to 99
    //tau1->3/tau1->2 in bins of tau1->2
    h1CATfe2TfeSig[NBin] = new TH1F(Form("h1CATfe2TfeSig[%d]",NBin), Form("h1CATfe2TfeSig[%d]",NBin), 50, 0,4);//0,1 bin width
    fBins->Add(h1CATfe2TfeSig[NBin]);
    h1TauTfe2TfeSig[NBin] = new TH1F(Form("h1TauTfe2TfeSig[%d]",NBin), Form("h1TauTfe2TfeSig[%d]",NBin), 50, 0,4);
    fBins->Add(h1TauTfe2TfeSig[NBin]);
    h1CATfe2Tfe[NBin] = new TH1F(Form("h1CATfe2Tfe[%d]",NBin), Form("h1CATfe2Tfe[%d]",NBin), 50, 0,4);//0,1 bin width
    fBins->Add(h1CATfe2Tfe[NBin]);
    h1TauTfe2Tfe[NBin] = new TH1F(Form("h1TauTfe2Tfe[%d]",NBin), Form("h1TauTfe2Tfe[%d]",NBin), 50, 0,4);
    fBins->Add(h1TauTfe2Tfe[NBin]);
    h1CATfe2TfeCs[NBin] = new TH1F(Form("h1CATfe2TfeCs[%d]",NBin), Form("h1CATfe2TfeCs[%d]",NBin), 50, 0,4);
    fBins->Add(h1CATfe2TfeCs[NBin]);
    h1TauTfe2TfeCs[NBin] = new TH1F(Form("h1TauTfe2TfeCs[%d]",NBin), Form("h1TauTfe2TfeCs[%d]",NBin), 50, 0,4);
    fBins->Add(h1TauTfe2TfeCs[NBin]);
    h1CATfe2TfeCsUg[NBin] = new TH1F(Form("h1CATfe2TfeCsUg[%d]",NBin), Form("h1CATfe2TfeCsUg[%d]",NBin), 50, 0,4);
    fBins->Add(h1CATfe2TfeCsUg[NBin]);
    h1TauTfe2TfeCsUg[NBin] = new TH1F(Form("h1TauTfe2TfeCsUg[%d]",NBin), Form("h1TauTfe2TfeCsUg[%d]",NBin), 50, 0,4);
    fBins->Add(h1TauTfe2TfeCsUg[NBin]);
    //(tau1->3-tau1->2)/tau1->2 in bins of tau1->2
    hdiff1CATfe2TfeSig[NBin] = new TH1F(Form("hdiff1CATfe2TfeSig[%d]",NBin), Form("hdiff1CATfe2TfeSig[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1CATfe2TfeSig[NBin]);
    hdiff1TauTfe2TfeSig[NBin] = new TH1F(Form("hdiff1TauTfe2TfeSig[%d]",NBin), Form("hdiff1TauTfe2TfeSig[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1TauTfe2TfeSig[NBin]);
    hdiff1CATfe2Tfe[NBin] = new TH1F(Form("hdiff1CATfe2Tfe[%d]",NBin), Form("hdiff1CATfe2Tfe[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1CATfe2Tfe[NBin]);
    hdiff1TauTfe2Tfe[NBin] = new TH1F(Form("hdiff1TauTfe2Tfe[%d]",NBin), Form("hdiff1TauTfe2Tfe[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1TauTfe2Tfe[NBin]);
    hdiff1CATfe2TfeCs[NBin] = new TH1F(Form("hdiff1CATfe2TfeCs[%d]",NBin), Form("hdiff1CATfe2TfeCs[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1CATfe2TfeCs[NBin]);
    hdiff1TauTfe2TfeCs[NBin] = new TH1F(Form("hdiff1TauTfe2TfeCs[%d]",NBin), Form("hdiff1TauTfe2TfeCs[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1TauTfe2TfeCs[NBin]);
    hdiff1CATfe2TfeCsUg[NBin] = new TH1F(Form("hdiff1CATfe2TfeCsUg[%d]",NBin), Form("hdiff1CATfe2TfeCsUg[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1CATfe2TfeCsUg[NBin]);
    hdiff1TauTfe2TfeCsUg[NBin] = new TH1F(Form("hdiff1TauTfe2TfeCsUg[%d]",NBin), Form("hdiff1TauTfe2TfeCsUg[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1TauTfe2TfeCsUg[NBin]);
    //tau1->2/tau1->3 in bins of tau1->3
    h1CATfeTfe2Sig[NBin] = new TH1F(Form("h1CATfeTfe2Sig[%d]",NBin), Form("h1CATfeTfe2Sig[%d]",NBin), 50, 0,4);
    fBins->Add(h1CATfeTfe2Sig[NBin]);
    h1TauTfeTfe2Sig[NBin] = new TH1F(Form("h1TauTfeTfe2Sig[%d]",NBin), Form("h1TauTfeTfe2Sig[%d]",NBin), 50, 0,4);
    fBins->Add(h1TauTfeTfe2Sig[NBin]);
    h1CATfeTfe2[NBin] = new TH1F(Form("h1CATfeTfe2[%d]",NBin), Form("h1CATfeTfe2[%d]",NBin), 50, 0,4);
    fBins->Add(h1CATfeTfe2[NBin]);
    h1TauTfeTfe2[NBin] = new TH1F(Form("h1TauTfeTfe2[%d]",NBin), Form("h1TauTfeTfe2[%d]",NBin), 50, 0,4);
    fBins->Add(h1TauTfeTfe2[NBin]);
    h1CATfeTfe2Cs[NBin] = new TH1F(Form("h1CATfeTfe2Cs[%d]",NBin), Form("h1CATfeTfe2Cs[%d]",NBin), 50, 0,4);
    fBins->Add(h1CATfeTfe2Cs[NBin]);
    h1TauTfeTfe2Cs[NBin] = new TH1F(Form("h1TauTfeTfe2Cs[%d]",NBin), Form("h1TauTfeTfe2Cs[%d]",NBin), 50, 0,4);
    fBins->Add(h1TauTfeTfe2Cs[NBin]);
    h1CATfeTfe2CsUg[NBin] = new TH1F(Form("h1CATfeTfe2CsUg[%d]",NBin), Form("h1CATfeTfe2CsUg[%d]",NBin), 50, 0,4);
    fBins->Add(h1CATfeTfe2CsUg[NBin]);
    h1TauTfeTfe2CsUg[NBin] = new TH1F(Form("h1TauTfeTfe2CsUg[%d]",NBin), Form("h1TauTfeTfe2CsUg[%d]",NBin), 50, 0,4);
    fBins->Add(h1TauTfeTfe2CsUg[NBin]);
    //(tau1->2-tau1->3)/tau1->3 in bins of tau1->3
    hdiff1CATfeTfe2Sig[NBin] = new TH1F(Form("hdiff1CATfeTfe2Sig[%d]",NBin), Form("hdiff1CATfeTfe2Sig[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1CATfeTfe2Sig[NBin]);
    hdiff1TauTfeTfe2Sig[NBin] = new TH1F(Form("hdiff1TauTfeTfe2Sig[%d]",NBin), Form("hdiff1TauTfeTfe2Sig[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1TauTfeTfe2Sig[NBin]);
    hdiff1CATfeTfe2[NBin] = new TH1F(Form("hdiff1CATfeTfe2[%d]",NBin), Form("hdiff1CATfeTfe2[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1CATfeTfe2[NBin]);
    hdiff1TauTfeTfe2[NBin] = new TH1F(Form("hdiff1TauTfeTfe2[%d]",NBin), Form("hdiff1TauTfeTfe2[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1TauTfeTfe2[NBin]);
    hdiff1CATfeTfe2Cs[NBin] = new TH1F(Form("hdiff1CATfeTfe2Cs[%d]",NBin), Form("hdiff1CATfeTfe2Cs[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1CATfeTfe2Cs[NBin]);
    hdiff1TauTfeTfe2Cs[NBin] = new TH1F(Form("hdiff1TauTfeTfe2Cs[%d]",NBin), Form("hdiff1TauTfeTfe2Cs[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1TauTfeTfe2Cs[NBin]);
    hdiff1CATfeTfe2CsUg[NBin] = new TH1F(Form("hdiff1CATfeTfe2CsUg[%d]",NBin), Form("hdiff1CATfeTfe2CsUg[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1CATfeTfe2CsUg[NBin]);
    hdiff1TauTfeTfe2CsUg[NBin] = new TH1F(Form("hdiff1TauTfeTfe2CsUg[%d]",NBin), Form("hdiff1TauTfeTfe2CsUg[%d]",NBin), 50, -2,2);
    fBins->Add(hdiff1TauTfeTfe2CsUg[NBin]);
   }

   //Jet substructure -- loop for 0th and first
   for (int nsplit = 0; nsplit < 2; nsplit++){
    //fill 0 with 1st split, 1 with 2nd and 3 with all.
    h1ZgSigCA[nsplit] = new TH1F(Form("h1ZgSigCA[%d]",nsplit), Form("h1ZgSigCA[%d]",nsplit),30,0.,1);
    fSplitList->Add(h1ZgSigCA[nsplit]);
    h1ZgSigTau[nsplit] = new TH1F(Form("h1ZgSigTau[%d]", nsplit), Form("h1ZgSigTau[%d]",nsplit),30,0.,1);
    fSplitList->Add(h1ZgSigTau[nsplit]);
    //angles of the algorithm
    h1RgCASig[nsplit] = new TH1F(Form("h1RgCASig[%d]", nsplit),Form("h1RgCASig[%d]", nsplit),15,0.,0.5);
    fSplitList->Add(h1RgCASig[nsplit]);
    h1RgTauSig[nsplit] = new TH1F(Form("h1RgTauSig[%d]", nsplit),Form("h1RgTauSig[%d]", nsplit),15,0.,0.5);
    fSplitList->Add(h1RgTauSig[nsplit]);
    //energy of sum of the two branches at each node
    h1EradSigCA[nsplit] = new TH1F(Form("h1EradSigCA[%d]", nsplit), Form("h1EradSigCA[%d]", nsplit),50,0,1000);
    fSplitList->Add(h1EradSigCA[nsplit]);
    h1EradSigTau[nsplit] = new TH1F(Form("h1EradSigTau[%d]", nsplit),Form("h1EradSigTau[%d]", nsplit),50,0,1000);
    fSplitList->Add(h1EradSigTau[nsplit]);
    //Log(1/angle) in the algorithm
    h1LogDr12SigCA[nsplit] = new TH1F(Form("h1LogDr12SigCA[%d]", nsplit),Form("h1LogDr12SigCA[%d]", nsplit),25,-5,10);
    fSplitList->Add(h1LogDr12SigCA[nsplit]);
    h1LogDr12SigTau[nsplit] = new TH1F(Form("h1LogDr12SigTau[%d]", nsplit), Form("h1LogDr12SigTau[%d]", nsplit),25,-5,10);
    fSplitList->Add(h1LogDr12SigTau[nsplit]);
    // Log(z*angle) in the algorithm
    h1LogZthetaSigCA[nsplit] = new TH1F(Form("h1LogZthetaSigCA[%d]", nsplit),Form("h1LogZthetaSigCA[%d]", nsplit),25,-10,10);
    fSplitList->Add(h1LogZthetaSigCA[nsplit]);
    h1LogZthetaSigTau[nsplit] = new TH1F(Form("h1LogZthetaSigTau[%d]", nsplit), Form("h1LogZthetaSigTau[%d]", nsplit),25,-10,10);
    fSplitList->Add(h1LogZthetaSigTau[nsplit]);
    // formation time
    h1TfSigTau[nsplit] = new TH1F(Form("h1TfSigTau[%d]", nsplit), Form("h1TfSigTau[%d]", nsplit),50,0,10);
    fSplitList->Add(h1TfSigTau[nsplit]);
    h1TfSigCA[nsplit] = new TH1F(Form("h1TfSigCA[%d]", nsplit), Form("h1TfSigCA[%d]", nsplit),50,0,10);
    fSplitList->Add(h1TfSigCA[nsplit]);
   }


   bool filled = false;
   int jet = 0;
   Long64_t nbytes = 0, nb = 0;
   if(reduced==true){
    std::cout<< " cutting tree "<<std::endl;
   }
   //std::cout << "before first loop " << std::endl;
   double prev_evWeight = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {//nentries
    //std::cout << "+++++++++++++++++++++++++jentry: " << jentry << std::endl;
    hEntries->Fill(1);
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
    //if (Cut(ientry) < 0) continue;
    double evWeight = eventWeight->at(0);
    //std::cout << "+++++++++++++++++++++++ ientry: " << ientry << std::endl;
    //std::cout << "event weight in first loop: " << evWeight << std::endl;
    double maxPt = 0;
    for(int i = 0; i<sigJetPt->size(); ++i) {
      //std::cout<<"je pt "<<sigJetPt->at(i)<<" eta "<<sigJetEta->at(i)<<" phi "<<sigJetPhi->at(i)<<std::endl;           

      double sigPt = sigJetPt->at(i);// is this 0 or i ? for the first jet in each evet
      double sigPhi = sigJetPhi->at(i);
      double sigEta = sigJetEta->at(i);

      if(sigPt>maxPt){
        maxPt = sigPt;
        //std::cout<<"ientry: "<<ientry <<" i: "<<i<<" maxPt "<<maxPt<<std::endl;
      } 
      //double sigPt = sigPtVec[0];
      if (i > 0) continue;//for lead jet only
      if(sigPt<60. || abs(sigJetEta->at(i)) > (0.9 - R)) continue;//thermal events are set from -0.8 to 0.8 - merged files from 70 to 110 
      jet = jet + 1; 
      //std::cout<<"jet "<<i<<std::endl;
      //vectors - formation time
      //sigJets
      std::vector<double> sigTfs = sigJetRecur_tf->at(i);//Johanna 
      std::vector<double> sigTfsTau = sigJetRecurTau_tf->at(i);//Johanna
      std::vector<double> sigTfsZ = sigJetRecurZcut_tf->at(i);//Johanna
      std::vector<double> sigTfsZTau = sigJetRecurTauZcut_tf->at(i);//Johanna
      std::vector<double> sigTfes = sigJetRecur_tfe->at(i);//Morgan
      std::vector<double> sigTfesTau = sigJetRecurTau_tfe->at(i);//Morgan
      std::vector<double> sigTfesZ = sigJetRecurZcut_tfe->at(i);//Morgan
      std::vector<double> sigTfesZTau = sigJetRecurTauZcut_tfe->at(i);//Morgan

      //vectors - momentum fraction
      std::vector<double> sigZgsCA = sigJetRecur_z->at(i);
      std::vector<double> sigZgsTau = sigJetRecurTau_z->at(i);
      std::vector<double> sigZgsCAZ = sigJetRecurZcut_z->at(i);
      std::vector<double> sigZgsTauZ = sigJetRecurTauZcut_z->at(i);
      // Rg
      std::vector<double> sigRgsCA = sigJetRecur_dr12->at(i);
      std::vector<double> sigRgsTau = sigJetRecurTau_dr12->at(i);
      std::vector<double> sigRgsCAZ = sigJetRecurZcut_dr12->at(i);//sig after c/a 0.1 grooming
      std::vector<double> sigRgsTauZ = sigJetRecurTauZcut_dr12->at(i);//sig tau c/a 0.1 grooming
      // Substructure 
      std::vector<double> sigErads = sigJetRecurZcut_erad->at(i);      
      std::vector<double> sigEradsTau = sigJetRecurTauZcut_erad->at(i);
      std::vector<double> sigLogdr12 = sigJetRecurZcut_logdr12->at(i);//log(1/dr)
      std::vector<double> sigLogdr12Tau = sigJetRecurTauZcut_logdr12->at(i);
      std::vector<double> sigLogZtheta = sigJetRecurZcut_logztheta->at(i);        
      std::vector<double> sigLogZthetaTau = sigJetRecurTauZcut_logztheta->at(i);

      std::vector<double> sigKtsTauZ = sigJetRecurTauZcut_kt->at(i);
      std::vector<double> sigKtsCAZ = sigJetRecurTauZcut_kt->at(i);
      std::vector<double> sigKtsTau = sigJetRecurTau_kt->at(i);
      std::vector<double> sigKtsCA = sigJetRecurTau_kt->at(i);

      //droppedBranchInformation
      std::vector<double> brDropDrCA = sigJetRecurZcut_droppedBFS->at(i);// yields dr for every emission before the zg threshold passed (before 1st real split)
      std::vector<double> brDropDrTau = sigJetRecurTauZcut_droppedBFS->at(i);// yields dr for every emission before the zg threshold passed (before 1st real split)
      std::vector<double> brDropLogDrCA = sigJetRecurZcut_droppedLog1drBFS->at(i);// yields log(1/dr) for every emission before the zg threshold passed (before 1st real split)
      std::vector<double> brDropLogDrTau = sigJetRecurTauZcut_droppedLog1drBFS->at(i);// yields log(1/dr) for every emission before the zg threshold passed (before 1st real split)
      
      std::vector<double> brDropPtCA = sigJetRecurZcut_droppedPts->at(i);// yields pt for every emission before the zg threshold passed (before 1st real split)
      std::vector<double> brDropPtTau = sigJetRecurTauZcut_droppedPts->at(i);// yields pt for every emission before the zg threshold passed (before 1st real split)
      std::vector<double> brDropKtCA = sigJetRecurZcut_droppedKts->at(i);// yields kt = min(PT1, PT2)*DeltaR for every emission before the zg threshold passed (before 1st real split)
      std::vector<double> brDropKtTau = sigJetRecurTauZcut_droppedKts->at(i);// yields kt = min(PT1, PT2)*DeltaR for every emission before the zg threshold passed (before 1st real split)
      

      //csFullJets
      std::vector<double> csTfes = csFullJetRecur_tfe->at(i);
      std::vector<double> csTfs = csFullJetRecur_tf->at(i);//Johanna
      std::vector<double> csTauTfes = csFullJetRecurTau_tfe->at(i);
      std::vector<double> csTauTfs = csFullJetRecurTau_tf->at(i);//Johanna

      std::vector<double> csRgsCA = csFullJetRecur_dr12->at(i);
      std::vector<double> csRgsTau = csFullJetRecurTau_dr12->at(i);
      std::vector<double> csKtsCA = csFullJetRecur_kt->at(i);
      std::vector<double> csKtsTau = csFullJetRecurTau_kt->at(i);

      std::vector<double> csTfesZ = csFullJetRecurZcut_tfe->at(i);
      std::vector<double> csTfsZ = csFullJetRecurZcut_tf->at(i);//Johanna
      std::vector<double> csTauTfesZ = csFullJetRecurTauZcut_tfe->at(i);
      std::vector<double> csTauTfsZ = csFullJetRecurTauZcut_tf->at(i);//Johanna



      //analysis
      h1JPtW->Fill(sigPt);// add jet norm im plotting - from pt spectrun
      h1JPt->Fill(sigPt);
      //first splits - matched to signal jets
      if(generator == "PYTHIA"){
        int sigJetID = sigJetRecur_partonMatchID->at(i);
        if(sigJetID == 3) sigJetID = 1;// because ID is 0 or 3 for -1 status codes but splits are stored at 0 or 1
        double drparton = sigJetRecur_partonMatchDr->at(i);
        double tfeSplit = -1;
        double zgSplit = -1; 
        double drSplit = -1;
        double ktSplit = -1;
        if(drparton < 0.05) { //only consider partons close to jet - this we can adjust to get proper stats and/or better correlations
          if(tfesplit->size()>sigJetID){
            tfeSplit = tfesplit->at(sigJetID);//all sigJetID
            zgSplit = zgsplit->at(sigJetID);
            drSplit = drsplit->at(sigJetID);
            ktSplit = ktsplit->at(sigJetID);
            //std::cout<<"tfesplit diff: "<<tfeSplit-sigTfs[0]<<" ptDiff "<<partonWd1d2pt->at(sigJetID)-sigJetPt->at(i)<<std::endl;
            //std::cout<<" drparton "<<drparton<<" zgparton-zgJetCA "<<zgsplit->at(sigJetID)-sigZgsCA[0]<<std::endl;
            //std::cout<<"pt "<<partonsPt->at(sigJetID)<<" eta "<<partonsEta->at(sigJetID)<<" phi "<<partonsPhi->at(sigJetID)<<std::endl;           
            if(sigTfs.size()>0 && tfeSplit > 0) {
              h2PartonSigTf->Fill(tfeSplit,sigTfs[0]);
            }
            if(sigTfsTau.size()>0 && tfeSplit > 0) {
              h2PartonSigTfTau->Fill(tfeSplit, sigTfsTau[0]);
            }
            if(sigTfsZ.size()>0 && tfeSplit > 0) {
              h2PartonSigTfZcut->Fill(tfeSplit, sigTfsZ[0]);
              h1LogTfSigToPartonZ->Fill(log(sigTfsZ[0]));// sd grooming
  
            }
            if(sigTfsZTau.size()>0 && tfeSplit > 0) {
              h2PartonSigTfZcutTau->Fill(tfeSplit, sigTfsZTau[0]);
              h1LogTfSigToPartonZTau->Fill(log(sigTfsZTau[0]));// sd grooming
  
            }
  
            if(sigTfesZ.size()>0 && tfeSplit > 0) {
              h2PartonSigTfeZcut->Fill(tfeSplit, sigTfesZ[0]);
            }
            if(sigTfesZTau.size()>0 && tfeSplit > 0) {
              h2PartonSigTfeZcutTau->Fill(tfeSplit, sigTfesZTau[0]);
            }
            // look at first parton LP - in production on stbc for R = 0.4
            h2PartonLogTfKt->Fill(log(tfeSplit), log(ktSplit));
            h2PartonLogDrKt->Fill(log(1/drSplit), log(ktSplit));
            h2PartonLundPlaneNT->Fill(log(R/drSplit), log(zgSplit*drSplit/R));
            if(zgSplit >= 0.1){
              h2PartonLogTfKtZcut->Fill(log(tfeSplit), log(ktSplit));
              h2PartonLogDrKtZcut->Fill(log(1/drSplit), log(ktSplit));
            }
          }// if enough entries in split vector
        }else{
         // std::cout<<"-----------------------------------rejected parton to jetdistance drparton:: "<<drparton<<std::endl;
        }

      }
 
      // Check unmatched info
      // Substructure observables
      // Zg
      if(sigZgsCAZ.size()>0) h1ZgSigCA[0]->Fill(sigZgsCAZ[0]);
      if(sigZgsCAZ.size()>1) h1ZgSigCA[1]->Fill(sigZgsCAZ[1]);
      
      if(sigZgsTauZ.size()>0) h1ZgSigTau[0]->Fill(sigZgsTauZ[0]);
      if(sigZgsTauZ.size()>1) h1ZgSigTau[1]->Fill(sigZgsTauZ[1]);
            
      if(sigRgsCAZ.size()>0) h1RgCASig[0]->Fill(sigRgsCAZ[0]);
      if(sigRgsCAZ.size()>1) h1RgCASig[1]->Fill(sigRgsCAZ[1]);
      if(sigRgsTauZ.size()>0) h1RgTauSig[0]->Fill(sigRgsTauZ[0]);
      if(sigRgsTauZ.size()>0) h1RgTauSig[1]->Fill(sigRgsTauZ[1]);
      // erad 
      if(sigErads.size()>0) h1EradSigCA[0]->Fill(sigErads[0]);
      if(sigErads.size()>1) h1EradSigCA[1]->Fill(sigErads[1]);
      if(sigEradsTau.size()>0) h1EradSigTau[0]->Fill(sigEradsTau[0]);
      if(sigEradsTau.size()>1) h1EradSigTau[1]->Fill(sigEradsTau[1]);
      // log(1/angle)
      if(sigLogdr12.size()>0) h1LogDr12SigCA[0]->Fill(sigLogdr12[0]);
      if(sigLogdr12.size()>1) h1LogDr12SigCA[1]->Fill(sigLogdr12[1]);  
      if(sigLogdr12Tau.size()>0) h1LogDr12SigTau[0]->Fill(sigLogdr12Tau[0]);
      if(sigLogdr12Tau.size()>1) h1LogDr12SigTau[1]->Fill(sigLogdr12Tau[1]);

      // log(z*angle)
      if(sigLogZtheta.size()>0) h1LogZthetaSigCA[0]->Fill(sigLogZtheta[0]);
      if(sigLogZtheta.size()>1) h1LogZthetaSigCA[1]->Fill(sigLogZtheta[1]);        
      if(sigLogZthetaTau.size()>0) h1LogZthetaSigTau[0]->Fill(sigLogZthetaTau[0]);
      if(sigLogZthetaTau.size()>1) h1LogZthetaSigTau[1]->Fill(sigLogZthetaTau[1]);
 
      // check (unmatched) sigJets correlations tf to tfe
      if(sigTfes.size()>0 && sigTfs.size() > 0) h2SigTfTfe->Fill(sigTfs[0],sigTfes[0]);
      if(sigTfesTau.size()>0 && sigTfsTau.size() > 0) h2SigTfTfeTau->Fill(sigTfsTau[0], sigTfesTau[0]);
      if(sigTfesZ.size()>0 && sigTfsZ.size() > 0) h2SigTfTfeZcut->Fill(sigTfsZ[0], sigTfesZ[0]);
      if(sigTfesZTau.size()>0 && sigTfsZTau.size() > 0) h2SigTfTfeZcutTau->Fill(sigTfsZTau[0], sigTfesZTau[0]);
      //Log tau_form
      if(sigTfs.size()>0) h1LogTfSig->Fill(log(sigTfs[0]));
      if(sigTfsZ.size()>0) h1LogTfSigZ->Fill(log(sigTfsZ[0]));// sd grooming
      if(sigTfsTau.size()>0) h1LogTfSigTau->Fill(log(sigTfsTau[0]));
      if(sigTfsZTau.size()>0) h1LogTfSigZTau->Fill(log(sigTfsZTau[0]));// sd grooming
      
      // formation time unmatched
      if(sigTfsZ.size()>0){
        h1SigNCaTfZcut->Fill(sigTfsZ.size());
        h1TfSigCA[0]->Fill(sigTfsZ[0]);
      }
      if(sigTfsZ.size()>1) h1TfSigCA[1]->Fill(sigTfsZ[1]);
      if(sigTfsZTau.size()>0){
        h1SigNTauTfZcut->Fill(sigTfsZTau.size());
        h1TfSigTau[0]->Fill(sigTfsZTau[0]);
      } 
      if(sigTfsZTau.size()>1) h1TfSigTau[1]->Fill(sigTfsZTau[1]);

      // Lund Plane Novel Tools Fig3 
      for(int i = 0; i<sigRgsCA.size(); i++){
        h2CALundPlaneNT->Fill(log(R/sigRgsCA[i]), log(sigZgsCA[i]*sigRgsCA[i]/R));//all emissions
      }
      for(int i = 0; i<sigRgsTau.size(); i++){
        h2TauLundPlaneNT->Fill(log(R/sigRgsTau[i]), log(sigZgsTau[i]*sigRgsTau[i]/R));//all emissions
      }
      // checking zg's and rg's correlations from tau to C/A
      //2D corr 
      if(sigRgsTauZ.size()>0 && sigRgsCAZ.size()>0) h2SigRgTauRgZcut->Fill(sigRgsCAZ[0], sigRgsTauZ[0]);
      if(sigZgsTauZ.size()>0 && sigZgsCAZ.size()>0) h2SigZgTauZgZcut->Fill(sigZgsCAZ[0], sigZgsTauZ[0]);

      //ungroomed correlation theta zg
      if(sigRgsCA.size()>0 && sigZgsCA.size()>0) h2ThetaZgCA->Fill(sigZgsCA[0], sigRgsCA[0]/R);
      if(sigRgsTau.size()>0 && sigZgsTau.size()>0) h2ThetaZgTau->Fill(sigZgsTau[0], sigRgsTau[0]/R);

      //groomed correlations theta tau
      if(sigZgsCAZ.size()>0 && sigRgsCAZ.size()>0) h2ThetaZgCAZcut->Fill(sigZgsCAZ[0], sigRgsCAZ[0]/R);
      if(sigZgsTauZ.size()>0 && sigRgsTauZ.size()>0) h2ThetaZgTauZcut->Fill(sigZgsTauZ[0], sigRgsTauZ[0]/R);

      //check what is dropped before first accepted emission
      if(brDropDrCA.size()>0){//CA loop
        int n = 0;
        for(int j = 0; j< brDropDrCA.size(); j++){
          if(brDropDrCA[j] < 99){
            n = n+1;
            h2DroppedDrZcut->Fill(n, brDropDrCA[j]);
            h2DroppedPtZcut->Fill(n, brDropPtCA[j]);
            h2DroppedKtZcut->Fill(n, brDropKtCA[j]);
            h2DroppedLPZcut->Fill(brDropLogDrCA[j], log(brDropKtCA[j]));
          }
        }
      }

      if(brDropDrTau.size()>0){//Tau loop
        int n = 0;
        for(int j = 0; j< brDropDrTau.size(); j++){
          if(brDropDrTau[j] < 99){
            n = n+1;
            h2DroppedDrTauZcut->Fill(n, brDropDrTau[j]);
            h2DroppedPtTauZcut->Fill(n, brDropPtTau[j]);
            h2DroppedKtTauZcut->Fill(n, brDropKtTau[j]);
            h2DroppedLPTauZcut->Fill(brDropLogDrTau[j], log(brDropKtTau[j]));
          }
        }
      }
     //LP's for 1 jet
     //only one jet - all splits - ungroomed
     if( jet == 11){
       for(int i = 0; i<sigRgsCA.size(); i++){
        //std::cout<<"============================================================================= pt "<<sigPt<<" sigPhi "<<sigPhi <<" eta "<<sigEta<<std::endl;
         h2CALundPlaneSigSingle->Fill(log(R/sigRgsCA[i]), log(sigKtsCA[i]), i+1);//fill with weight
         filled = true;
       }
       for(int i = 0; i<sigRgsTau.size(); i++){
        h2TauLundPlaneSigSingle->Fill(log(R/sigRgsTau[i]), log(sigKtsTau[i]), i+1);
       }
       for(int i = 0; i<sigTfs.size(); i++){
         h2CALundPlaneSigTfeSingle->Fill(log(sigTfs[i]), log(sigKtsCA[i]), i+1);
       }
       for(int i = 0; i<sigTfsTau.size(); i++){
         h2TauLundPlaneSigTfeSingle->Fill(log(sigTfsTau[i]), log(sigKtsTau[i]), i+1);
       }
     }
     //all jets splits in separet canvas
      if (sigKtsCAZ.size() > 0){
        for(int i = 0; i<sigRgsCAZ.size(); i++){
          if(sigRgsCAZ.size()>0) h3CALundPlaneSigSingleZ[0]->Fill(log(R/sigRgsCAZ[0]), log(sigKtsCAZ[0]), 0);
          if(sigRgsCAZ.size()>1) h3CALundPlaneSigSingleZ[1]->Fill(log(R/sigRgsCAZ[1]), log(sigKtsCAZ[1]), 1);
          if(sigRgsCAZ.size()>2) h3CALundPlaneSigSingleZ[2]->Fill(log(R/sigRgsCAZ[2]), log(sigKtsCAZ[2]), 2);
          if(sigRgsCAZ.size()>3) h3CALundPlaneSigSingleZ[3]->Fill(log(R/sigRgsCAZ[3]), log(sigKtsCAZ[3]), 3);
          h3CALundPlaneSigSingleInclZ->Fill(log(R/sigRgsCAZ[i]), log(sigKtsCAZ[i]), i);//store number of emission for every jet separately
        }
        for(int i = 0; i<sigTfsZ.size(); i++){
            if(sigTfsZ.size()>0) h3CALundPlaneSigTfeSingleZ[0]->Fill(log(sigTfsZ[0]), log(sigKtsCAZ[0]), 0);
            if(sigTfsZ.size()>1) h3CALundPlaneSigTfeSingleZ[1]->Fill(log(sigTfsZ[1]), log(sigKtsCAZ[1]), 1);
            if(sigTfsZ.size()>2) h3CALundPlaneSigTfeSingleZ[2]->Fill(log(sigTfsZ[2]), log(sigKtsCAZ[2]), 2);
            if(sigTfsZ.size()>3) h3CALundPlaneSigTfeSingleZ[3]->Fill(log(sigTfsZ[3]), log(sigKtsCAZ[3]), 3);
            h3CALundPlaneSigTfeSingleInclZ->Fill(log(sigTfsZ[i]), log(sigKtsCAZ[i]), i);
          }
      }
      if (sigKtsTauZ.size() > 0){
          for(int i = 0; i<sigRgsTauZ.size(); i++){
            if(sigRgsTauZ.size()>0) h3TauLundPlaneSigSingleZ[0]->Fill(log(R/sigRgsTauZ[0]), log(sigKtsTauZ[0]), 0);
            if(sigRgsTauZ.size()>1) h3TauLundPlaneSigSingleZ[1]->Fill(log(R/sigRgsTauZ[1]), log(sigKtsTauZ[1]), 1);
            if(sigRgsTauZ.size()>2) h3TauLundPlaneSigSingleZ[2]->Fill(log(R/sigRgsTauZ[2]), log(sigKtsTauZ[2]), 2);
            if(sigRgsTauZ.size()>3) h3TauLundPlaneSigSingleZ[3]->Fill(log(R/sigRgsTauZ[3]), log(sigKtsTauZ[3]), 3);
            h3TauLundPlaneSigSingleInclZ->Fill(log(R/sigRgsTauZ[i]), log(sigKtsTauZ[i]), i);
          }
          for(int i = 0; i<sigTfsZTau.size(); i++){
            if(sigTfsZTau.size()>0) h3TauLundPlaneSigTfeSingleZ[0]->Fill(log(sigTfsZTau[0]), log(sigKtsTauZ[0]), 0);
            if(sigTfsZTau.size()>1) h3TauLundPlaneSigTfeSingleZ[1]->Fill(log(sigTfsZTau[1]), log(sigKtsTauZ[1]), 1);
            if(sigTfsZTau.size()>2) h3TauLundPlaneSigTfeSingleZ[2]->Fill(log(sigTfsZTau[2]), log(sigKtsTauZ[2]), 2);
            if(sigTfsZTau.size()>3) h3TauLundPlaneSigTfeSingleZ[3]->Fill(log(sigTfsZTau[3]), log(sigKtsTauZ[3]), 3);
            h3TauLundPlaneSigTfeSingleInclZ->Fill(log(sigTfsZTau[i]), log(sigKtsTauZ[i]), i);
          }
      }

      // Here without z cut
      if (sigKtsCA.size() > 0){
        for(int i = 0; i<sigRgsCA.size(); i++){
          if(sigRgsCA.size()>0) h3CALundPlaneSigSingle[0]->Fill(log(R/sigRgsCA[0]), log(sigKtsCA[0]), 0);
          if(sigRgsCA.size()>1) h3CALundPlaneSigSingle[1]->Fill(log(R/sigRgsCA[1]), log(sigKtsCA[1]), 1);
          if(sigRgsCA.size()>2) h3CALundPlaneSigSingle[2]->Fill(log(R/sigRgsCA[2]), log(sigKtsCA[2]), 2);
          if(sigRgsCA.size()>3) h3CALundPlaneSigSingle[3]->Fill(log(R/sigRgsCA[3]), log(sigKtsCA[3]), 3);
          h2CALundPlaneSig->Fill(log(R/sigRgsCA[i]), log(sigKtsCA[i]));
          h3CALundPlaneSigSingleIncl->Fill(log(R/sigRgsCA[i]), log(sigKtsCA[i]), i);//store number of emission for every jet separately
        }
        for(int i = 0; i<sigTfs.size(); i++){
          if(sigTfs.size()>0) h3CALundPlaneSigTfeSingle[0]->Fill(log(sigTfs[0]), log(sigKtsCA[0]), 0);
          if(sigTfs.size()>1) h3CALundPlaneSigTfeSingle[1]->Fill(log(sigTfs[1]), log(sigKtsCA[1]), 1);
          if(sigTfs.size()>2) h3CALundPlaneSigTfeSingle[2]->Fill(log(sigTfs[2]), log(sigKtsCA[2]), 2);
          if(sigTfs.size()>3) h3CALundPlaneSigTfeSingle[3]->Fill(log(sigTfs[3]), log(sigKtsCA[3]), 3);
          h2CALundPlaneSigTfe->Fill(log(sigTfs[i]), log(sigKtsCA[i]));
          h3CALundPlaneSigTfeSingleIncl->Fill(log(sigTfs[i]), log(sigKtsCA[i]), i);
        }
      }
      if (sigKtsTau.size() > 0){
        for(int i = 0; i<sigRgsTau.size(); i++){
          if(sigRgsTau.size()>0) h3TauLundPlaneSigSingle[0]->Fill(log(R/sigRgsTau[0]), log(sigKtsTau[0]), 0);
          if(sigRgsTau.size()>1) h3TauLundPlaneSigSingle[1]->Fill(log(R/sigRgsTau[1]), log(sigKtsTau[1]), 1);
          if(sigRgsTau.size()>2) h3TauLundPlaneSigSingle[2]->Fill(log(R/sigRgsTau[2]), log(sigKtsTau[2]), 2);
          if(sigRgsTau.size()>3) h3TauLundPlaneSigSingle[3]->Fill(log(R/sigRgsTau[3]), log(sigKtsTau[3]), 3);
          h2TauLundPlaneSig->Fill(log(R/sigRgsTau[i]), log(sigKtsTau[i]));
          h3TauLundPlaneSigSingleIncl->Fill(log(R/sigRgsTau[i]), log(sigKtsTau[i]), i);
        }
        for(int i = 0; i<sigTfsTau.size(); i++){
          if(sigTfsTau.size()>0) h3TauLundPlaneSigTfeSingle[0]->Fill(log(sigTfsTau[0]), log(sigKtsTau[0]), 0);
          if(sigTfsTau.size()>1) h3TauLundPlaneSigTfeSingle[1]->Fill(log(sigTfsTau[1]), log(sigKtsTau[1]), 1);
          if(sigTfsTau.size()>2) h3TauLundPlaneSigTfeSingle[2]->Fill(log(sigTfsTau[2]), log(sigKtsTau[2]), 2);
          if(sigTfsTau.size()>3) h3TauLundPlaneSigTfeSingle[3]->Fill(log(sigTfsTau[3]), log(sigKtsTau[3]), 3);
          h2TauLundPlaneSigTfe->Fill(log(sigTfsTau[i]), log(sigKtsTau[i]));
          h3TauLundPlaneSigTfeSingleIncl->Fill(log(sigTfsTau[i]), log(sigKtsTau[i]), i);
        }
      }

      ///////////////////////////////////////////////////////
      ///               CS studies                        ///   
      ///                                                ///
      //////////////////////////////////////////////////////
      //first splits - matched to CS 
      if(generator == "PYTHIA"){
        int csFullID = csFullJetRecur_partonMatchID->at(i);  
        if(csFullID == 3) csFullID = 1;// because ID is 0 or 3 for -1 status codes but splits are stored at 0 or 1
        double drpartonCS = csFullJetRecur_partonMatchDr->at(i);
        double tfeSplit = -1;
        double drSplit = -1;
        double ktSplit = -1;
        if(drpartonCS < 0.05) {
          if(tfesplit->size()>0 && tfesplit->size()>csFullID){
            tfeSplit = tfesplit->at(csFullID);
            drSplit = drsplit->at(csFullID);
            ktSplit = ktsplit->at(csFullID);
            if(csTfes.size() > 0) {
              h2CsPartonTfe->Fill(tfeSplit,csTfs[0]);
            }
            if(csTauTfs.size() > 0) {
              h2CsPartonTfeTau->Fill(tfeSplit,csTauTfs[0]);
            }
            if(sigTfsZ.size()>0 && csTfs.size() > 0) {
              h2CsSigTfeZcut->Fill(csTfs[0], sigTfsZ[0]);
            }
            if(sigTfsZTau.size()>0 && csTauTfs.size() > 0) {
              h2CsSigTfeZcutTau->Fill(csTauTfes[0], sigTfsZTau[0]);
            }

          }
        }else{
        //std::cout<<"Missing/Unmatched Jet nr:"<<i<<" with partonMatchID: "<<sigJetID<<" csFullpartonmatchID: "<<csFullID<<" drParton: "<<drpartonCS<<std::endl;
        }
      }
      if (csKtsCA.size() > 0){
        for(int i = 0; i<csRgsCA.size(); i++){
          h2CALundPlaneCs->Fill(log(R/csRgsCA[i]), log(csKtsCA[i]));
        }
        for(int i = 0; i<csTfs.size(); i++){
          h2CALundPlaneCsTfe->Fill(log(csTfs[i]), log(csKtsCA[i]));
        }
      }
      if (csKtsTau.size() > 0){
        for(int i = 0; i<csRgsTau.size(); i++){
          h2TauLundPlaneCs->Fill(log(R/csRgsTau[i]), log(csKtsTau[i]));
        }
        for(int i = 0; i<csTauTfs.size(); i++){
          h2TauLundPlaneCsTfe->Fill(log(csTauTfs[i]), log(csKtsTau[i]));
        }
      }

  }// end of jet entry loop (jentry)


  /*
      std::vector<double> csTauPtVec = csFullJetRecurTau_jetpt->at(i);
      std::vector<double> csCAPtVec = csFullJetRecur_jetpt->at(i);
      std::vector<double> sigTauPtVec = sigJetRecurTau_jetpt->at(i);
      std::vector<double> sigCAPtVec = sigJetRecur_jetpt->at(i);
      std::vector<double> sigTauEtaVec = sigJetRecurTau_jeteta->at(i);
      std::vector<double> sigCAEtaVec = sigJetRecur_jeteta->at(i);
      std::vector<double> csTauEtaVec = csFullJetRecurTau_jeteta->at(i);
      std::vector<double> csCAEtaVec = csFullJetRecur_jeteta->at(i);

      if(sigCAPtVec.size()>0){
        h1umsigCAPt->Fill(sigCAPtVec[0]);
      }
      if(sigCAEtaVec.size()>0){
        h1umsigCAEta->Fill(sigCAEtaVec[0]);      
      }
      if(sigTauPtVec.size()>0){
        h1umsigTauPt->Fill(sigTauPtVec[0]);
      }
      if(sigTauEtaVec.size()>0){
        h1umsigTauEta->Fill(sigTauEtaVec[0]);
      }
      //Pythia info from leading branch
      jet = jet+1;
      
      //fill unmatched sig histo
      //Nsubjettiness (ca/tau)
      double tauCa21 = sigJetRecurZcut_tau21->at(i);
      double tauTau21 = sigJetRecurTauZcut_tau21->at(i);
      double tauCa32 = sigJetRecurZcut_tau32->at(i);
      double tauTau32 = sigJetRecurTauZcut_tau32->at(i);
      if( tauCa21>-99 && tauTau21>-99) h2CaTauNsubj21->Fill(tauCa21, tauTau21);
      if( tauCa32>-99 && tauTau32>-99) h2CaTauNsubj32->Fill(tauCa32, tauTau32);
      if(sigJetID == csFullID){
        //fill matched sig histo
        if(sigTauEtaVec.size()>0){
          h1matchedsigCAEta->Fill(sigCAEtaVec[0]);
        }
        if(sigTauEtaVec.size()>0){
          h1matchedsigTauEta->Fill(sigTauEtaVec[0]);
        }
        if(sigCAPtVec.size()>0){
          h1matchedsigCAPt->Fill(sigCAPtVec[0]);
          h1Ca->Fill(sigPt-sigCAPtVec[0]);
          h2Ca->Fill(sigPt,sigCAPtVec[0]);
        }
        if(sigTauPtVec.size()>0){
          h1matchedsigTauPt->Fill(sigTauPtVec[0]);
          h1Tau->Fill(sigPt-sigTauPtVec[0]);
          h2Tau->Fill(sigPt,sigTauPtVec[0]);
        }

        if(reduced==true && generator=="CA"){
          if(!(sigTfesZ.size()>0)) continue;
          if(duration == "short" && (sigTfesZ[0] >= 1)) continue;
          if(duration == "supershort" && (sigTfesZ[0] >= 0.5)) continue;
          if(duration == "large" && (sigTfesZ[0] <= 3)) continue;
          if(duration == "intermed" && ((sigTfesZ[0] < 1) || (sigTfesZ[0] > 3))) continue;
        }
        if(reduced==true && generator=="Tau"){
          if(!(sigTfesZTau.size()>0)) continue;
          if(duration == "short" && (sigTfesZTau[0] >= 1)) continue;
          if(duration == "supershort" && (sigTfesZTau[0] >= 0.5)) continue;
          if(duration == "large" && (sigTfesZTau[0] <= 3)) continue;
          if(duration == "intermed" && ((sigTfesZTau[0] < 1) || (sigTfesZTau[0] > 3))) continue;
        }
        if(reduced==true && generator=="Soft"){
          if(!(sigTfesZTau.size()>0)) continue;
          if(duration == "short" && (sigTfsZTau[0] >= 1)) continue;
          if(duration == "supershort" && (sigTfsZTau[0] >= 0.5)) continue;
          if(duration == "large" && (sigTfsZTau[0] <= 3)) continue;
          if(duration == "intermed" && ((sigTfsZTau[0] < 1) || (sigTfsZTau[0] > 3))) continue;
        }
        

        std::vector<double> sigCAPt = sigJetRecurZcut_jetpt->at(i);//sig after c/a 
        if(sigCAPt.size()>0){
          h1JsigCAPtW->Fill(sigCAPt[0], evWeight);
          h1JsigCAPt->Fill(sigCAPt[0]);
        } 

        std::vector<double> sigTauPt = sigJetRecurTauZcut_jetpt->at(i);//sig after tau
        if(sigTauPt.size()>0){
          h1JsigTauPtW->Fill(sigTauPt[0], evWeight);
          h1JsigTauPt->Fill(sigTauPt[0]);
        } 

        if(sigTfes.size()>0 && sigTfesTau.size()>0) h2SigTfeTauTfe->Fill(sigTfes[0], sigTfesTau[0]);
        if(sigTfesZ.size()>0 && sigTfesZTau.size()>0) h2SigTfeTauTfeZcut->Fill(sigTfesZ[0], sigTfesZTau[0]);

        //Number of split times calculated via CA algorithm 
        if(sigTfes.size()>0) h1SigNCaTf->Fill(sigTfes.size());

        //Number of split times calculated via tau algorithm
        if(sigTfesTau.size()>0) h1SigNTauTf->Fill(sigTfesTau.size());

        //std::cout<<"before dropped"<<std::endl;
        std::vector<double> branchDroppedT = sigJetRecurTauZcut_droppedBFS->at(i);// yields dr for every emission before the zg threshold passed (before 1st real split)
        std::vector<double> branchDroppedTTf = sigJetRecurTauZcut_droppedTfBFS->at(i);// yields dr for every emission before the zg threshold passed (before 1st real split)
        std::vector<double> branchDroppedTTfe = sigJetRecurTauZcut_droppedTfeBFS->at(i);// yields dr for every emission before the zg threshold passed (before 1st real split)
        
            dropped.push_back(DeltaR);
            droppedTf.push_back(2./(zg*CurrentJet.perp()*GeVtofm*(DeltaR/r0_)*(DeltaR/r0_)));
            droppedTfe.push_back(1./(2.*zg*(1-zg)*CurrentJet.perp()*GeVtofm*(1-cos(DeltaR/r0_))));
            droppedPt.push_back(CurrentJet.perp());
            droppedKt.push_back(min(PT1, PT2)*DeltaR);
            droppedLog1dr.push_back(log(1/DeltaR));


        if(branchDroppedT.size()>0){
          int n = 0;
          for(int j = 0; j< branchDroppedT.size(); j++){
            if(branchDroppedT[j] < 99){
              n = n+1;
              }
          }
          if(sigTfesZTau.size() > 0 && sigTfsZTau.size() > 0){
            h2DroppedLogTfTauZcut->Fill(n, log(sigTfsZTau[0])); 
            h2DroppedLogTfeTauZcut->Fill(n, log(sigTfesZTau[0])); 
          }
          if(sigTfesZTau.size() > 1 && sigTfsZTau.size() > 0){
            h2DroppedLog2TfTauZcut->Fill(n, log(sigTfsZTau[1])); 
            h2DroppedLog2TfeTauZcut->Fill(n, log(sigTfesZTau[1])); 
          }
        } 
        
        std::vector<double> branchDropped = sigJetRecurZcut_droppedBFS->at(i);// yields dr for every emission before the zg threshold passed (before 1st real split)
        std::vector<double> branchDroppedTf = sigJetRecurZcut_droppedTfBFS->at(i);// yields dr for every emission before the zg threshold passed (before 1st real split)
        std::vector<double> branchDroppedTfe = sigJetRecurZcut_droppedTfeBFS->at(i);// yields dr for every emission before the zg threshold passed (before 1st real split)        if(branchDropped.size()>0){
        if(branchDropped.size()>0){
          int n = 0;
          for(int j = 0; j< branchDropped.size(); j++){
            if(branchDropped[j] < 99){
              n = n+1;
              }
          }
          if(sigRgs.size()>0) h2DroppedDrZcut->Fill(n, sigRgs[0]); 
          if(sigTfesZ.size() > 0 && sigTfsZ.size() > 0){
            h2DroppedLogTfZcut->Fill(n, log(sigTfsZ[0])); 
            h2DroppedLogTfeZcut->Fill(n, log(sigTfesZ[0])); 
          }
          if(sigTfesZ.size() > 1 && sigTfsZ.size() > 1){
            h2DroppedLog2TfZcut->Fill(n, log(sigTfsZ[1])); 
            h2DroppedLog2TfeZcut->Fill(n, log(sigTfesZ[1])); 
          }
        } 

        // Lund Planes ungroomed
        //log (kt) = log(pt*dr/R) vs log(R/dr)
        std::vector<double> ugKtsCA = sigJetRecur_kt->at(i);
        std::vector<double> ugKtsTau = sigJetRecurTau_kt->at(i);

        //Full Cs
        std::vector<double> ugKtsCACs = csFullJetRecur_kt->at(i);
        std::vector<double> ugKtsTauCs = csFullJetRecurTau_kt->at(i);
        std::vector<double> ugRgsCACs = csFullJetRecur_dr12->at(i);
        std::vector<double> ugRgsTauCs = csFullJetRecurTau_dr12->at(i);
        std::vector<double> ugTfesCACs = csFullJetRecur_tfe->at(i);
        std::vector<double> ugTfesTauCs = csFullJetRecurTau_tfe->at(i);

        //log (kt) = log(pt*dr/R) vs log(tfe)   
        std::vector<double> ugTfesCA = sigJetRecur_tfe->at(i);
        std::vector<double> ugTfesTau = sigJetRecurTau_tfe->at(i);




        Double_t kcut = 5;// because log(10)==1
        Double_t LPcut = 2;
        //all jets splits in separet canvas
        //if( finsublled == false){//use nsubjettiness as discriminator to see how they change in terms of hardness
        if (ugKtsCA.size() > 0){
        //  if( ugKtsCA[0] > kcut){
          //if( sqrt(ugKtsCA[0]*ugKtsCA[0] + (R/ugRgsCA[0])*(R/ugRgsCA[0]) ) < LPcut){
          for(int i = 0; i<ugRgsCA.size(); i++){
            if(ugRgsCA.size()<2) h3CALundPlaneSigSingle[0]->Fill(log(R/ugRgsCA[i]), log(ugKtsCA[i]), i);
            if(ugRgsCA.size()<3) h3CALundPlaneSigSingle[1]->Fill(log(R/ugRgsCA[i]), log(ugKtsCA[i]), i);
            if(ugRgsCA.size()<4) h3CALundPlaneSigSingle[2]->Fill(log(R/ugRgsCA[i]), log(ugKtsCA[i]), i);
            if(ugRgsCA.size()<5) h3CALundPlaneSigSingle[3]->Fill(log(R/ugRgsCA[i]), log(ugKtsCA[i]), i);
            if(ugRgsCA.size()<6) h3CALundPlaneSigSingle[4]->Fill(log(R/ugRgsCA[i]), log(ugKtsCA[i]), i);
            if(ugRgsCA.size()<7) h3CALundPlaneSigSingle[5]->Fill(log(R/ugRgsCA[i]), log(ugKtsCA[i]), i);
            if(ugRgsCA.size()<8) h3CALundPlaneSigSingle[6]->Fill(log(R/ugRgsCA[i]), log(ugKtsCA[i]), i);
            h3CALundPlaneSigSingleIncl->Fill(log(R/ugRgsCA[i]), log(ugKtsCA[i]), i);//store number of emission for every jet separately
          }
          for(int i = 0; i<ugTfesCA.size(); i++){
            if(ugTfesCA.size()<2) h3CALundPlaneSigTfeSingle[0]->Fill(log(ugTfesCA[i]), log(ugKtsCA[i]), i);
            if(ugTfesCA.size()<3) h3CALundPlaneSigTfeSingle[1]->Fill(log(ugTfesCA[i]), log(ugKtsCA[i]), i);
            if(ugTfesCA.size()<4) h3CALundPlaneSigTfeSingle[2]->Fill(log(ugTfesCA[i]), log(ugKtsCA[i]), i);
            if(ugTfesCA.size()<5) h3CALundPlaneSigTfeSingle[3]->Fill(log(ugTfesCA[i]), log(ugKtsCA[i]), i);
            if(ugTfesCA.size()<6) h3CALundPlaneSigTfeSingle[4]->Fill(log(ugTfesCA[i]), log(ugKtsCA[i]), i);
            if(ugTfesCA.size()<7) h3CALundPlaneSigTfeSingle[5]->Fill(log(ugTfesCA[i]), log(ugKtsCA[i]), i);
            if(ugTfesCA.size()<8) h3CALundPlaneSigTfeSingle[6]->Fill(log(ugTfesCA[i]), log(ugKtsCA[i]), i);
            h3CALundPlaneSigTfeSingleIncl->Fill(log(ugTfesCA[i]), log(ugKtsCA[i]), i);
          }
       // }//testing the hard emissions
      }
      if (ugKtsTau.size() > 0){
        //if( ugKtsTau[0] > kcut){
        //if( sqrt(ugKtsTau[0]*ugKtsTau[0] + (R/ugRgsTau[0])*(R/ugRgsTau[0]) ) < LPcut){
          for(int i = 0; i<ugRgsTau.size(); i++){
            if(ugRgsTau.size()<2) h3TauLundPlaneSigSingle[0]->Fill(log(R/ugRgsTau[i]), log(ugKtsTau[i]), i);
            if(ugRgsTau.size()<3) h3TauLundPlaneSigSingle[1]->Fill(log(R/ugRgsTau[i]), log(ugKtsTau[i]), i);
            if(ugRgsTau.size()<4) h3TauLundPlaneSigSingle[2]->Fill(log(R/ugRgsTau[i]), log(ugKtsTau[i]), i);
            if(ugRgsTau.size()<5) h3TauLundPlaneSigSingle[3]->Fill(log(R/ugRgsTau[i]), log(ugKtsTau[i]), i);
            if(ugRgsTau.size()<6) h3TauLundPlaneSigSingle[4]->Fill(log(R/ugRgsTau[i]), log(ugKtsTau[i]), i);
            if(ugRgsTau.size()<7) h3TauLundPlaneSigSingle[5]->Fill(log(R/ugRgsTau[i]), log(ugKtsTau[i]), i);
            if(ugRgsTau.size()<8) h3TauLundPlaneSigSingle[6]->Fill(log(R/ugRgsTau[i]), log(ugKtsTau[i]), i);
            h3TauLundPlaneSigSingleIncl->Fill(log(R/ugRgsTau[i]), log(ugKtsTau[i]), i);
          }
          for(int i = 0; i<ugTfesTau.size(); i++){
            if(ugTfesTau.size()<2) h3TauLundPlaneSigTfeSingle[0]->Fill(log(ugTfesTau[i]), log(ugKtsTau[i]), i);
            if(ugTfesTau.size()<3) h3TauLundPlaneSigTfeSingle[1]->Fill(log(ugTfesTau[i]), log(ugKtsTau[i]), i);
            if(ugTfesTau.size()<4) h3TauLundPlaneSigTfeSingle[2]->Fill(log(ugTfesTau[i]), log(ugKtsTau[i]), i);
            if(ugTfesTau.size()<5) h3TauLundPlaneSigTfeSingle[3]->Fill(log(ugTfesTau[i]), log(ugKtsTau[i]), i);
            if(ugTfesTau.size()<6) h3TauLundPlaneSigTfeSingle[4]->Fill(log(ugTfesTau[i]), log(ugKtsTau[i]), i);
            if(ugTfesTau.size()<7) h3TauLundPlaneSigTfeSingle[5]->Fill(log(ugTfesTau[i]), log(ugKtsTau[i]), i);
            if(ugTfesTau.size()<8) h3TauLundPlaneSigTfeSingle[6]->Fill(log(ugTfesTau[i]), log(ugKtsTau[i]), i);
            h3TauLundPlaneSigTfeSingleIncl->Fill(log(ugTfesTau[i]), log(ugKtsTau[i]), i);
          }
      //  }//testing the 'hard emissions'
      }

      // Here we are !
        std::vector<double> sigCAKt = sigJetRecurZcut_kt->at(i);
        std::vector<double> sigTauKt = sigJetRecurTauZcut_kt->at(i);
        //check for the case that we have a Zcut:
        //all jets splits in separet canvas
        //if( finsublled == false){//use nsubjettiness as discriminator to see how they change in terms of hardness
        if (sigCAKt.size() > 0){
          //if( sigCAKt[0] > kcut){
          //if( sqrt(sigCAKt[0]*sigCAKt[0] + (R/sigRgs[0])*(R/sigRgs[0]) ) < LPcut){
          for(int i = 0; i<sigRgs.size(); i++){
            if(sigRgs.size()<2) h3CALundPlaneSigSingleZ[0]->Fill(log(R/sigRgs[i]), log(sigCAKt[i]), i);
            if(sigRgs.size()<3) h3CALundPlaneSigSingleZ[1]->Fill(log(R/sigRgs[i]), log(sigCAKt[i]), i);
            if(sigRgs.size()<4) h3CALundPlaneSigSingleZ[2]->Fill(log(R/sigRgs[i]), log(sigCAKt[i]), i);
            if(sigRgs.size()<5) h3CALundPlaneSigSingleZ[3]->Fill(log(R/sigRgs[i]), log(sigCAKt[i]), i);
            if(sigRgs.size()<6) h3CALundPlaneSigSingleZ[4]->Fill(log(R/sigRgs[i]), log(sigCAKt[i]), i);
            if(sigRgs.size()<7) h3CALundPlaneSigSingleZ[5]->Fill(log(R/sigRgs[i]), log(sigCAKt[i]), i);
            if(sigRgs.size()<8) h3CALundPlaneSigSingleZ[6]->Fill(log(R/sigRgs[i]), log(sigCAKt[i]), i);
            h3CALundPlaneSigSingleInclZ->Fill(log(R/sigRgs[i]), log(sigCAKt[i]), i);//store number of emission for every jet separately
          }
          for(int i = 0; i<sigTfesZ.size(); i++){
            if(sigTfesZ.size()<2) h3CALundPlaneSigTfeSingleZ[0]->Fill(log(sigTfesZ[i]), log(sigCAKt[i]), i);
            if(sigTfesZ.size()<3) h3CALundPlaneSigTfeSingleZ[1]->Fill(log(sigTfesZ[i]), log(sigCAKt[i]), i);
            if(sigTfesZ.size()<4) h3CALundPlaneSigTfeSingleZ[2]->Fill(log(sigTfesZ[i]), log(sigCAKt[i]), i);
            if(sigTfesZ.size()<5) h3CALundPlaneSigTfeSingleZ[3]->Fill(log(sigTfesZ[i]), log(sigCAKt[i]), i);
            if(sigTfesZ.size()<6) h3CALundPlaneSigTfeSingleZ[4]->Fill(log(sigTfesZ[i]), log(sigCAKt[i]), i);
            if(sigTfesZ.size()<7) h3CALundPlaneSigTfeSingleZ[5]->Fill(log(sigTfesZ[i]), log(sigCAKt[i]), i);
            if(sigTfesZ.size()<8) h3CALundPlaneSigTfeSingleZ[6]->Fill(log(sigTfesZ[i]), log(sigCAKt[i]), i);
            h3CALundPlaneSigTfeSingleInclZ->Fill(log(sigTfesZ[i]), log(sigCAKt[i]), i);
          }
        //}//testing the hard emissions
      }
      if (sigTauKt.size() > 0){
        //if( sigTauKt[0] > kcut){
        //if( sqrt(sigTauKt[0]*sigTauKt[0] + (R/sigRgsTau[0])*(R/sigRgsTau[0]) ) < LPcut){
          for(int i = 0; i<sigRgsTau.size(); i++){
            if(sigRgsTau.size()<2) h3TauLundPlaneSigSingleZ[0]->Fill(log(R/sigRgsTau[i]), log(sigTauKt[i]), i);
            if(sigRgsTau.size()<3) h3TauLundPlaneSigSingleZ[1]->Fill(log(R/sigRgsTau[i]), log(sigTauKt[i]), i);
            if(sigRgsTau.size()<4) h3TauLundPlaneSigSingleZ[2]->Fill(log(R/sigRgsTau[i]), log(sigTauKt[i]), i);
            if(sigRgsTau.size()<5) h3TauLundPlaneSigSingleZ[3]->Fill(log(R/sigRgsTau[i]), log(sigTauKt[i]), i);
            if(sigRgsTau.size()<6) h3TauLundPlaneSigSingleZ[4]->Fill(log(R/sigRgsTau[i]), log(sigTauKt[i]), i);
            if(sigRgsTau.size()<7) h3TauLundPlaneSigSingleZ[5]->Fill(log(R/sigRgsTau[i]), log(sigTauKt[i]), i);
            if(sigRgsTau.size()<8) h3TauLundPlaneSigSingleZ[6]->Fill(log(R/sigRgsTau[i]), log(sigTauKt[i]), i);
            h3TauLundPlaneSigSingleInclZ->Fill(log(R/sigRgsTau[i]), log(sigTauKt[i]), i);
          }
          for(int i = 0; i<sigTfesZTau.size(); i++){
            if(sigTfesZTau.size()<2) h3TauLundPlaneSigTfeSingleZ[0]->Fill(log(sigTfesZTau[i]), log(sigTauKt[i]), i);
            if(sigTfesZTau.size()<3) h3TauLundPlaneSigTfeSingleZ[1]->Fill(log(sigTfesZTau[i]), log(sigTauKt[i]), i);
            if(sigTfesZTau.size()<4) h3TauLundPlaneSigTfeSingleZ[2]->Fill(log(sigTfesZTau[i]), log(sigTauKt[i]), i);
            if(sigTfesZTau.size()<5) h3TauLundPlaneSigTfeSingleZ[3]->Fill(log(sigTfesZTau[i]), log(sigTauKt[i]), i);
            if(sigTfesZTau.size()<6) h3TauLundPlaneSigTfeSingleZ[4]->Fill(log(sigTfesZTau[i]), log(sigTauKt[i]), i);
            if(sigTfesZTau.size()<7) h3TauLundPlaneSigTfeSingleZ[5]->Fill(log(sigTfesZTau[i]), log(sigTauKt[i]), i);
            if(sigTfesZTau.size()<8) h3TauLundPlaneSigTfeSingleZ[6]->Fill(log(sigTfesZTau[i]), log(sigTauKt[i]), i);
            h3TauLundPlaneSigTfeSingleInclZ->Fill(log(sigTfesZTau[i]), log(sigTauKt[i]), i);
          }
        //}//testing the 'hard emissions'
        }
        // LundPlane - for sigJet
        //std::cout<<"LP"<<std::endl;
        for(int i = 0; i<sigRgs.size(); i++){
          if(sigRgs.size()>0 && sigCAKt.size()>0) h2LogLundPlane->Fill(log(R/sigRgs[i]), log(sigCAKt[i]));//1st emissions only or all ?
          if(sigTfsZ.size()>0 && sigCAKt.size()>0) h2LogLundPlaneTf->Fill(log(sigTfsZ[i]), log(sigCAKt[i]));
          if(sigTfesZ.size()>0 && sigCAKt.size()>0) h2LogLundPlaneTfe->Fill(log(sigTfesZ[i]), log(sigCAKt[i]));
        }
        
        for(int i = 0; i<sigRgsTau.size(); i++){
          if(sigRgsTau.size()>0 && sigTauKt.size()>0) h2LogLundPlaneTau->Fill(log(R/sigRgsTau[i]), log(sigTauKt[i]));
          if(sigTfsZTau.size()>0 && sigTauKt.size()>0) h2LogLundPlaneTauTf->Fill(log(sigTfsZTau[i]), log(sigTauKt[i]));
          if(sigTfesZTau.size()>0 && sigTauKt.size()>0) h2LogLundPlaneTauTfe->Fill(log(sigTfesZTau[i]), log(sigTauKt[i]));
        }
        std::vector<double> droppedKtTau = sigJetRecurTauZcut_droppedKts->at(i);
        for(int i = 0; i<branchDroppedT.size(); i++){
          if(branchDroppedT.size()>0 && droppedKtTau.size()>0) h2LogLundPlaneTauBFS->Fill(log(R/branchDroppedT[i]), log(droppedKtTau[i]));
        }
        std::vector<double> droppedKt = sigJetRecurZcut_droppedKts->at(i);
        for(int i = 0; i<branchDropped.size(); i++){
          if(branchDropped.size()>0 && droppedKt.size()>0) h2LogLundPlaneBFS->Fill(log(R/branchDropped[i]), log(droppedKt[i])); 
        }
        // LundPlanes - for csFullJets
        //std::cout<<"LP csFull"<<std::endl;
        std::vector<double> csFullCAKt = csFullJetRecurZcut_kt->at(i);
        std::vector<double> csFullCARg = csFullJetRecurZcut_dr12->at(i);
        std::vector<double> csFullTauKt = csFullJetRecurTauZcut_kt->at(i);
        std::vector<double> csFullTauRg = csFullJetRecurTauZcut_dr12->at(i);
        //1st splits only
        if(sigRgs.size()>0 && sigCAKt.size()>0) h2Log1stLundPlane->Fill(log(R/sigRgs[0]), log(sigCAKt[0]));//1st emissions only or all ?
        if(sigTfsZ.size()>0 && sigCAKt.size()>0) h2Log1stLundPlaneTf->Fill(log(sigTfsZ[0]), log(sigCAKt[0]));
        if(sigTfesZ.size()>0 && sigCAKt.size()>0) h2Log1stLundPlaneTfe->Fill(log(sigTfesZ[0]), log(sigCAKt[0]));
        if(sigRgsTau.size()>0 && sigTauKt.size()>0) h2Log1stLundPlaneTau->Fill(log(R/sigRgsTau[0]), log(sigTauKt[0]));
        if(sigTfsZTau.size()>0 && sigTauKt.size()>0) h2Log1stLundPlaneTauTf->Fill(log(sigTfsZTau[0]), log(sigTauKt[0]));
        if(sigTfesZTau.size()>0 && sigTauKt.size()>0) h2Log1stLundPlaneTauTfe->Fill(log(sigTfesZTau[0]), log(sigTauKt[0]));
        if(branchDroppedT.size()>0 && droppedKtTau.size()>0) h2Log1stLundPlaneTauBFS->Fill(log(R/branchDroppedT[0]), log(droppedKtTau[0]));
        if(branchDropped.size()>0 && droppedKt.size()>0) h2Log1stLundPlaneBFS->Fill(log(R/branchDropped[0]), log(droppedKt[0])); 

        // compare 1->2 with 1->3 high E 
        std::vector<double> sigTfes2ZTau = sigJetRecurTauZcut_tfe2->at(i);
        std::vector<double> sigTfes2ZCA = sigJetRecurZcut_tfe2->at(i);
        std::vector<double> csTfesZCA = csFullJetRecurZcut_tfe->at(i);
        std::vector<double> csTfesZTau = csFullJetRecurTauZcut_tfe->at(i);
        std::vector<double> csTfes2ZTau = csFullJetRecurTauZcut_tfe2->at(i);
        std::vector<double> csTfes2ZCA = csFullJetRecurZcut_tfe2->at(i);

        if(sigTfesZ.size()>1 && sigTfes2ZCA.size()>0){
          h2SumCATfe2Tfe->Fill(sigTfesZ[0]+sigTfesZ[1], sigTfes2ZCA[0]);
        }
        if(sigTfesZTau.size()>1 && sigTfes2ZTau.size()>0){
          h2SumTauTfe2Tfe->Fill(sigTfesZTau[0]+sigTfesZTau[1], sigTfes2ZTau[0]);        
        }

        if(sigTfesZ.size()>0 && sigTfes2ZCA.size()>0){
          h2CATfe2Tfe->Fill(sigTfesZ[0], sigTfes2ZCA[0]);
          //tau1->3/tau1->2 in bins of tau1->2
          int bin = returnBinNumber(sigTfesZ[0],100,0,9);//for final tauBinning 0,9 in 0.1 -> 100Bins
          if( (bin >= 0) && (bin <= 99) ) h1CATfe2Tfe[bin]->Fill(sigTfes2ZCA[0]/sigTfesZ[0]);
          if( (bin >= 0) && (bin <= 99) ) hdiff1CATfe2Tfe[bin]->Fill((sigTfes2ZCA[0]-sigTfesZ[0])/sigTfesZ[0]);
          //tau1->2/tau1->3 in bins of tau1->3
          bin = returnBinNumber(sigTfes2ZCA[0],100,0,9);
          if( (bin >= 0) && (bin <= 99) ) h1CATfeTfe2[bin]->Fill(sigTfesZ[0]/sigTfes2ZCA[0]);
          if( (bin >= 0) && (bin <= 99) ) hdiff1CATfeTfe2[bin]->Fill((sigTfesZ[0]-sigTfes2ZCA[0])/sigTfes2ZCA[0]);
        }
        if(sigTfesZTau.size()>0 && sigTfes2ZTau.size()>0){
          h2TauTfe2Tfe->Fill(sigTfesZTau[0], sigTfes2ZTau[0]);
          //tau1->3/tau1->2 in bins of tau1->2
          int bin = returnBinNumber(sigTfesZTau[0],100,0,9);
          if( (bin >= 0) && (bin <= 99) ) h1TauTfe2Tfe[bin]->Fill(sigTfes2ZTau[0]/sigTfesZTau[0]);
          if( (bin >= 0) && (bin <= 99) ) hdiff1TauTfe2Tfe[bin]->Fill((sigTfes2ZTau[0]-sigTfesZTau[0])/sigTfesZTau[0]);
          //tau1->2/tau1->3 in bins of tau1->3
          bin = returnBinNumber(sigTfes2ZTau[0],100,0,9);
          if( (bin >= 0) && (bin <= 99) ) h1TauTfeTfe2[bin]->Fill(sigTfesZTau[0]/sigTfes2ZTau[0]);
          if( (bin >= 0) && (bin <= 99) ) hdiff1TauTfeTfe2[bin]->Fill((sigTfesZTau[0]-sigTfes2ZTau[0])/sigTfes2ZTau[0]);
        } 

        //ungroomed timing
        std::vector<double> sigTfesTauUg = sigJetRecurTau_tfe->at(i);
        std::vector<double> sigTfesCAUg = sigJetRecur_tfe->at(i);
        std::vector<double> sigTfes2TauUg = sigJetRecurTau_tfe2->at(i);
        std::vector<double> sigTfes2CAUg = sigJetRecur_tfe2->at(i);
        std::vector<double> csTfesCAUg = csFullJetRecur_tfe->at(i);
        std::vector<double> csTfesTauUg = csFullJetRecurTau_tfe->at(i);
        std::vector<double> csTfes2TauUg = csFullJetRecurTau_tfe2->at(i);
        std::vector<double> csTfes2CAUg = csFullJetRecur_tfe2->at(i);
        if(sigTfes2ZCA.size()>0){
          h1Ntf2ZCA->Fill(sigTfes2ZCA[0]);
        }
        if(sigTfes2ZTau.size()>0){
          h1Ntf2ZTau->Fill(sigTfes2ZTau[0]);
        }
        if(sigTfes2CAUg.size()>0){
          h1Ntf2CA->Fill(sigTfes2CAUg[0]);
        }
        if(sigTfes2TauUg.size()>0){
          h1Ntf2Tau->Fill(sigTfes2TauUg[0]);
        }
        if(sigTfesZ.size()>0){
          h1NtfZCA->Fill(sigTfesZ[0]);
        }
        if(sigTfesZTau.size()>0){
          h1NtfZTau->Fill(sigTfesZTau[0]);
        }
        if(sigTfesCAUg.size()>0){
          h1NtfCA->Fill(sigTfesCAUg[0]);
        }
        if(sigTfesTauUg.size()>0){
          h1NtfTau->Fill(sigTfesTauUg[0]);
        }

        if(sigTfesCAUg.size()>0 && sigTfes2CAUg.size()>0){
          h2CATfe2TfeSig->Fill(sigTfesCAUg[0], sigTfes2CAUg[0]);
          //tau1->3/tau1->2 in bins of tau1->2
          int bin = returnBinNumber(sigTfesCAUg[0],100,0,9);//for final tauBinning 0,9 in 0.1 -> 100Bins
          if( (bin >= 0) && (bin <= 99) ) h1CATfe2TfeSig[bin]->Fill(sigTfes2CAUg[0]/sigTfesCAUg[0]);
          if( (bin >= 0) && (bin <= 99) ) hdiff1CATfe2TfeSig[bin]->Fill((sigTfes2CAUg[0]-sigTfesCAUg[0])/sigTfesCAUg[0]);
          //tau1->2/tau1->3 in bins of tau1->3
          bin = returnBinNumber(sigTfes2CAUg[0],100,0,9);
          if( (bin >= 0) && (bin <= 99) ) h1CATfeTfe2Sig[bin]->Fill(sigTfesCAUg[0]/sigTfes2CAUg[0]);
          if( (bin >= 0) && (bin <= 99) ) hdiff1CATfeTfe2Sig[bin]->Fill((sigTfesCAUg[0]-sigTfes2CAUg[0])/sigTfes2CAUg[0]);
        }
        if(sigTfesTauUg.size()>0 && sigTfes2TauUg.size()>0){
          h2TauTfe2TfeSig->Fill(sigTfesTauUg[0], sigTfes2TauUg[0]);
          //tau1->3/tau1->2 in bins of tau1->2
          int bin = returnBinNumber(sigTfesTauUg[0],100,0,9);
          if( (bin >= 0) && (bin <= 99) ) h1TauTfe2TfeSig[bin]->Fill(sigTfes2TauUg[0]/sigTfesTauUg[0]);
          if( (bin >= 0) && (bin <= 99) ) hdiff1TauTfe2TfeSig[bin]->Fill((sigTfes2TauUg[0]-sigTfesTauUg[0])/sigTfesTauUg[0]);
          //tau1->2/tau1->3 in bins of tau1->3
          bin = returnBinNumber(sigTfes2TauUg[0],100,0,9);
          if( (bin >= 0) && (bin <= 99) ) h1TauTfeTfe2Sig[bin]->Fill(sigTfesTauUg[0]/sigTfes2TauUg[0]);
          if( (bin >= 0) && (bin <= 99) ) hdiff1TauTfeTfe2Sig[bin]->Fill((sigTfesTauUg[0]-sigTfes2TauUg[0])/sigTfes2TauUg[0]);
        } 
        //jetMatching original position
          if(csTauPtVec.size()>0){
            h1SigTauCs->Fill(sigPt-csTauPtVec[0]);
            h2SigTauCs->Fill(sigPt,csTauPtVec[0]);
            if(sigPt-csTauPtVec[0] > 70){
              std::cout<<"Sig pt: "<<sigPt<<" eta: "<<sigJetEta->at(i)<<" phi: "<<sigJetPhi->at(i)<<" m: "<<sigJetM->at(i)<<std::endl;
              std::cout<<"Cs pt: "<<csTauPtVec[0]<<std::endl;

            }
          }
          if(csCAPtVec.size()>0){
            h1SigCaCs->Fill(sigPt-csCAPtVec[i]);
            h2SigCaCs->Fill(sigPt, csCAPtVec[i]);
          }
          if(csCAPtVec.size()>0 && csTauPtVec.size()>0){
            h1CaTauCs->Fill(csCAPtVec[i]-csTauPtVec[i]);
            h2CaTauCs->Fill(csCAPtVec[i],csTauPtVec[i]);
          }


          
          for(int i = 0; i<csFullCARg.size(); i++){
            if(csFullCARg.size()>0 && csFullCAKt.size()>0) h2LogLundPlaneCs->Fill(log(R/csFullCARg[i]), log(csFullCAKt[i]));
            if(csTfs.size()>0 && csFullCAKt.size()>0) h2LogLundPlaneTfCs->Fill(log(csTfs[i]), log(csFullCAKt[i]));
            if(csTfes.size()>0 && csFullCAKt.size()>0) h2LogLundPlaneTfeCs->Fill(log(csTfes[i]), log(csFullCAKt[i]));
          }
          for(int i = 0; i<csFullTauRg.size(); i++){
            if(csFullTauRg.size()>0 && csFullTauKt.size()>0) h2LogLundPlaneTauCs->Fill(log(R/csFullTauRg[i]), log(csFullTauKt[i]));
            if(csTauTfs.size()>0 && csFullTauKt.size()>0) h2LogLundPlaneTauTfCs->Fill(log(csTauTfs[i]), log(csFullTauKt[i]));
            if(csTauTfes.size()>0 && csFullTauKt.size()>0) h2LogLundPlaneTauTfeCs->Fill(log(csTauTfes[i]), log(csFullTauKt[i]));
          }

          if(csFullCARg.size()>0 && csFullCAKt.size()>0) h2Log1stLundPlaneCs->Fill(log(R/csFullCARg[0]), log(csFullCAKt[0]));
          if(csTfs.size()>0 && csFullCAKt.size()>0) h2Log1stLundPlaneTfCs->Fill(log(csTfs[0]), log(csFullCAKt[0]));
          if(csTfes.size()>0 && csFullCAKt.size()>0) h2Log1stLundPlaneTfeCs->Fill(log(csTfes[0]), log(csFullCAKt[0]));
          if(csFullTauRg.size()>0 && csFullTauKt.size()>0) h2Log1stLundPlaneTauCs->Fill(log(R/csFullTauRg[0]), log(csFullTauKt[0]));
          if(csTauTfs.size()>0 && csFullTauKt.size()>0) h2Log1stLundPlaneTauTfCs->Fill(log(csTauTfs[0]), log(csFullTauKt[0]));
          if(csTauTfes.size()>0 && csFullTauKt.size()>0) h2Log1stLundPlaneTauTfeCs->Fill(log(csTauTfes[0]), log(csFullTauKt[0]));

          if(csTfesZCA.size()>1 && csTfes2ZCA.size()>0){
            h2SumCATfe2TfeCs->Fill(csTfesZCA[0]+csTfesZCA[1], csTfes2ZCA[0]);
          }
          if(csTfesZTau.size()>1 && csTfes2ZTau.size()>0){
            h2SumTauTfe2TfeCs->Fill(csTfesZTau[0]+csTfesZTau[1], csTfes2ZTau[0]);
          }  
          if(csTfesZCA.size()>0 && csTfes2ZCA.size()>0){
            h2CATfe2TfeCs->Fill(csTfesZCA[0], csTfes2ZCA[0]);
            //tau1->3/tau1->2 in bins of tau1->2
            int bin = returnBinNumber(csTfesZCA[0],100,0,9);
            if( (bin >= 0) && (bin <= 99) ) h1CATfe2TfeCs[bin]->Fill(csTfes2ZCA[0]/csTfesZCA[0]);
            if( (bin >= 0) && (bin <= 99) ) hdiff1CATfe2TfeCs[bin]->Fill((csTfes2ZCA[0]-csTfesZCA[0])/csTfesZCA[0]);
            //tau1->2/tau1->3 in bins of tau1->3
            
            bin = returnBinNumber(csTfes2ZCA[0],100,0,9);
            if( (bin >= 0) && (bin <= 99) ) h1CATfeTfe2Cs[bin]->Fill(csTfesZCA[0]/csTfes2ZCA[0]);
            if( (bin >= 0) && (bin <= 99) ) hdiff1CATfeTfe2Cs[bin]->Fill((csTfesZCA[0]-csTfes2ZCA[0])/csTfes2ZCA[0]);
          } 
          if(csTfesZTau.size()>0 && csTfes2ZTau.size()>0){
            h2TauTfe2TfeCs->Fill(csTfesZTau[0], csTfes2ZTau[0]);
            //tau1->3/tau1->2 in bins of tau1->2
            int bin = returnBinNumber(csTfesZTau[0],100,0,9);
            if( (bin >= 0) && (bin <= 99) ) h1TauTfe2TfeCs[bin]->Fill(csTfes2ZTau[0]/csTfesZTau[0]);
            if( (bin >= 0) && (bin <= 99) ) hdiff1TauTfe2TfeCs[bin]->Fill((csTfes2ZTau[0]-csTfesZTau[0])/csTfesZTau[0]);
            //tau1->2/tau1->3 in bins of tau1->3
            bin = returnBinNumber(csTfes2ZTau[0],100,0,9);
            if( (bin >= 0) && (bin <= 99) ) h1TauTfeTfe2Cs[bin]->Fill(csTfesZTau[0]/csTfes2ZTau[0]);
            if( (bin >= 0) && (bin <= 99) ) hdiff1TauTfeTfe2Cs[bin]->Fill((csTfesZTau[0]-csTfes2ZTau[0])/csTfes2ZTau[0]);
          }
          //timing
          //and for csFullJets
          if(csTfes2ZCA.size()>0){
            h1Ncstf2ZCA->Fill(csTfes2ZCA[0]);
          }
          if(csTfes2ZTau.size()>0){
            h1Ncstf2ZTau->Fill(csTfes2ZTau[0]);
          }
          if(csTfes2CAUg.size()>0){
            h1Ncstf2CA->Fill(csTfes2CAUg[0]);
          }
          if(csTfes2TauUg.size()>0){
            h1Ncstf2Tau->Fill(csTfes2TauUg[0]);
          }
          if(csTfesZCA.size()>0){
            h1NcstfZCA->Fill(csTfesZCA[0]);
          }
          if(csTfesZTau.size()>0){
            h1NcstfZTau->Fill(csTfesZTau[0]);
          }
          if(csTfesCAUg.size()>0){
            h1NcstfCA->Fill(csTfesCAUg[0]);
          }
          if(csTfesTauUg.size()>0){
            h1NcstfTau->Fill(csTfesTauUg[0]);
          }
          //Response csFull/sigJet
          //no grooming
          //1->2
          if(csTfesCAUg.size()>0 && sigTfesCAUg.size()>0){
            h2TfeRespCA->Fill(sigTfesCAUg[0], csTfesCAUg[0]/sigTfesCAUg[0]);
          }
          if(csTfesTauUg.size()>0 && sigTfesTauUg.size()>0){
            h2TfeRespTau->Fill(sigTfesTauUg[0], csTfesTauUg[0]/sigTfesTauUg[0]);
          }
          //1->3
          if(csTfes2CAUg.size()>0 && sigTfes2CAUg.size()>0){
            h2Tfe2RespCA->Fill(sigTfes2CAUg[0], csTfes2CAUg[0]/sigTfes2CAUg[0]);
          }
          if(csTfes2TauUg.size()>0 && sigTfes2TauUg.size()>0){
            h2Tfe2RespTau->Fill(sigTfes2TauUg[0], csTfes2TauUg[0]/sigTfes2TauUg[0]);
          }
          //1->2
          if(csTfesZCA.size()>0 && sigTfesZ.size()>0){
            h2TfeRespCAZcut->Fill(sigTfesZ[0], csTfesZCA[0]/sigTfesZ[0]);
          }
          if(csTfesZTau.size()>0 && sigTfesZTau.size()>0){
            h2TfeRespTauZcut->Fill(sigTfesZTau[0], csTfesZTau[0]/sigTfesZTau[0]);
          }
          //1->3
          if(csTfes2ZCA.size()>0 && sigTfes2ZCA.size()>0){
            h2Tfe2RespCAZcut->Fill(sigTfes2ZCA[0], csTfes2ZCA[0]/sigTfes2ZCA[0]);
          }
          if(csTfes2ZTau.size()>0 && sigTfes2ZTau.size()>0){
            h2Tfe2RespTauZcut->Fill(sigTfes2ZTau[0], csTfes2ZTau[0]/sigTfes2ZTau[0]);
          }
          //2D tfe and ratio/difference
          if(csTfesCAUg.size()>0 && csTfes2CAUg.size()>0){
            h2CATfe2TfeCsUg->Fill(csTfesCAUg[0], csTfes2CAUg[0]);
            //tau1->3/tau1->2 in bins of tau1->2
            int bin = returnBinNumber(csTfesCAUg[0],100,0,9);
            if( (bin >= 0) && (bin <= 99) ) h1CATfe2TfeCsUg[bin]->Fill(csTfes2CAUg[0]/csTfesCAUg[0]);
            if( (bin >= 0) && (bin <= 99) ) hdiff1CATfe2TfeCsUg[bin]->Fill((csTfes2CAUg[0]-csTfesCAUg[0])/csTfesCAUg[0]);
            //tau1->2/tau1->3 in bins of tau1->3
            bin = returnBinNumber(csTfes2CAUg[0],100,0,9);
            if( (bin >= 0) && (bin <= 99) ) h1CATfeTfe2CsUg[bin]->Fill(csTfesCAUg[0]/csTfes2CAUg[0]);
            if( (bin >= 0) && (bin <= 99) ) hdiff1CATfeTfe2CsUg[bin]->Fill((csTfesCAUg[0]-csTfes2CAUg[0])/csTfes2CAUg[0]);
          } 
          if(csTfesTauUg.size()>0 && csTfes2TauUg.size()>0){
            h2TauTfe2TfeCsUg->Fill(csTfesTauUg[0], csTfes2TauUg[0]);
            //tau1->3/tau1->2 in bins of tau1->2
            int bin = returnBinNumber(csTfesTauUg[0],100,0,9);
            if( (bin >= 0) && (bin <= 99) ) h1TauTfe2TfeCsUg[bin]->Fill(csTfes2TauUg[0]/csTfesTauUg[0]);
            if( (bin >= 0) && (bin <= 99) ) hdiff1TauTfe2TfeCsUg[bin]->Fill((csTfes2TauUg[0]-csTfesTauUg[0])/csTfesTauUg[0]);
            //tau1->2/tau1->3 in bins of tau1->3
            bin = returnBinNumber(csTfes2TauUg[0],100,0,9);
            if( (bin >= 0) && (bin <= 99) ) h1TauTfeTfe2CsUg[bin]->Fill(csTfesTauUg[0]/csTfes2TauUg[0]);
            if( (bin >= 0) && (bin <= 99) ) hdiff1TauTfeTfe2CsUg[bin]->Fill((csTfesTauUg[0]-csTfes2TauUg[0])/csTfes2TauUg[0]);
          }
        
      }//end of sigJetPtLoop
    for(int i = 0; i<partonWd1Wod2pt->size(); i++){
      h1PartonDaughter1Pt->Fill(partonWd1Wod2pt->at(i));
    }
    for(int i = 0; i<partonWod1Wd2pt->size(); i++){
      h1PartonDaughter2Pt->Fill(partonWod1Wd2pt->at(i));
    }
    for(int i = 0; i<partonWd1d2pt->size(); i++){
      h1Parton2DaughtersPt->Fill(partonWd1d2pt->at(i));
    }
    for(int i = 0; i<partonWod1d2pt->size(); i++){
      h1PartonNoDaughterPt->Fill(partonWod1d2pt->at(i));
    }
   */ 
   
   }
   TFile *fout;
   if(reduced == true ){
    fout = new TFile(fin+generator+"_"+duration+".root","RECREATE");
    //fout = new TFile("Johanna/Jewel/JetToyHITimeClusHistosJewel_simple_"+generator+"_"+duration+".root","RECREATE");
   }else{
    std::cout<<"create file"<<std::endl;
    fout = new TFile(fin+"_Histos.root","RECREATE");
    //fout = new TFile("Johanna/Jewel/JetToyHITimeClusHistosJewel_simple.root","RECREATE");
  }

  //add new histos to write !
   std::cout << "write content" << std::endl;
   //event info
   hEntries->Write();
   //parton info
   h2PartonSigTf->Write();
   h2PartonSigTfTau->Write();
   h2PartonSigTfZcut->Write();
   h2PartonSigTfZcutTau->Write();

   h2PartonSigTfeZcutTau->Write();
   h2PartonSigTfeZcut->Write();

   h2SigTfTfe->Write();
   h2SigTfTfeTau->Write();
   h2SigTfTfeZcut->Write();
   h2SigTfTfeZcutTau->Write();

   h2CsPartonTfe->Write();
   h2CsPartonTfeTau->Write();
   h2CsSigTfeZcut->Write();
   h2CsSigTfeZcutTau->Write();
   
   h1LogTfSig->Write();
   h1LogTfSigTau->Write();
   h1LogTfSigZ->Write();
   h1LogTfSigZTau->Write();
   h1LogTfSigToPartonZ->Write();
   h1LogTfSigToPartonZTau->Write();

   //sig jets theta zg
   h2ThetaZgCAZcut->Write();
   h2ThetaZgTauZcut->Write();
   h2ThetaZgCA->Write();
   h2ThetaZgTau->Write();
   h2CATfe2Tfe->Write();
   h2TauTfe2Tfe->Write();
   h2CATfe2TfeCs->Write();
   h2TauTfe2TfeCs->Write();
   h2CATfe2TfeSig->Write();
   h2TauTfe2TfeSig->Write();
   h2CATfe2TfeCsUg->Write();
   h2TauTfe2TfeCsUg->Write();

   h2SumCATfe2Tfe->Write();
   h2SumTauTfe2Tfe->Write();
   h2SumCATfe2TfeCs->Write();
   h2SumTauTfe2TfeCs->Write();
   //Number of formation times
   h1Ntf2ZCA->Write();//1->3
   h1Ntf2ZTau->Write();
   h1Ntf2CA->Write();
   h1Ntf2Tau->Write();
   h1NtfZCA->Write();//1->2
   h1NtfZTau->Write();
   h1NtfCA->Write();
   h1NtfTau->Write();
   //cs
   h1Ncstf2ZCA->Write();//1->3
   h1Ncstf2ZTau->Write();
   h1Ncstf2CA->Write();
   h1Ncstf2Tau->Write();
   h1NcstfZCA->Write();//1->2
   h1NcstfZTau->Write();
   h1NcstfCA->Write();
   h1NcstfTau->Write();

   //Lund Planes Novel Tools Fig. 3
   h2CALundPlaneNT->Write();
   h2TauLundPlaneNT->Write();
   // parton level (0-th split) LP's
   h2PartonLogTfKt->Write();
   h2PartonLogDrKt->Write();
   h2PartonLundPlaneNT->Write();
   h2PartonLogTfKtZcut->Write();
   h2PartonLogDrKtZcut->Write();
   // single jet LP
   h3CALundPlaneSigTfeSingleIncl->Write();
   h3TauLundPlaneSigTfeSingleIncl->Write();   
   h3CALundPlaneSigTfeSingleInclZ->Write();
   h3TauLundPlaneSigTfeSingleInclZ->Write();
   h2CALundPlaneSigTfeSingle->Write();
   h2TauLundPlaneSigTfeSingle->Write();
   //log (kt) = log(pt*dr/R) vs log(R/dr)
   h3CALundPlaneSigSingleIncl->Write();
   h3TauLundPlaneSigSingleIncl->Write();
   h3CALundPlaneSigSingleInclZ->Write();
   h3TauLundPlaneSigSingleInclZ->Write();
   h2CALundPlaneSigSingle->Write();
   h2TauLundPlaneSigSingle->Write();
   // ungroomed sigJets
   h2CALundPlaneSig->Write();
   h2TauLundPlaneSig->Write();
   h2CALundPlaneSigTfe->Write();
   h2TauLundPlaneSigTfe->Write();
   // ungroomed Cs
   h2CALundPlaneCs->Write();
   h2TauLundPlaneCs->Write();
   h2CALundPlaneCsTfe->Write();
   h2TauLundPlaneCsTfe->Write();
   //Groomed sigJets
   h2LogLundPlane->Write();//LD C/A
   h2LogLundPlaneTau->Write();//LP tau recl
   h2LogLundPlaneTf->Write();
   h2LogLundPlaneTauTf->Write();
   h2LogLundPlaneTfe->Write();
   h2LogLundPlaneTauTfe->Write();
   h2LogLundPlaneBFS->Write();//LD C/A
   h2LogLundPlaneTauBFS->Write();//LP tau recl
   //Groomed embedded 
   h2LogLundPlaneCs->Write();//LD C/A - constituent subtracted
   h2LogLundPlaneTfCs->Write();
   h2LogLundPlaneTfeCs->Write();
   h2LogLundPlaneTauCs->Write();//LP tau recl - constituent subtracted
   h2LogLundPlaneTauTfeCs->Write();
   h2LogLundPlaneTauTfCs->Write();
      
   h2Log1stLundPlane->Write();//LD C/A
   h2Log1stLundPlaneTau->Write();//LP tau recl
   h2Log1stLundPlaneTf->Write();
   h2Log1stLundPlaneTauTf->Write();
   h2Log1stLundPlaneTfe->Write();
   h2Log1stLundPlaneTauTfe->Write();
   h2Log1stLundPlaneBFS->Write();//LD C/A
   h2Log1stLundPlaneTauBFS->Write();//LP tau recl
   h2Log1stLundPlaneCs->Write();//LD C/A - constituent subtracted
   h2Log1stLundPlaneTfCs->Write();
   h2Log1stLundPlaneTfeCs->Write();
   h2Log1stLundPlaneTauCs->Write();//LP tau recl - constituent subtracted
   h2Log1stLundPlaneTauTfeCs->Write();
   h2Log1stLundPlaneTauTfCs->Write();

   h1SigTauCs->Write();
   h1SigCaCs->Write();
   h1CaTauCs->Write();
   h1Ca->Write();
   h1Tau->Write();
   
   h2SigTauCs->Write();
   h2SigCaCs->Write();
   h2CaTauCs->Write();
   h2Ca->Write();
   h2Tau->Write();

   h2SigTfeTauTfe->Write();
   h2SigTfeTauTfeZcut->Write();
   h1SigNCaTf->Write();
   h1SigNCaTfZcut->Write();
   h1SigNTauTf->Write();
   h1SigNTauTfZcut->Write();

   h2DroppedDrTauZcut->Write();
   h2DroppedPtTauZcut->Write();
   h2DroppedKtTauZcut->Write();
   h2DroppedLPTauZcut->Write();
   h2DroppedDrZcut->Write();
   h2DroppedPtZcut->Write();
   h2DroppedKtZcut->Write();
   h2DroppedLPZcut->Write();



   h2SigRgTauRgZcut->Write();
   h2SigZgTauZgZcut->Write();
  
   h1umsigTauPt->Write();
   h1umsigTauEta->Write();
   h1umsigCAEta->Write();
   h1umsigCAPt->Write();
   h1matchedsigTauPt->Write();
   h1matchedsigTauEta->Write();
   h1matchedsigCAEta->Write();
   h1matchedsigCAPt->Write();

  
   h1JPtW->Write();
   h1JPt->Write();
   h1JsigTauPtW->Write();
   h1JsigTauPt->Write();
   h1JsigCAPtW->Write();
   h1JsigCAPt->Write();

   //Parton pt's givien certain daughter constallations
   h1PartonNoDaughterPt->Write();
   h1Parton2DaughtersPt->Write();
   h1PartonDaughter1Pt->Write();
   h1PartonDaughter2Pt->Write();

   //Pythia information
   h3tfs->Write();
   h3tfes->Write();
   h3drs->Write();
   h3tftsfesdrs->Write();
   h3jetsplittfs->Write();
   h3jetsplittfes->Write();
   h3jetsplitdrs->Write();
   h1DroppedBr->Write();
   h1DroppedBrTau->Write();
   h2Nsubj21->Write();
   h2Nsubj32->Write();
   h2TauNsubj21->Write();

   h2CaTauNsubj21->Write();
   h2CaTauNsubj32->Write();

  for (int NBin = 0; NBin < NBins; NBin++){
    //tau1->3/tau1->2 in bins of tau1->2
    h1CATfe2Tfe[NBin]->Write();
    h1TauTfe2Tfe[NBin]->Write();
    h1CATfe2TfeCs[NBin]->Write();
    h1TauTfe2TfeCs[NBin]->Write();
    h1CATfe2TfeSig[NBin]->Write();
    h1TauTfe2TfeSig[NBin]->Write();
    h1CATfe2TfeCsUg[NBin]->Write();
    h1TauTfe2TfeCsUg[NBin]->Write();
    //(tau1->3 - tau1->2) /tau1->2 in bins of tau1->2
    hdiff1CATfe2Tfe[NBin]->Write();
    hdiff1TauTfe2Tfe[NBin]->Write();
    hdiff1CATfe2TfeCs[NBin]->Write();
    hdiff1TauTfe2TfeCs[NBin]->Write();
    hdiff1CATfe2TfeSig[NBin]->Write();
    hdiff1TauTfe2TfeSig[NBin]->Write();
    hdiff1CATfe2TfeCsUg[NBin]->Write();
    hdiff1TauTfe2TfeCsUg[NBin]->Write();
    //tau1->2/tau1->3 in bins of tau1->3
    h1CATfeTfe2[NBin]->Write();
    h1TauTfeTfe2[NBin]->Write();
    h1CATfeTfe2Cs[NBin]->Write();
    h1TauTfeTfe2Cs[NBin]->Write();
    h1CATfeTfe2Sig[NBin]->Write();
    h1TauTfeTfe2Sig[NBin]->Write();
    h1CATfeTfe2CsUg[NBin]->Write();
    h1TauTfeTfe2CsUg[NBin]->Write();
    //(tau1->2 - tau1->3) /tau1->3 in bins of tau1->3
    hdiff1CATfeTfe2[NBin]->Write();
    hdiff1TauTfeTfe2[NBin]->Write();
    hdiff1CATfeTfe2Cs[NBin]->Write();
    hdiff1TauTfeTfe2Cs[NBin]->Write();
    hdiff1CATfeTfe2Sig[NBin]->Write();
    hdiff1TauTfeTfe2Sig[NBin]->Write();
    hdiff1CATfeTfe2CsUg[NBin]->Write();
    hdiff1TauTfeTfe2CsUg[NBin]->Write();
   }

   //Jet substructure 
   for (int nsplit = 0; nsplit < Nsplit; nsplit++){
    //fill 0 with 1st split, 1 with 2nd and 3 with all.
    h1TfSigTau[nsplit]->Write();
    h1TfSigCA[nsplit]->Write();
    h1ZgSigCA[nsplit]->Write();
    h1ZgSigTau[nsplit]->Write();
    h1RgCASig[nsplit]->Write();
    h1RgTauSig[nsplit]->Write();
    h1EradSigCA[nsplit]->Write();
    h1EradSigTau[nsplit]->Write();
    h1LogDr12SigCA[nsplit]->Write();
    h1LogDr12SigTau[nsplit]->Write();
    h1LogZthetaSigCA[nsplit]->Write();
    h1LogZthetaSigTau[nsplit]->Write();
    
   }

   for (int nsplit = 0; nsplit < Nemission; nsplit++){//LundPlane evolution
    h3CALundPlaneSigTfeSingle[nsplit]->Write();
    h3TauLundPlaneSigTfeSingle[nsplit]->Write();
    h3CALundPlaneSigSingle[nsplit]->Write();
    h3TauLundPlaneSigSingle[nsplit]->Write();

    h3CALundPlaneSigTfeSingleZ[nsplit]->Write();
    h3TauLundPlaneSigTfeSingleZ[nsplit]->Write();
    h3CALundPlaneSigSingleZ[nsplit]->Write();
    h3TauLundPlaneSigSingleZ[nsplit]->Write();
   }

   h2TfeRespCA->Write();
   h2TfeRespTau->Write();
   h2Tfe2RespCA->Write();
   h2Tfe2RespTau->Write();
   //with zcut
   h2TfeRespCAZcut->Write();
   h2TfeRespTauZcut->Write();
   h2Tfe2RespCAZcut->Write();
   h2Tfe2RespTauZcut->Write();
   
   //fSplitList->Write();
   fout->Write();
   fout->Close();

}
