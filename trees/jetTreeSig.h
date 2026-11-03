//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar  2 14:47:15 2021 by ROOT version 6.12/06
// from TTree jetTreeSig/JetToyHI tree
// found on file: JetToyHIResultTimeClusBkg_Merged.root
//////////////////////////////////////////////////////////

#ifndef jetTreeSig_h
#define jetTreeSig_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class jetTreeSig {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   // Parton info only in pythia.
   vector<double>  *eventWeight;
   vector<double>  *partonsPt;
   vector<double>  *partonsEta;
   vector<double>  *partonsPhi;
   vector<double>  *partonsM;
   vector<int>     *partonsPDG;
   vector<double>  *partonsFirstSplitPt;
   vector<double>  *partonsFirstSplitEta;
   vector<double>  *partonsFirstSplitPhi;
   vector<double>  *partonsFirstSplitM;
   vector<int>     *partonsFirstSplitPDG;

   vector<double>  *drsplit;
   vector<double>  *ktsplit;
   vector<double>  *tfsplit;
   vector<double>  *tfesplit;
   vector<double>  *partonWd1Wod2pt;
   vector<double>  *partonWod1Wd2pt;
   vector<double>  *partonWd1d2pt;
   vector<double>  *partonWod1d2pt;
   vector<double>  *zgsplit;
   vector<double>  *sigJetPt;
   vector<double>  *sigJetEta;
   vector<double>  *sigJetPhi;
   vector<double>  *sigJetM;
   vector<double>  *sigJetArea;

   vector<vector<double>> *drsplits;
   vector<vector<double>> *tfsplits;
   vector<vector<double>> *tfesplits;

   vector<vector<double>> *jetCollectionSigDYEPt;
   vector<double>  *jetCollectionSigDYEEta;

   vector<double>  *sigJetRecurZcut_nSD;
   vector<double>  *sigJetRecurZcut_zSD;
   vector<double>  *sigJetRecur_nSD;
   vector<double>  *sigJetRecur_partonMatchDr;
   vector<double>  *sigJetRecur_zSD;
   //add vectors for dynGrooming
   vector<vector<double> > *sigJetDYE_z;         //zg of groomed jets
   vector<vector<double> > *sigJetDYE_dr12;       //distance between the two subjet
   vector<vector<double> > *sigJetDYE_kt;         //kt of groomed jets
   vector<vector<double> > *sigJetDYE_tau21;      //n-subjettiness ratio 21
   vector<vector<double> > *sigJetDYE_tau32;      //n-subjettiness ratio 32
   vector<vector<double> > *sigJetDYE_tf;     	// formation time tf = 2 / (z*(1-z)*pt*theta^2) (soft limit)
   vector<vector<double> > *sigJetDYE_tfe;     	// formation time tf = 1 / 2*(z*(1-z)*E*(1-cos(theta))) (hard limit)
   vector<vector<double> > *sigJetkappaDYE;      //kappa of groomed jets
   vector<vector<double> > *sigJetDYE_drBranches;      //number of dropped branches
   vector<vector<double> > *sigJetTauDYE_z;
   vector<vector<double> > *sigJetTauDYE_dr12;
   vector<vector<double> > *sigJetTauDYE_kt;
   vector<vector<double> > *sigJetTauDYE_tau21;
   vector<vector<double> > *sigJetTauDYE_tau32;
   vector<vector<double> > *sigJetTauDYE_tf;
   vector<vector<double> > *sigJetTauDYE_tfe;
   vector<vector<double> > *sigJetTaukappaDYE;
   vector<vector<double> > *sigJetTauDYE_drBranches;

   vector<vector<double> > *sigJetRecurZcut_dr12;
   vector<vector<double> > *sigJetRecurZcut_erad;
   vector<vector<double> > *sigJetRecurZcut_jetpt;
   vector<vector<double> > *sigJetRecurZcut_jeteta;
   vector<vector<double> > *sigJetRecurZcut_logdr12;
   vector<vector<double> > *sigJetRecurZcut_logztheta;
   vector<vector<double> > *sigJetRecurZcut_tf;
   vector<vector<double> > *sigJetRecurZcut_tfe;
   vector<vector<double> > *sigJetRecurZcut_tfe2;
   vector<vector<double> > *sigJetRecurZcut_z;
   vector<vector<double> > *sigJetRecurZcut_kt;
   vector<vector<double> > *sigJetRecurZcut_droppedBFS;
   vector<vector<double> > *sigJetRecurZcut_droppedTfBFS;
   vector<vector<double> > *sigJetRecurZcut_droppedTfeBFS;
   vector<vector<double> > *sigJetRecurZcut_droppedPts;
   vector<vector<double> > *sigJetRecurZcut_droppedKts;
   vector<vector<double> > *sigJetRecurZcut_droppedLog1drBFS;
   vector<double> *sigJetRecurZcut_tau21;
   vector<double> *sigJetRecurZcut_tau32;


   vector<vector<double> > *sigJetRecur_dr12;
   vector<vector<double> > *sigJetRecur_erad;
   vector<vector<double> > *sigJetRecur_jetpt;
   vector<vector<double> > *sigJetRecur_jeteta;
   vector<vector<double> > *sigJetRecur_logdr12;
   vector<vector<double> > *sigJetRecur_logztheta;
   vector<vector<double> > *sigJetRecur_tf;
   vector<vector<double> > *sigJetRecur_tfe;
   vector<vector<double> > *sigJetRecur_tfe2;
   vector<vector<double> > *sigJetRecur_kt;
   vector<vector<double> > *sigJetRecur_z;
   vector<vector<double> > *sigJetRecur_droppedBFS;
   vector<int>     *sigJetRecur_partonMatchID;
   vector<int>     *csFullJetRecur_partonMatchID;
   vector<int>     *csFullJetRecur_partonMatchDr;

   vector<double>  *sigJetRecurTauZcut_nSD;
   vector<double>  *sigJetRecurTauZcut_zSD;
   vector<double>  *sigJetRecurTau_nSD;
   vector<double>  *sigJetRecurTau_partonMatchDr;
   vector<double>  *sigJetRecurTau_zSD;
   vector<vector<double> > *sigJetRecurTauZcut_dr12;
   vector<vector<double> > *sigJetRecurTauZcut_erad;
   vector<vector<double> > *sigJetRecurTauZcut_jetpt;
   vector<vector<double> > *sigJetRecurTauZcut_jeteta;
   vector<vector<double> > *sigJetRecurTauZcut_logdr12;
   vector<vector<double> > *sigJetRecurTauZcut_logztheta;
   vector<vector<double> > *sigJetRecurTauZcut_tf;
   vector<vector<double> > *sigJetRecurTauZcut_tfe;
   vector<vector<double> > *sigJetRecurTauZcut_tfe2;
   vector<vector<double> > *sigJetRecurTauZcut_z;
   vector<vector<double> > *sigJetRecurTauZcut_kt;
   vector<vector<double> > *sigJetRecurTauZcut_droppedBFS;
   vector<vector<double> > *sigJetRecurTauZcut_droppedTfBFS;
   vector<vector<double> > *sigJetRecurTauZcut_droppedTfeBFS;
   vector<vector<double> > *sigJetRecurTauZcut_droppedPts;
   vector<vector<double> > *sigJetRecurTauZcut_droppedKts;
   vector<vector<double> > *sigJetRecurTauZcut_droppedLog1drBFS;
   vector<double> *sigJetRecurTauZcut_tau21;
   vector<double> *sigJetRecurTauZcut_tau32;

   vector<vector<double> > *sigJetRecurTau_dr12;
   vector<vector<double> > *sigJetRecurTau_erad;
   vector<vector<double> > *sigJetRecurTau_jetpt;
   vector<vector<double> > *sigJetRecurTau_jeteta;
   vector<vector<double> > *sigJetRecurTau_logdr12;
   vector<vector<double> > *sigJetRecurTau_logztheta;
   vector<vector<double> > *sigJetRecurTau_tf;
   vector<vector<double> > *sigJetRecurTau_tfe;
   vector<vector<double> > *sigJetRecurTau_tfe2;
   vector<vector<double> > *sigJetRecurTau_kt;
   vector<vector<double> > *sigJetRecurTau_z;
   vector<vector<double> > *sigJetRecurTau_droppedBFS;
   vector<int>     *sigJetRecurTau_partonMatchID;

   vector<vector<double> > *csFullJetRecurZcut_dr12;
   vector<vector<double> > *csFullJetRecurZcut_erad;
   vector<vector<double> > *csFullJetRecurZcut_jetpt;
   vector<vector<double> > *csFullJetRecurZcut_jeteta;
   vector<vector<double> > *csFullJetRecurZcut_logdr12;
   vector<vector<double> > *csFullJetRecurZcut_logztheta;
   vector<vector<double> > *csFullJetRecurZcut_tf;
   vector<vector<double> > *csFullJetRecurZcut_tfe;
   vector<vector<double> > *csFullJetRecurZcut_tfe2;
   vector<vector<double> > *csFullJetRecurZcut_z;
   vector<vector<double> > *csFullJetRecurZcut_kt;
   vector<vector<double> > *csFullJetRecurZcut_droppedBFS;
   vector<vector<double> > *csFullJetRecurZcut_droppedTfBFS;
   vector<vector<double> > *csFullJetRecurZcut_droppedTfeBFS;
   vector<vector<double> > *csFullJetRecurZcut_droppedPts;
   vector<vector<double> > *csFullJetRecurZcut_droppedKts;
   vector<vector<double> > *csFullJetRecurZcut_droppedLog1drBFS;

   vector<vector<double> > *csFullJetRecurTauZcut_dr12;
   vector<vector<double> > *csFullJetRecurTauZcut_erad;
   vector<vector<double> > *csFullJetRecurTauZcut_jetpt;
   vector<vector<double> > *csFullJetRecurTauZcut_jeteta;
   vector<vector<double> > *csFullJetRecurTauZcut_logdr12;
   vector<vector<double> > *csFullJetRecurTauZcut_logztheta;
   vector<vector<double> > *csFullJetRecurTauZcut_tf;
   vector<vector<double> > *csFullJetRecurTauZcut_tfe;
   vector<vector<double> > *csFullJetRecurTauZcut_tfe2;
   vector<vector<double> > *csFullJetRecurTauZcut_z;
   vector<vector<double> > *csFullJetRecurTauZcut_kt;
   vector<vector<double> > *csFullJetRecurTauZcut_droppedBFS;
   vector<vector<double> > *csFullJetRecurTauZcut_droppedTfBFS;
   vector<vector<double> > *csFullJetRecurTauZcut_droppedTfeBFS;
   vector<vector<double> > *csFullJetRecurTauZcut_droppedPts;
   vector<vector<double> > *csFullJetRecurTauZcut_droppedKts;
   vector<vector<double> > *csFullJetRecurTauZcut_droppedLog1drBFS;

   vector<vector<double> > *csFullJetRecur_dr12;
   vector<vector<double> > *csFullJetRecur_erad;
   vector<vector<double> > *csFullJetRecur_jetpt;
   vector<vector<double> > *csFullJetRecur_jeteta;
   vector<vector<double> > *csFullJetRecur_logdr12;
   vector<vector<double> > *csFullJetRecur_logztheta;
   vector<vector<double> > *csFullJetRecur_tf;
   vector<vector<double> > *csFullJetRecur_tfe;
   vector<vector<double> > *csFullJetRecur_tfe2;
   vector<vector<double> > *csFullJetRecur_z;
   vector<vector<double> > *csFullJetRecur_kt;

   vector<vector<double> > *csFullJetRecurTau_dr12;
   vector<vector<double> > *csFullJetRecurTau_erad;
   vector<vector<double> > *csFullJetRecurTau_jetpt;
   vector<vector<double> > *csFullJetRecurTau_jeteta;
   vector<vector<double> > *csFullJetRecurTau_logdr12;
   vector<vector<double> > *csFullJetRecurTau_logztheta;
   vector<vector<double> > *csFullJetRecurTau_tf;
   vector<vector<double> > *csFullJetRecurTau_tfe;
   vector<vector<double> > *csFullJetRecurTau_tfe2;
   vector<vector<double> > *csFullJetRecurTau_z;
   vector<vector<double> > *csFullJetRecurTau_kt;

   // List of branches
   TBranch        *b_eventWeight;   //!
   TBranch        *b_partonsPt;   //!
   TBranch        *b_partonsEta;   //!
   TBranch        *b_partonsPhi;   //!
   TBranch        *b_partonsM;   //!
   TBranch        *b_partonsPDG;   //!
   TBranch        *b_partonsFirstSplitPt;   //!
   TBranch        *b_partonsFirstSplitEta;   //!
   TBranch        *b_partonsFirstSplitPhi;   //!
   TBranch        *b_partonsFirstSplitM;   //!
   TBranch        *b_partonsFirstSplitPDG;   //!
   TBranch        *b_drsplit;   //!
   TBranch        *b_ktsplit;   //!
   TBranch        *b_tfsplit;   //!
   TBranch        *b_tfesplit;   //!
   TBranch        *b_partonWd1Wod2pt;//!
   TBranch        *b_partonWod1Wd2pt;
   TBranch        *b_partonWd1d2pt;
   TBranch        *b_partonWod1d2pt;
   TBranch        *b_zgsplit;   //!
   TBranch        *b_sigJetPt;   //!
   TBranch        *b_sigJetEta;   //!
   TBranch        *b_sigJetPhi;   //!
   TBranch        *b_sigJetM;   //!
   TBranch        *b_sigJetArea;   //!

   TBranch        *b_drsplits;
   TBranch        *b_tfsplits;
   TBranch        *b_tfesplits;

   TBranch        *b_sigJetRecurZcut_nSD;   //!
   TBranch        *b_sigJetRecurZcut_zSD;   //!
   TBranch        *b_sigJetRecur_nSD;   //!
   TBranch        *b_sigJetRecur_partonMatchDr;   //!
   TBranch        *b_sigJetRecur_zSD;   //!
   TBranch        *b_sigJetRecurZcut_dr12;   //!
   TBranch        *b_sigJetRecurZcut_erad;   //!
   TBranch        *b_sigJetRecurZcut_jetpt;   //!
   TBranch        *b_sigJetRecurZcut_jeteta;   //!
   TBranch        *b_sigJetRecurZcut_logdr12;   //!
   TBranch        *b_sigJetRecurZcut_logztheta;   //!
   TBranch        *b_sigJetRecurZcut_tf;   //!
   TBranch        *b_sigJetRecurZcut_tfe;   //!
   TBranch        *b_sigJetRecurZcut_tfe2;   //!
   TBranch        *b_sigJetRecurZcut_kt;   //!
   TBranch        *b_sigJetRecurZcut_z;   //!
   TBranch        *b_sigJetRecurZcut_droppedBFS;
   TBranch        *b_sigJetRecurZcut_droppedTfBFS;
   TBranch        *b_sigJetRecurZcut_droppedTfeBFS;
   TBranch        *b_sigJetRecurZcut_droppedPts;
   TBranch        *b_sigJetRecurZcut_droppedKts;
   TBranch        *b_sigJetRecurZcut_droppedLog1drBFS;

   TBranch        *b_sigJetRecurZcut_tau21;
   TBranch        *b_sigJetRecurZcut_tau32;
   TBranch        *b_sigJetRecur_dr12;   //!
   TBranch        *b_sigJetRecur_erad;   //!
   TBranch        *b_sigJetRecur_jetpt;   //!
   TBranch        *b_sigJetRecur_jeteta;   //!
   TBranch        *b_sigJetRecur_logdr12;   //!
   TBranch        *b_sigJetRecur_logztheta;   //!
   TBranch        *b_sigJetRecur_tf;   //!
   TBranch        *b_sigJetRecur_tfe;   //!
   TBranch        *b_sigJetRecur_tfe2;   //!
   TBranch        *b_sigJetRecur_kt;   //!
   TBranch        *b_sigJetRecur_z;   //!
   TBranch        *b_sigJetRecur_droppedBFS;
   TBranch        *b_sigJetRecur_partonMatchID;   //!
   TBranch        *b_csFullJetRecur_partonMatchID;   //!
   TBranch        *b_csFullJetRecur_partonMatchDr;   //!

   TBranch        *b_sigJetRecurTauZcut_nSD;   //!
   TBranch        *b_sigJetRecurTauZcut_zSD;   //!
   TBranch        *b_sigJetRecurTau_nSD;   //!
   TBranch        *b_sigJetRecurTau_partonMatchDr;   //!
   TBranch        *b_sigJetRecurTau_zSD;   //!
   TBranch        *b_sigJetRecurTauZcut_dr12;   //!
   TBranch        *b_sigJetRecurTauZcut_erad;   //!
   TBranch        *b_sigJetRecurTauZcut_jetpt;   //!
   TBranch        *b_sigJetRecurTauZcut_jeteta;   //!
   TBranch        *b_sigJetRecurTauZcut_logdr12;   //!
   TBranch        *b_sigJetRecurTauZcut_logztheta;   //!
   TBranch        *b_sigJetRecurTauZcut_tf;   //!
   TBranch        *b_sigJetRecurTauZcut_tfe;   //!
   TBranch        *b_sigJetRecurTauZcut_tfe2;   //!
   TBranch        *b_sigJetRecurTauZcut_z;   //!
   TBranch        *b_sigJetRecurTauZcut_kt;   //!
   TBranch        *b_sigJetRecurTauZcut_droppedBFS;
   TBranch        *b_sigJetRecurTauZcut_droppedTfBFS;
   TBranch        *b_sigJetRecurTauZcut_droppedTfeBFS;
   TBranch        *b_sigJetRecurTauZcut_droppedPts;
   TBranch        *b_sigJetRecurTauZcut_droppedKts;
   TBranch        *b_sigJetRecurTauZcut_droppedLog1drBFS;

   TBranch        *b_sigJetRecurTauZcut_tau21;
   TBranch        *b_sigJetRecurTauZcut_tau32;
   TBranch        *b_sigJetRecurTau_dr12;   //!
   TBranch        *b_sigJetRecurTau_erad;   //!
   TBranch        *b_sigJetRecurTau_jetpt;   //!
   TBranch        *b_sigJetRecurTau_jeteta;   //!
   TBranch        *b_sigJetRecurTau_logdr12;   //!
   TBranch        *b_sigJetRecurTau_logztheta;   //!
   TBranch        *b_sigJetRecurTau_tf;   //!
   TBranch        *b_sigJetRecurTau_tfe;   //!
   TBranch        *b_sigJetRecurTau_tfe2;   //!
   TBranch        *b_sigJetRecurTau_kt;   //!
   TBranch        *b_sigJetRecurTau_z;   //!
   TBranch        *b_sigJetRecurTau_droppedBFS;
   TBranch        *b_sigJetRecurTau_partonMatchID;   //!

   TBranch         *b_sigJetDYE_z; 
   TBranch         *b_sigJetDYE_dr12; 
   TBranch         *b_sigJetDYE_kt; 
   TBranch         *b_sigJetDYE_tau21; 
   TBranch         *b_sigJetDYE_tau32; 
   TBranch         *b_sigJetDYE_tf; 
   TBranch         *b_sigJetDYE_tfe; 
   TBranch         *b_sigJetkappaDYE; 
   TBranch         *b_sigJetDYE_drBranches;
   TBranch         *b_sigJetTauDYE_z; 
   TBranch         *b_sigJetTauDYE_dr12; 
   TBranch         *b_sigJetTauDYE_kt; 
   TBranch         *b_sigJetTauDYE_tau21; 
   TBranch         *b_sigJetTauDYE_tau32; 
   TBranch         *b_sigJetTauDYE_tf; 
   TBranch         *b_sigJetTauDYE_tfe; 
   TBranch         *b_sigJetTaukappaDYE; 
   TBranch         *b_sigJetTauDYE_drBranches;

   TBranch         *b_csFullJetRecurZcut_dr12;
   TBranch         *b_csFullJetRecurZcut_erad;
   TBranch         *b_csFullJetRecurZcut_jetpt;
   TBranch         *b_csFullJetRecurZcut_jeteta;
   TBranch         *b_csFullJetRecurZcut_logdr12;
   TBranch         *b_csFullJetRecurZcut_logztheta;
   TBranch         *b_csFullJetRecurZcut_tf;
   TBranch         *b_csFullJetRecurZcut_tfe;
   TBranch         *b_csFullJetRecurZcut_tfe2;
   TBranch         *b_csFullJetRecurZcut_z;
   TBranch         *b_csFullJetRecurZcut_kt;
   TBranch         *b_csFullJetRecurZcut_droppedBFS;
   TBranch         *b_csFullJetRecurZcut_droppedTfBFS;
   TBranch         *b_csFullJetRecurZcut_droppedTfeBFS;
   TBranch         *b_csFullJetRecurZcut_droppedPts;
   TBranch         *b_csFullJetRecurZcut_droppedKts;
   TBranch         *b_csFullJetRecurZcut_droppedLog1drBFS;

   TBranch         *b_csFullJetRecurTauZcut_dr12;
   TBranch         *b_csFullJetRecurTauZcut_erad;
   TBranch         *b_csFullJetRecurTauZcut_jetpt;
   TBranch         *b_csFullJetRecurTauZcut_jeteta;
   TBranch         *b_csFullJetRecurTauZcut_logdr12;
   TBranch         *b_csFullJetRecurTauZcut_logztheta;
   TBranch         *b_csFullJetRecurTauZcut_tf;
   TBranch         *b_csFullJetRecurTauZcut_tfe;
   TBranch         *b_csFullJetRecurTauZcut_tfe2;
   TBranch         *b_csFullJetRecurTauZcut_z;
   TBranch         *b_csFullJetRecurTauZcut_kt;
   TBranch         *b_csFullJetRecurTauZcut_droppedBFS;
   TBranch         *b_csFullJetRecurTauZcut_droppedTfBFS;
   TBranch         *b_csFullJetRecurTauZcut_droppedTfeBFS;
   TBranch         *b_csFullJetRecurTauZcut_droppedPts;
   TBranch         *b_csFullJetRecurTauZcut_droppedKts;
   TBranch         *b_csFullJetRecurTauZcut_droppedLog1drBFS;

   TBranch         *b_csFullJetRecur_dr12;
   TBranch         *b_csFullJetRecur_erad;
   TBranch         *b_csFullJetRecur_jetpt;
   TBranch         *b_csFullJetRecur_jeteta;
   TBranch         *b_csFullJetRecur_logdr12;
   TBranch         *b_csFullJetRecur_logztheta;
   TBranch         *b_csFullJetRecur_tf;
   TBranch         *b_csFullJetRecur_tfe;
   TBranch         *b_csFullJetRecur_tfe2;
   TBranch         *b_csFullJetRecur_z;
   TBranch         *b_csFullJetRecur_kt;

   TBranch         *b_csFullJetRecurTau_dr12;
   TBranch         *b_csFullJetRecurTau_erad;
   TBranch         *b_csFullJetRecurTau_jetpt;
   TBranch         *b_csFullJetRecurTau_jeteta;
   TBranch         *b_csFullJetRecurTau_logdr12;
   TBranch         *b_csFullJetRecurTau_logztheta;
   TBranch         *b_csFullJetRecurTau_tf;
   TBranch         *b_csFullJetRecurTau_tfe;
   TBranch         *b_csFullJetRecurTau_tfe2;
   TBranch         *b_csFullJetRecurTau_z;
   TBranch         *b_csFullJetRecurTau_kt;

   jetTreeSig(TTree *tree=0);
   virtual ~jetTreeSig();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
  // virtual TTree*   CopyTree(const char* selection);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Bool_t reduced, TString generator, TString duration, TString fin, double R);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef jetTreeSig_cxx
jetTreeSig::jetTreeSig(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) { //new TFile("Johanna/Jewel/JetToyHITimeClusHistosJewel_vac_"+generator+"_"+duration+".root","RECREATE");
      std::cout<<"tree 0"<<std::endl;
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/dcache/alice/jlomker/PYTHIA/LHC_13p6TeV/analysisResults/JetToyHIResultTimeClusBkg_13p6_200.root");
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/dcache/alice/jlomker/Jewel/vac/analysisResults/JetToyHIResultTimeClusJewel_200.root");
      if (!f || !f->IsOpen()) {
	 std::cout<<"file not open"<<std::endl;
         //f = new TFile("/dcache/alice/jlomker/PYTHIA/LHC_13p6TeV/analysisResults/JetToyHIResultTimeClusBkg_13p6_200.root");
         f = new TFile("/dcache/alice/jlomker/Jewel/vac/analysisResults/JetToyHIResultTimeClusJewel_200.root");
	 //f = new TFile("JetToyHIResultTimeClusBkg_13p6.root");
      }
      f->GetObject("jetTreeSig",tree);

   }
   Init(tree);
}

jetTreeSig::~jetTreeSig()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t jetTreeSig::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t jetTreeSig::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void jetTreeSig::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   eventWeight = 0;
   partonsPt = 0;
   partonsEta = 0;
   partonsPhi = 0;
   partonsM = 0;
   partonsPDG = 0;
   partonsFirstSplitPt = 0;
   partonsFirstSplitEta = 0;
   partonsFirstSplitPhi = 0;
   partonsFirstSplitM = 0;
   partonsFirstSplitPDG = 0;
   drsplit = 0;
   ktsplit = 0;
   tfsplit = 0;
   tfesplit = 0;
   partonWd1Wod2pt = 0;
   partonWod1Wd2pt = 0;
   partonWd1d2pt = 0;
   partonWod1d2pt = 0;
   zgsplit = 0;
   drsplits = 0;
   tfsplits = 0;
   tfesplits = 0;
   sigJetPt = 0;
   //jetCollectionSigDYEPt = 0;
   sigJetEta = 0;
   sigJetPhi = 0;
   sigJetM = 0;
   sigJetArea = 0;
   // Reclustering effects
   sigJetRecurZcut_nSD = 0;
   sigJetRecurZcut_zSD = 0;
   sigJetRecur_nSD = 0;
   sigJetRecur_partonMatchDr = 0;
   sigJetRecur_zSD = 0;
   sigJetRecurZcut_dr12 = 0;
   sigJetRecurZcut_erad = 0;
   sigJetRecurZcut_jetpt = 0;
   sigJetRecurZcut_jeteta = 0;
   sigJetRecurZcut_logdr12 = 0;
   sigJetRecurZcut_logztheta = 0;
   sigJetRecurZcut_tf = 0;
   sigJetRecurZcut_tfe = 0;
   sigJetRecurZcut_tfe2 = 0;
   sigJetRecurZcut_kt = 0;
   sigJetRecurZcut_z = 0;
   sigJetRecurZcut_droppedBFS = 0;
   sigJetRecurZcut_droppedTfBFS = 0;
   sigJetRecurZcut_droppedTfeBFS = 0;
   sigJetRecurZcut_droppedPts = 0;
   sigJetRecurZcut_droppedKts = 0;
   sigJetRecurZcut_droppedLog1drBFS = 0;
   
   sigJetRecurZcut_tau21 = 0;
   sigJetRecurZcut_tau32 = 0;
   sigJetRecur_dr12 = 0;
   sigJetRecur_erad = 0;
   sigJetRecur_jetpt = 0;
   sigJetRecur_jeteta = 0;
   sigJetRecur_logdr12 = 0;
   sigJetRecur_logztheta = 0;
   sigJetRecur_tf = 0;
   sigJetRecur_tfe = 0;
   sigJetRecur_tfe2 = 0;
   sigJetRecur_kt = 0;
   sigJetRecur_z = 0;
   sigJetRecur_droppedBFS = 0;
   sigJetRecur_partonMatchID = 0;
   csFullJetRecur_partonMatchID = 0;
   csFullJetRecur_partonMatchDr=0;
   
   sigJetRecurTauZcut_nSD = 0;
   sigJetRecurTauZcut_zSD = 0;
   sigJetRecurTau_nSD = 0;
   sigJetRecurTau_partonMatchDr = 0;
   sigJetRecurTau_zSD = 0;
   sigJetRecurTauZcut_dr12 = 0;
   sigJetRecurTauZcut_erad = 0;
   sigJetRecurTauZcut_jetpt = 0;
   sigJetRecurTauZcut_jeteta = 0;
   sigJetRecurTauZcut_logdr12 = 0;
   sigJetRecurTauZcut_logztheta = 0;
   sigJetRecurTauZcut_tf = 0;
   sigJetRecurTauZcut_tfe = 0;
   sigJetRecurTauZcut_tfe2 = 0;
   sigJetRecurTauZcut_z = 0;
   sigJetRecurTauZcut_kt = 0;
   sigJetRecurTauZcut_droppedBFS = 0;
   sigJetRecurTauZcut_droppedTfBFS = 0;
   sigJetRecurTauZcut_droppedTfeBFS = 0;
   sigJetRecurTauZcut_droppedPts = 0;
   sigJetRecurTauZcut_droppedKts = 0;
   sigJetRecurTauZcut_droppedLog1drBFS = 0;
   sigJetRecurTauZcut_tau21 = 0;
   sigJetRecurTauZcut_tau32 = 0;
   sigJetRecurTau_dr12 = 0;
   sigJetRecurTau_erad = 0;
   sigJetRecurTau_jetpt = 0;
   sigJetRecurTau_jeteta = 0;
   sigJetRecurTau_logdr12 = 0;
   sigJetRecurTau_logztheta = 0;
   sigJetRecurTau_tf = 0;
   sigJetRecurTau_tfe = 0;
   sigJetRecurTau_tfe2 = 0;
   sigJetRecurTau_kt = 0;
   sigJetRecurTau_z = 0;
   sigJetRecurTau_droppedBFS = 0;

   sigJetDYE_z = 0;         //zg of groomed jets
   sigJetDYE_dr12 = 0;       //distance between the two subjet
   sigJetDYE_kt = 0;         //kt of groomed jets
   sigJetDYE_tau21 = 0;      //n-subjettiness ratio 21
   sigJetDYE_tau32 = 0;      //n-subjettiness ratio 32
   sigJetDYE_tf = 0;     	// formation time tf = 2 / (z*(1-z)*pt*theta^2) (soft limit)
   sigJetDYE_tfe = 0;     	// formation time tf = 1 / 2*(z*(1-z)*E*(1-cos(theta))) (hard limit)
   sigJetkappaDYE = 0;      //kappa of groomed jets
   sigJetDYE_drBranches = 0;
   sigJetTauDYE_z = 0;
   sigJetTauDYE_dr12 = 0;
   sigJetTauDYE_kt = 0;
   sigJetTauDYE_tau21 = 0;
   sigJetTauDYE_tau32 = 0;
   sigJetTauDYE_tf = 0;
   sigJetTauDYE_tfe = 0;
   sigJetTaukappaDYE = 0;
   sigJetTauDYE_drBranches = 0;

   csFullJetRecurZcut_dr12 = 0;
   csFullJetRecurZcut_erad = 0;
   csFullJetRecurZcut_jetpt = 0;
   csFullJetRecurZcut_jeteta = 0;
   csFullJetRecurZcut_logdr12 = 0;
   csFullJetRecurZcut_logztheta = 0;
   csFullJetRecurZcut_tf = 0;
   csFullJetRecurZcut_tfe = 0;
   csFullJetRecurZcut_tfe2 = 0;
   csFullJetRecurZcut_z = 0;
   csFullJetRecurZcut_kt = 0;
   csFullJetRecurZcut_droppedBFS = 0;
   csFullJetRecurZcut_droppedTfBFS = 0;
   csFullJetRecurZcut_droppedTfeBFS = 0;
   csFullJetRecurZcut_droppedPts = 0;
   csFullJetRecurZcut_droppedKts = 0;
   csFullJetRecurZcut_droppedLog1drBFS = 0;

   csFullJetRecurTauZcut_jetpt = 0;
   csFullJetRecurTauZcut_dr12 = 0;
   csFullJetRecurTauZcut_erad = 0;
   csFullJetRecurTauZcut_jetpt = 0;
   csFullJetRecurTauZcut_jeteta = 0;
   csFullJetRecurTauZcut_logdr12 = 0;
   csFullJetRecurTauZcut_logztheta = 0;
   csFullJetRecurTauZcut_tf = 0;
   csFullJetRecurTauZcut_tfe = 0;
   csFullJetRecurTauZcut_tfe2 = 0;
   csFullJetRecurTauZcut_z = 0;
   csFullJetRecurTauZcut_kt = 0;
   csFullJetRecurTauZcut_droppedBFS = 0;
   csFullJetRecurTauZcut_droppedTfBFS = 0;
   csFullJetRecurTauZcut_droppedTfeBFS = 0;
   csFullJetRecurTauZcut_droppedPts = 0;
   csFullJetRecurTauZcut_droppedKts = 0;
   csFullJetRecurTauZcut_droppedLog1drBFS = 0;

   csFullJetRecur_dr12 = 0;
   csFullJetRecur_erad = 0;
   csFullJetRecur_jetpt = 0;
   csFullJetRecur_jeteta = 0;
   csFullJetRecur_logdr12 = 0;
   csFullJetRecur_logztheta = 0;
   csFullJetRecur_tf = 0;
   csFullJetRecur_tfe = 0;
   csFullJetRecur_tfe2 = 0;
   csFullJetRecur_z = 0;
   csFullJetRecur_kt = 0;

   csFullJetRecurTau_jetpt = 0;
   csFullJetRecurTau_dr12 = 0;
   csFullJetRecurTau_erad = 0;
   csFullJetRecurTau_jetpt = 0;
   csFullJetRecurTau_jeteta = 0;
   csFullJetRecurTau_logdr12 = 0;
   csFullJetRecurTau_logztheta = 0;
   csFullJetRecurTau_tf = 0;
   csFullJetRecurTau_tfe = 0;
   csFullJetRecurTau_tfe2 = 0;
   csFullJetRecurTau_z = 0;
   csFullJetRecurTau_kt = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventWeight", &eventWeight, &b_eventWeight);
   // True parton info
   fChain->SetBranchAddress("partonsPt", &partonsPt, &b_partonsPt);
   fChain->SetBranchAddress("partonsEta", &partonsEta, &b_partonsEta);
   fChain->SetBranchAddress("partonsPhi", &partonsPhi, &b_partonsPhi);
   fChain->SetBranchAddress("partonsM", &partonsM, &b_partonsM);
   fChain->SetBranchAddress("partonsPDG", &partonsPDG, &b_partonsPDG);
   fChain->SetBranchAddress("partonsFirstSplitPt", &partonsFirstSplitPt, &b_partonsFirstSplitPt);
   fChain->SetBranchAddress("partonsFirstSplitEta", &partonsFirstSplitEta, &b_partonsFirstSplitEta);
   fChain->SetBranchAddress("partonsFirstSplitPhi", &partonsFirstSplitPhi, &b_partonsFirstSplitPhi);
   fChain->SetBranchAddress("partonsFirstSplitM", &partonsFirstSplitM, &b_partonsFirstSplitM);
   fChain->SetBranchAddress("partonsFirstSplitPDG", &partonsFirstSplitPDG, &b_partonsFirstSplitPDG);
   fChain->SetBranchAddress("drsplit", &drsplit, &b_drsplit);
   fChain->SetBranchAddress("ktsplit", &ktsplit, &b_ktsplit);
   fChain->SetBranchAddress("tfsplit", &tfsplit, &b_tfsplit);
   fChain->SetBranchAddress("tfesplit", &tfesplit, &b_tfesplit);
   fChain->SetBranchAddress("partonWd1Wod2pt", &partonWd1Wod2pt, &b_partonWd1Wod2pt);
   fChain->SetBranchAddress("partonWod1Wd2pt", &partonWod1Wd2pt, &b_partonWod1Wd2pt);
   fChain->SetBranchAddress("partonWd1d2pt", &partonWd1d2pt, &b_partonWd1d2pt);
   fChain->SetBranchAddress("partonWod1d2pt", &partonWod1d2pt, &b_partonWod1d2pt);
   fChain->SetBranchAddress("zgsplit", &zgsplit, &b_zgsplit);
   fChain->SetBranchAddress("drsplits", &drsplits, &b_drsplits);
   fChain->SetBranchAddress("tfsplits", &tfsplits, &b_tfsplits);
   fChain->SetBranchAddress("tfesplits", &tfesplits, &b_tfesplits);
   // True jet info (anti-kt clustering)
   fChain->SetBranchAddress("sigJetPt", &sigJetPt, &b_sigJetPt);
   fChain->SetBranchAddress("sigJetEta", &sigJetEta, &b_sigJetEta);
   fChain->SetBranchAddress("sigJetPhi", &sigJetPhi, &b_sigJetPhi);
   fChain->SetBranchAddress("sigJetM", &sigJetM, &b_sigJetM);
   fChain->SetBranchAddress("sigJetArea", &sigJetArea, &b_sigJetArea);
   //pTDSig and widthSig (from lambda angularity)
   // Declustering info C/A with grooming
   fChain->SetBranchAddress("sigJetRecurZcut_nSD", &sigJetRecurZcut_nSD, &b_sigJetRecurZcut_nSD);
   fChain->SetBranchAddress("sigJetRecurZcut_zSD", &sigJetRecurZcut_zSD, &b_sigJetRecurZcut_zSD);
   fChain->SetBranchAddress("sigJetRecur_partonMatchDr", &sigJetRecur_partonMatchDr, &b_sigJetRecur_partonMatchDr);
   fChain->SetBranchAddress("sigJetRecurZcut_dr12", &sigJetRecurZcut_dr12, &b_sigJetRecurZcut_dr12);
   fChain->SetBranchAddress("sigJetRecurZcut_erad", &sigJetRecurZcut_erad, &b_sigJetRecurZcut_erad);
   fChain->SetBranchAddress("sigJetRecurZcut_jetpt", &sigJetRecurZcut_jetpt, &b_sigJetRecurZcut_jetpt);
   fChain->SetBranchAddress("sigJetRecurZcut_jeteta", &sigJetRecurZcut_jeteta, &b_sigJetRecurZcut_jeteta);
   fChain->SetBranchAddress("sigJetRecurZcut_logdr12", &sigJetRecurZcut_logdr12, &b_sigJetRecurZcut_logdr12);
   fChain->SetBranchAddress("sigJetRecurZcut_logztheta", &sigJetRecurZcut_logztheta, &b_sigJetRecurZcut_logztheta);
   fChain->SetBranchAddress("sigJetRecurZcut_tf", &sigJetRecurZcut_tf, &b_sigJetRecurZcut_tf);
   fChain->SetBranchAddress("sigJetRecurZcut_tfe", &sigJetRecurZcut_tfe, &b_sigJetRecurZcut_tfe);
   fChain->SetBranchAddress("sigJetRecurZcut_tfe2", &sigJetRecurZcut_tfe2, &b_sigJetRecurZcut_tfe2);
   fChain->SetBranchAddress("sigJetRecurZcut_kt", &sigJetRecurZcut_kt, &b_sigJetRecurZcut_kt);
   fChain->SetBranchAddress("sigJetRecurZcut_z", &sigJetRecurZcut_z, &b_sigJetRecurZcut_z);
   fChain->SetBranchAddress("sigJetRecurZcut_droppedBFS", &sigJetRecurZcut_droppedBFS, &b_sigJetRecurZcut_droppedBFS);
   fChain->SetBranchAddress("sigJetRecurZcut_droppedTfBFS", &sigJetRecurZcut_droppedTfBFS, &b_sigJetRecurZcut_droppedTfBFS);
   fChain->SetBranchAddress("sigJetRecurZcut_droppedTfeBFS", &sigJetRecurZcut_droppedTfeBFS, &b_sigJetRecurZcut_droppedTfeBFS);
   fChain->SetBranchAddress("sigJetRecurZcut_droppedPts", &sigJetRecurZcut_droppedPts, &b_sigJetRecurZcut_droppedPts);
   fChain->SetBranchAddress("sigJetRecurZcut_droppedKts", &sigJetRecurZcut_droppedKts, &b_sigJetRecurZcut_droppedKts);
   fChain->SetBranchAddress("sigJetRecurZcut_droppedLog1drBFS", &sigJetRecurZcut_droppedLog1drBFS, &b_sigJetRecurZcut_droppedLog1drBFS);
   fChain->SetBranchAddress("sigJetRecurZcut_tau21", &sigJetRecurZcut_tau21, &b_sigJetRecurZcut_tau21);
   fChain->SetBranchAddress("sigJetRecurZcut_tau32", &sigJetRecurZcut_tau32, &b_sigJetRecurTauZcut_tau32);
   // Declustering info C/A no grooming
   fChain->SetBranchAddress("sigJetRecur_erad", &sigJetRecur_erad, &b_sigJetRecur_erad);
   fChain->SetBranchAddress("sigJetRecur_jetpt", &sigJetRecur_jetpt, &b_sigJetRecur_jetpt);
   fChain->SetBranchAddress("sigJetRecur_jeteta", &sigJetRecur_jeteta, &b_sigJetRecur_jeteta);
   fChain->SetBranchAddress("sigJetRecur_logdr12", &sigJetRecur_logdr12, &b_sigJetRecur_logdr12);
   fChain->SetBranchAddress("sigJetRecur_dr12", &sigJetRecur_dr12, &b_sigJetRecur_dr12);
   fChain->SetBranchAddress("sigJetRecur_logztheta", &sigJetRecur_logztheta, &b_sigJetRecur_logztheta);
   fChain->SetBranchAddress("sigJetRecur_tf", &sigJetRecur_tf, &b_sigJetRecur_tf);
   fChain->SetBranchAddress("sigJetRecur_tfe", &sigJetRecur_tfe, &b_sigJetRecur_tfe);
   fChain->SetBranchAddress("sigJetRecur_tfe2", &sigJetRecur_tfe2, &b_sigJetRecur_tfe2);
   fChain->SetBranchAddress("sigJetRecur_kt", &sigJetRecur_kt, &b_sigJetRecur_kt);
   fChain->SetBranchAddress("sigJetRecur_nSD", &sigJetRecur_nSD, &b_sigJetRecur_nSD);
   fChain->SetBranchAddress("sigJetRecur_zSD", &sigJetRecur_zSD, &b_sigJetRecur_zSD);
   fChain->SetBranchAddress("sigJetRecur_z", &sigJetRecur_z, &b_sigJetRecur_z);
   fChain->SetBranchAddress("sigJetRecur_droppedBFS", &sigJetRecur_droppedBFS, &b_sigJetRecur_droppedBFS);
   fChain->SetBranchAddress("sigJetRecur_partonMatchID", &sigJetRecur_partonMatchID, &b_sigJetRecur_partonMatchID);
   // Declustering info tau with grooming
   fChain->SetBranchAddress("sigJetRecurTauZcut_nSD", &sigJetRecurTauZcut_nSD, &b_sigJetRecurTauZcut_nSD);
   fChain->SetBranchAddress("sigJetRecurTauZcut_zSD", &sigJetRecurTauZcut_zSD, &b_sigJetRecurTauZcut_zSD);
   //fChain->SetBranchAddress("sigJetRecurTauZcut_partonMatchDr", &sigJetRecurTau_partonMatchDr, &b_sigJetRecurTau_partonMatchDr);
   fChain->SetBranchAddress("sigJetRecurTauZcut_dr12", &sigJetRecurTauZcut_dr12, &b_sigJetRecurTauZcut_dr12);
   fChain->SetBranchAddress("sigJetRecurTauZcut_erad", &sigJetRecurTauZcut_erad, &b_sigJetRecurTauZcut_erad);
   fChain->SetBranchAddress("sigJetRecurTauZcut_jetpt", &sigJetRecurTauZcut_jetpt, &b_sigJetRecurTauZcut_jetpt);
   fChain->SetBranchAddress("sigJetRecurTauZcut_jeteta", &sigJetRecurTauZcut_jeteta, &b_sigJetRecurTauZcut_jeteta);
   fChain->SetBranchAddress("sigJetRecurTauZcut_logdr12", &sigJetRecurTauZcut_logdr12, &b_sigJetRecurTauZcut_logdr12);
   fChain->SetBranchAddress("sigJetRecurTauZcut_logztheta", &sigJetRecurTauZcut_logztheta, &b_sigJetRecurTauZcut_logztheta);
   fChain->SetBranchAddress("sigJetRecurTauZcut_tf", &sigJetRecurTauZcut_tf, &b_sigJetRecurTauZcut_tf);
   fChain->SetBranchAddress("sigJetRecurTauZcut_tfe", &sigJetRecurTauZcut_tfe, &b_sigJetRecurTauZcut_tfe);
   fChain->SetBranchAddress("sigJetRecurTauZcut_tfe2", &sigJetRecurTauZcut_tfe2, &b_sigJetRecurTauZcut_tfe2);
   fChain->SetBranchAddress("sigJetRecurTauZcut_z", &sigJetRecurTauZcut_z, &b_sigJetRecurTauZcut_z);
   fChain->SetBranchAddress("sigJetRecurTauZcut_kt", &sigJetRecurTauZcut_kt, &b_sigJetRecurTauZcut_kt);
   fChain->SetBranchAddress("sigJetRecurTauZcut_droppedBFS", &sigJetRecurTauZcut_droppedBFS, &b_sigJetRecurTauZcut_droppedBFS);
   fChain->SetBranchAddress("sigJetRecurTauZcut_droppedTfBFS", &sigJetRecurTauZcut_droppedTfBFS, &b_sigJetRecurTauZcut_droppedTfBFS);
   fChain->SetBranchAddress("sigJetRecurTauZcut_droppedTfeBFS", &sigJetRecurTauZcut_droppedTfeBFS, &b_sigJetRecurTauZcut_droppedTfeBFS);
   fChain->SetBranchAddress("sigJetRecurTauZcut_droppedPts", &sigJetRecurTauZcut_droppedPts, &b_sigJetRecurTauZcut_droppedPts);
   fChain->SetBranchAddress("sigJetRecurTauZcut_droppedKts", &sigJetRecurTauZcut_droppedKts, &b_sigJetRecurTauZcut_droppedKts);
   fChain->SetBranchAddress("sigJetRecurTauZcut_droppedLog1drBFS", &sigJetRecurTauZcut_droppedLog1drBFS, &b_sigJetRecurTauZcut_droppedLog1drBFS);
   fChain->SetBranchAddress("sigJetRecurTauZcut_tau21", &sigJetRecurTauZcut_tau21, &b_sigJetRecurTauZcut_tau21);
   fChain->SetBranchAddress("sigJetRecurTauZcut_tau32", &sigJetRecurTauZcut_tau32, &b_sigJetRecurTauZcut_tau32);

   // Declustering info tau no grooming
   fChain->SetBranchAddress("sigJetRecurTau_erad", &sigJetRecurTau_erad, &b_sigJetRecurTau_erad);
   fChain->SetBranchAddress("sigJetRecurTau_jetpt", &sigJetRecurTau_jetpt, &b_sigJetRecurTau_jetpt);
   fChain->SetBranchAddress("sigJetRecurTau_jeteta", &sigJetRecurTau_jeteta, &b_sigJetRecurTau_jeteta);
   fChain->SetBranchAddress("sigJetRecurTau_logdr12", &sigJetRecurTau_logdr12, &b_sigJetRecurTau_logdr12);
   fChain->SetBranchAddress("sigJetRecurTau_dr12", &sigJetRecurTau_dr12, &b_sigJetRecurTau_dr12);
   fChain->SetBranchAddress("sigJetRecurTau_logztheta", &sigJetRecurTau_logztheta, &b_sigJetRecurTau_logztheta);
   fChain->SetBranchAddress("sigJetRecurTau_tf", &sigJetRecurTau_tf, &b_sigJetRecurTau_tf);
   fChain->SetBranchAddress("sigJetRecurTau_tfe", &sigJetRecurTau_tfe, &b_sigJetRecurTau_tfe);
   fChain->SetBranchAddress("sigJetRecurTau_tfe2", &sigJetRecurTau_tfe2, &b_sigJetRecurTau_tfe2);
   fChain->SetBranchAddress("sigJetRecurTau_kt", &sigJetRecurTau_kt, &b_sigJetRecurTau_kt);
   fChain->SetBranchAddress("sigJetRecurTau_z", &sigJetRecurTau_z, &b_sigJetRecurTau_z);
   fChain->SetBranchAddress("sigJetRecurTau_droppedBFS", &sigJetRecurTau_droppedBFS, &b_sigJetRecurTau_droppedBFS);
   fChain->SetBranchAddress("sigJetRecurTau_nSD", &sigJetRecurTau_nSD, &b_sigJetRecurTau_nSD);
   fChain->SetBranchAddress("sigJetRecurTau_zSD", &sigJetRecurTau_zSD, &b_sigJetRecurTau_zSD);
   //fChain->SetBranchAddress("sigJetRecurTau_partonMatchID", &sigJetRecurTau_partonMatchID, &b_sigJetRecurTau_partonMatchID); 
   // Effect of bkg subtraction removed here.
   fChain->SetBranchAddress("sigJetDYE_z", &sigJetDYE_z, &b_sigJetDYE_z);
   fChain->SetBranchAddress("sigJetDYE_dr12", &sigJetDYE_dr12, &b_sigJetDYE_dr12);
   fChain->SetBranchAddress("sigJetDYE_kt", &sigJetDYE_kt, &b_sigJetDYE_kt);
   fChain->SetBranchAddress("sigJetDYE_tau21", &sigJetDYE_tau21, &b_sigJetDYE_tau21);
   fChain->SetBranchAddress("sigJetDYE_tau32", &sigJetDYE_tau32, &b_sigJetDYE_tau32);
   fChain->SetBranchAddress("sigJetDYE_tf", &sigJetDYE_tf, &b_sigJetDYE_tf);
   fChain->SetBranchAddress("sigJetDYE_tfe", &sigJetDYE_tfe, &b_sigJetDYE_tfe);
   fChain->SetBranchAddress("sigJetkappaDYE", &sigJetkappaDYE, &b_sigJetkappaDYE);
   fChain->SetBranchAddress("sigJetDYE_drBranches", &sigJetDYE_drBranches, &b_sigJetDYE_drBranches);
   fChain->SetBranchAddress("sigJetTauDYE_z", &sigJetTauDYE_z, &b_sigJetTauDYE_z);
   fChain->SetBranchAddress("sigJetTauDYE_dr12", &sigJetTauDYE_dr12, &b_sigJetTauDYE_dr12);   
   fChain->SetBranchAddress("sigJetTauDYE_kt", &sigJetTauDYE_kt, &b_sigJetTauDYE_kt);
   fChain->SetBranchAddress("sigJetTauDYE_tau21", &sigJetTauDYE_tau21, &b_sigJetTauDYE_tau21);
   fChain->SetBranchAddress("sigJetTauDYE_tau32", &sigJetTauDYE_tau32, &b_sigJetTauDYE_tau32);
   fChain->SetBranchAddress("sigJetTauDYE_tf", &sigJetTauDYE_tf, &b_sigJetTauDYE_tf);
   fChain->SetBranchAddress("sigJetTauDYE_tfe", &sigJetTauDYE_tfe, &b_sigJetTauDYE_tfe);
   fChain->SetBranchAddress("sigJetTaukappaDYE", &sigJetTaukappaDYE, &b_sigJetTaukappaDYE);
   fChain->SetBranchAddress("sigJetTauDYE_drBranches", &sigJetTauDYE_drBranches, &b_sigJetTauDYE_drBranches);
   //constituent subtracted ones
   fChain->SetBranchAddress("csFullJetRecurZcut_dr12", &csFullJetRecurZcut_dr12, &b_csFullJetRecurZcut_dr12);
   fChain->SetBranchAddress("csFullJetRecurZcut_erad", &csFullJetRecurZcut_erad, &b_csFullJetRecurZcut_erad);
   fChain->SetBranchAddress("csFullJetRecurZcut_jetpt", &csFullJetRecurZcut_jetpt, &b_csFullJetRecurZcut_jetpt);
   fChain->SetBranchAddress("csFullJetRecurZcut_jeteta", &csFullJetRecurZcut_jeteta, &b_csFullJetRecurZcut_jeteta);
   fChain->SetBranchAddress("csFullJetRecurZcut_logdr12", &csFullJetRecurZcut_logdr12, &b_csFullJetRecurZcut_logdr12);
   fChain->SetBranchAddress("csFullJetRecurZcut_logztheta", &csFullJetRecurZcut_logztheta, &b_csFullJetRecurZcut_logztheta);
   fChain->SetBranchAddress("csFullJetRecurZcut_tf", &csFullJetRecurZcut_tf, &b_csFullJetRecurZcut_tf);
   fChain->SetBranchAddress("csFullJetRecurZcut_tfe", &csFullJetRecurZcut_tfe, &b_csFullJetRecurZcut_tfe);
   fChain->SetBranchAddress("csFullJetRecurZcut_tfe2", &csFullJetRecurZcut_tfe2, &b_csFullJetRecurZcut_tfe2);
   fChain->SetBranchAddress("csFullJetRecurZcut_z", &csFullJetRecurZcut_z, &b_csFullJetRecurZcut_z);
   fChain->SetBranchAddress("csFullJetRecurZcut_kt", &csFullJetRecurZcut_kt, &b_csFullJetRecurZcut_kt);
   fChain->SetBranchAddress("csFullJetRecurZcut_droppedBFS", &csFullJetRecurZcut_droppedBFS, &b_csFullJetRecurZcut_droppedBFS);
   fChain->SetBranchAddress("csFullJetRecurZcut_droppedTfBFS", &csFullJetRecurZcut_droppedTfBFS, &b_csFullJetRecurZcut_droppedTfBFS);   
   fChain->SetBranchAddress("csFullJetRecurZcut_droppedTfeBFS", &csFullJetRecurZcut_droppedTfeBFS, &b_csFullJetRecurZcut_droppedTfeBFS);
   fChain->SetBranchAddress("csFullJetRecurZcut_droppedTfeBFS", &csFullJetRecurZcut_droppedTfeBFS, &b_csFullJetRecurZcut_droppedTfeBFS);
   fChain->SetBranchAddress("csFullJetRecurZcut_droppedPts", &csFullJetRecurZcut_droppedPts, &b_csFullJetRecurZcut_droppedPts);
   fChain->SetBranchAddress("csFullJetRecurZcut_droppedKts", &csFullJetRecurZcut_droppedKts, &b_csFullJetRecurZcut_droppedKts);
   fChain->SetBranchAddress("csFullJetRecurZcut_droppedLog1drBFS", &csFullJetRecurZcut_droppedLog1drBFS, &b_csFullJetRecurZcut_droppedLog1drBFS);

   fChain->SetBranchAddress("csFullJetRecurTauZcut_jetpt", &csFullJetRecurTauZcut_jetpt, &b_csFullJetRecurTauZcut_jetpt);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_dr12", &csFullJetRecurTauZcut_dr12, &b_csFullJetRecurTauZcut_dr12);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_erad", &csFullJetRecurTauZcut_erad, &b_csFullJetRecurTauZcut_erad);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_jeteta", &csFullJetRecurTauZcut_jeteta, &b_csFullJetRecurTauZcut_jeteta);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_logdr12", &csFullJetRecurTauZcut_logdr12, &b_csFullJetRecurTauZcut_logdr12);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_logztheta", &csFullJetRecurTauZcut_logztheta, &b_csFullJetRecurTauZcut_logztheta);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_tf", &csFullJetRecurTauZcut_tf, &b_csFullJetRecurTauZcut_tf);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_tfe", &csFullJetRecurTauZcut_tfe, &b_csFullJetRecurTauZcut_tfe);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_tfe2", &csFullJetRecurTauZcut_tfe2, &b_csFullJetRecurTauZcut_tfe2);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_z", &csFullJetRecurTauZcut_z, &b_csFullJetRecurTauZcut_z);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_kt", &csFullJetRecurTauZcut_kt, &b_csFullJetRecurTauZcut_kt);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_droppedBFS", &csFullJetRecurTauZcut_droppedBFS, &b_csFullJetRecurTauZcut_droppedBFS);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_droppedTfBFS", &csFullJetRecurTauZcut_droppedTfBFS, &b_csFullJetRecurTauZcut_droppedTfBFS);   
   fChain->SetBranchAddress("csFullJetRecurTauZcut_droppedTfeBFS", &csFullJetRecurTauZcut_droppedTfeBFS, &b_csFullJetRecurTauZcut_droppedTfeBFS);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_droppedTfeBFS", &csFullJetRecurTauZcut_droppedTfeBFS, &b_csFullJetRecurTauZcut_droppedTfeBFS);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_droppedPts", &csFullJetRecurTauZcut_droppedPts, &b_csFullJetRecurTauZcut_droppedPts);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_droppedKts", &csFullJetRecurTauZcut_droppedKts, &b_csFullJetRecurTauZcut_droppedKts);
   fChain->SetBranchAddress("csFullJetRecurTauZcut_droppedLog1drBFS", &csFullJetRecurTauZcut_droppedLog1drBFS, &b_csFullJetRecurTauZcut_droppedLog1drBFS);

   fChain->SetBranchAddress("csFullJetRecur_dr12", &csFullJetRecur_dr12, &b_csFullJetRecur_dr12);
   fChain->SetBranchAddress("csFullJetRecur_erad", &csFullJetRecur_erad, &b_csFullJetRecur_erad);
   fChain->SetBranchAddress("csFullJetRecur_jetpt", &csFullJetRecur_jetpt, &b_csFullJetRecur_jetpt);
   fChain->SetBranchAddress("csFullJetRecur_jeteta", &csFullJetRecur_jeteta, &b_csFullJetRecur_jeteta);
   fChain->SetBranchAddress("csFullJetRecur_logdr12", &csFullJetRecur_logdr12, &b_csFullJetRecur_logdr12);
   fChain->SetBranchAddress("csFullJetRecur_logztheta", &csFullJetRecur_logztheta, &b_csFullJetRecur_logztheta);
   fChain->SetBranchAddress("csFullJetRecur_tf", &csFullJetRecur_tf, &b_csFullJetRecur_tf);
   fChain->SetBranchAddress("csFullJetRecur_tfe", &csFullJetRecur_tfe, &b_csFullJetRecur_tfe);
   fChain->SetBranchAddress("csFullJetRecur_tfe2", &csFullJetRecur_tfe2, &b_csFullJetRecur_tfe2);
   fChain->SetBranchAddress("csFullJetRecur_z", &csFullJetRecur_z, &b_csFullJetRecur_z);
   fChain->SetBranchAddress("csFullJetRecur_kt", &csFullJetRecur_kt, &b_csFullJetRecur_kt);

   fChain->SetBranchAddress("csFullJetRecurTau_jetpt", &csFullJetRecurTau_jetpt, &b_csFullJetRecurTau_jetpt);
   fChain->SetBranchAddress("csFullJetRecurTau_dr12", &csFullJetRecurTau_dr12, &b_csFullJetRecurTau_dr12);
   fChain->SetBranchAddress("csFullJetRecurTau_erad", &csFullJetRecurTau_erad, &b_csFullJetRecurTau_erad);
   fChain->SetBranchAddress("csFullJetRecurTau_jeteta", &csFullJetRecurTau_jeteta, &b_csFullJetRecurTau_jeteta);
   fChain->SetBranchAddress("csFullJetRecurTau_logdr12", &csFullJetRecurTau_logdr12, &b_csFullJetRecurTau_logdr12);
   fChain->SetBranchAddress("csFullJetRecurTau_logztheta", &csFullJetRecurTau_logztheta, &b_csFullJetRecurTau_logztheta);
   fChain->SetBranchAddress("csFullJetRecurTau_tf", &csFullJetRecurTau_tf, &b_csFullJetRecurTau_tf);
   fChain->SetBranchAddress("csFullJetRecurTau_tfe", &csFullJetRecurTau_tfe, &b_csFullJetRecurTau_tfe);
   fChain->SetBranchAddress("csFullJetRecurTau_tfe2", &csFullJetRecurTau_tfe2, &b_csFullJetRecurTau_tfe2);
   fChain->SetBranchAddress("csFullJetRecurTau_z", &csFullJetRecurTau_z, &b_csFullJetRecurTau_z);
   fChain->SetBranchAddress("csFullJetRecurTau_kt", &csFullJetRecurTau_kt, &b_csFullJetRecurTau_kt);

   fChain->SetBranchAddress("csFullJetRecur_partonMatchID", &csFullJetRecur_partonMatchID, &b_csFullJetRecur_partonMatchID);
   fChain->SetBranchAddress("csFullJetRecur_partonMatchDr", &csFullJetRecur_partonMatchDr, &b_csFullJetRecur_partonMatchDr);

   Notify();
}

Bool_t jetTreeSig::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void jetTreeSig::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t jetTreeSig::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef jetTreeSig_cxx
