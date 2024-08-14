#ifndef dyGroomer_h
#define dyGroomer_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/WrappedStructure.hh"
#include "fastjet/contrib/RecursiveSymmetryCutBase.hh"

#include "jetCollection.hh"
#include "jetMatcher.hh"
//---------------------------------------------------------------
// Description
// This class runs Dynamical grooming on a set of jets
// We include the possibility to have impose an SD condition
// Author: Alba Soto-Ontoso
//---------------------------------------------------------------

class dyGroomer {

private :
        int fRecursiveAlgo_;
	double a_;

  std::vector<fastjet::PseudoJet> fjInputs_;   //ungroomed jets
  std::vector<fastjet::PseudoJet> fjOutputs_;  //groomed jets
  std::vector<fastjet::PseudoJet> fjDaughters1_;  //daughter 1
  std::vector<fastjet::PseudoJet> fjDaughters2_;  //daughter 2
  std::vector<std::vector<double>> zg_;         //zg of groomed jets
  std::vector<std::vector<double>>                drBranches_; //dropped branches
  std::vector<std::vector<double>>             dr12_;       //distance between the two subjet
  std::vector<std::vector<double>>             kt_;         //kt of groomed jets
  std::vector<std::vector<double>>             tau21_;      //n-subjettiness ratio 21
  std::vector<std::vector<double>>             tau32_;      //n-subjettiness ratio 32
  //added for sdc  C/A comparison
  std::vector<std::vector<double>>		  tfs_;     	// formation time tf = 2 / (z*(1-z)*pt*theta^2) (soft limit)
  std::vector<std::vector<double>>		  tfes_;     	// formation time tf = 1 / 2*(z*(1-z)*E*(1-cos(theta))) (hard limit)

  std::vector<std::vector<double>>             kappa_;      //kappa of groomed jets

public :
  dyGroomer(double a = 1);
  void setInputJets(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> getGroomedJets() const;
  std::vector<fastjet::PseudoJet> getDaughters1() const;
  std::vector<fastjet::PseudoJet> getDaughters2() const;

  void setRecursiveAlgo(int r);
  std::vector<std::vector<double>> getZgs() const { return zg_; }
  std::vector<double>  calculateZg();
  std::vector<std::vector<double>> getDR12() const { return dr12_; }
  std::vector<std::vector<double>> getKts() const { return kt_; }
  std::vector<std::vector<double>>    getNDroppedSubjets() const { return drBranches_; }
  std::vector<std::vector<double>> getKappas() const { return kappa_; }
  std::vector<std::vector<double>> getTau21() const { return tau21_; }
  std::vector<std::vector<double>> getTau32() const { return tau32_; }
  //formation time information
  std::vector<std::vector<double>> getTfs() const { return tfs_; }
  std::vector<std::vector<double>> getTfes() const { return tfes_; }
  double getKappa(double pt, double theta, double z);
  std::vector<fastjet::PseudoJet> doGrooming(jetCollection &c);
  std::vector<fastjet::PseudoJet> doGrooming(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> doGrooming();
};

dyGroomer::dyGroomer(double a)
   : a_(a)
{
  fRecursiveAlgo_=0;
}

void dyGroomer::setInputJets(const std::vector<fastjet::PseudoJet> v)
{
   fjInputs_ = v;
}

void dyGroomer::setRecursiveAlgo(int r)
{
  fRecursiveAlgo_ = r;
}

std::vector<fastjet::PseudoJet> dyGroomer::getGroomedJets() const
{
   return fjOutputs_;
}

std::vector<fastjet::PseudoJet> dyGroomer::getDaughters1() const
{
   return fjDaughters1_;
}

std::vector<fastjet::PseudoJet> dyGroomer::getDaughters2() const
{
   return fjDaughters2_;
}
/*
std::vector<double> dyGroomer::getZgs() const
{
   return zg_;
}

std::vector<double> dyGroomer::getDR12() const
{
   return dr12_;
}

std::vector<double> dyGroomer::getKts() const
{
   return kt_;
}

std::vector<double> dyGroomer::getTau21() const
{
   return tau21_;
}

std::vector<double> dyGroomer::getTau32() const
{
   return tau32_;
}

std::vector<double> dyGroomer::getTfs() const
{
   return tfs_;
}

std::vector<double> dyGroomer::getTfes() const
{
   return tfes_;
}


std::vector<int> dyGroomer::getNDroppedSubjets() const
{
   return drBranches_;
}


std::vector<double> dyGroomer::getKappas() const
{
   return kappa_;
}
*/
std::vector<fastjet::PseudoJet> dyGroomer::doGrooming(jetCollection &c)
{
   return doGrooming(c.getJet());
}

std::vector<fastjet::PseudoJet> dyGroomer::doGrooming(std::vector<fastjet::PseudoJet> v)
{
   setInputJets(v);
   return doGrooming();
}

double dyGroomer::getKappa(double pt, double theta, double z)
{  
   double kappa = 1/(z * (1-z) * pt * pow(theta,a_));
  // double kappa = 1/(z * pt * pow(theta,a_));
  //ma double kappa = 1/(z * pt * pow(theta,a_));
   return kappa;
}

std::vector<fastjet::PseudoJet> dyGroomer::doGrooming()
{
   fjOutputs_.reserve(fjInputs_.size());
   fjDaughters1_.reserve(fjInputs_.size());
   fjDaughters2_.reserve(fjInputs_.size());
   zg_.reserve(fjInputs_.size());
   dr12_.reserve(fjInputs_.size());
   kt_.reserve(fjInputs_.size());
   drBranches_.reserve(fjInputs_.size());
   kappa_.reserve(fjInputs_.size());
   tau21_.reserve(fjInputs_.size());
   tau32_.reserve(fjInputs_.size());
   tfs_.reserve(fjInputs_.size());
   tfes_.reserve(fjInputs_.size());

   int ijet = -1;

   for(fastjet::PseudoJet& jet : fjInputs_) {
      ++ijet;

      if(!jet.has_constituents()) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         fjDaughters1_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         fjDaughters2_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(vector<double>());
         dr12_.push_back(vector<double>());
         kt_.push_back(vector<double>());
         drBranches_.push_back(vector<double>());
         kappa_.push_back(vector<double>());
         tau21_.push_back(vector<double>());
         tau32_.push_back(vector<double>());
	 tfs_.push_back(vector<double>());
         tfes_.push_back(vector<double>());
         continue;
      }
      
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
      // adjusted accodring to softDropCounter.hh
      // flexibility for tree ordering before grooming
      fastjet::JetAlgorithm jetalgo(fastjet::genkt_algorithm);
      float kpval=0;
      if(fRecursiveAlgo_==0){
        //CA
	kpval=0;
      }
      if(fRecursiveAlgo_==1){
        //kT
	kpval=-1;
      }
      if(fRecursiveAlgo_==2){
        //antikt
	kpval=1;
      }
      if(fRecursiveAlgo_==3){
        //tform-ordered
	kpval=0.5;
      }
      //std::cout<<"Recursive algo "<<fRecursiveAlgo_<<" kpval "<<kpval<<std::endl;
      fastjet::JetDefinition jet_def(jetalgo, fastjet::JetDefinition::max_allowable_R, kpval, static_cast<fastjet::RecombinationScheme>(0));
      fastjet::ClusterSequence cs(particles, jet_def);
      vector<fastjet::ClusterSequence::history_element> cs_history = cs.history();
      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());
      //std::cout<<"cs_history "<<cs_history<<std::endl;
      //std::cout<<"tempJets "<<tempJets<<std::endl;
      // To make the jetp_index work we need to define our jets like this
      std::vector<fastjet::PseudoJet> tempJets_two = cs.jets();

      if(tempJets.size()<1) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         fjDaughters1_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         fjDaughters2_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(vector<double>());
         dr12_.push_back(vector<double>());
         kt_.push_back(vector<double>());
         drBranches_.push_back(vector<double>());
         kappa_.push_back(vector<double>());
         tau21_.push_back(vector<double>());
         tau32_.push_back(vector<double>());
         tfs_.push_back(vector<double>());
         tfes_.push_back(vector<double>());
         continue;
      }

      fastjet::PseudoJet CurrentJet = tempJets[0];
      //std::cout<<"======================CurrentJet "<< CurrentJet <<std::endl;
      fastjet::PseudoJet piece1, piece2;
      double min_kappa = 1e8;
      double min_tf = 1e8;
      //double zg = -1.;
      //double deltaR = -1;
      //double pt = -1;
      double ndrop = 0.0;
      int current_in_ca_tree = -1; // (history) index of the current particle in the C/A or tau tree 


      std::vector<double> z;
      std::vector<double> dr12;
      std::vector<double> kt;
      std::vector<double> drBranches;
      std::vector<double> kappa;
      std::vector<double> tau21;
      std::vector<double> tau32;
      std::vector<double> tfs;
      std::vector<double> tfes;
      // now recurse into the jet's structure to find the maximum hardness/ minimum inverse hardness
  while (CurrentJet.has_parents(piece1, piece2)) {

     if (CurrentJet.pt2() <= 0) break;

    // if(piece1.pt() + piece2.pt() > 0 && piece1.E()>0. && piece2.E()>0. && piece1.m()>0. && piece2.m()>0.){
     if(piece1.pt() + piece2.pt() > 0 && piece1.E()>0. && piece2.E()>0.){
     double pt = piece1.pt() + piece2.pt();
     double zg = min(piece1.pt(), piece2.pt()) / pt;
     double deltaR = piece1.delta_R(piece2);
     double kap = getKappa(pt,deltaR,zg);

     double hbarc = 0.19732697;
     double GeVtofm = 1./hbarc; //~5.068;
     //for formation time veto
     double z1 = max(piece1.e(),piece2.e())/CurrentJet.e();
     double z2 = min(piece1.e(),piece2.e())/CurrentJet.e();
     double tf = (2./(zg*(1.-zg)*CurrentJet.perp()*GeVtofm*deltaR*deltaR));// actually this is missing: /r0_/r0_)); but r/r = 1 ...
     double tfe = (1./(2.*z1*z2*CurrentJet.e()*GeVtofm*(1-fastjet::cos_theta(piece1,piece2))));//sj1,sj2

     if(tf < min_tf) {// shouldn't this criterium be pt based (eg soft and hard limit transition)
        min_tf = tf; //tfe gives different values for the printout here and the one below - tf is consistent and also reduces pt while decreasing tf !
	//std::cout<<"kappa "<<kap<<std::endl;
        current_in_ca_tree = CurrentJet.cluster_hist_index();
	//std::cout<<"formatio time: "<<tf<<std::endl;
	//std::cout<<"zg: "<<zg<<std::endl;
	//std::cout<<"pt: "<<pt<<std::endl;
	//std::cout<<"current_in_ca_tree "<<current_in_ca_tree<<std::endl;
      }
      else{
        ndrop++;
       // std::cout<<"drop that.."<<std::endl;
      }
    }

    if(piece1.pt() > piece2.pt())
      CurrentJet = piece1;
    else
      CurrentJet = piece2;
    //std::cout<<"++ ++ ++ ++ ++ ++ ++ ++ ++ Make a piece the new Current jet"<<std::endl;

  }//extend the loop() ??
    //std::cout<<"no more pieces in this jet"<<std::endl;
    if(current_in_ca_tree >= 0){// generalized cluster tagg (not only ca anymore)
    fastjet::PseudoJet groomed_jet = tempJets_two[cs_history[current_in_ca_tree].jetp_index];//but here we define the tree based on the hardest split from above in the ca history
    //std::cout<<"current_in_ca_tree "<<current_in_ca_tree<<std::endl;
    // Obtain the (z,pT,theta) of the selected splitting
     fastjet::PseudoJet daughter1, daughter2;
     groomed_jet.has_parents(daughter1, daughter2);

     //  if(daughter1.pt() + daughter2.pt() > 0 && daughter1.E()>0. && daughter2.E()>0. && daughter1.m()>0. && daughter2.m()>0.){
     if (daughter1.pt() + daughter2.pt() > 0 && daughter1.E()>0. && daughter2.E()>0){
       double pt = daughter1.pt() + daughter2.pt();
       double zg = min(daughter1.pt(), daughter2.pt()) / pt;
       double deltaR = daughter1.delta_R(daughter2);
       double kT = min(daughter1.pt(), daughter2.pt()) * deltaR;
       double hbarc = 0.19732697;
       double GeVtofm = 1./hbarc; //~5.068;
       //formation time
       double z1 = max(daughter1.e(),daughter2.e())/groomed_jet.e();
       double z2 = min(daughter1.e(),daughter2.e())/groomed_jet.e();
       double tf = (2./(zg*(1.-zg)*groomed_jet.perp()*GeVtofm*deltaR*deltaR));// actually this is missing: /r0_/r0_)); but r/r = 1 ...
       double tfe = (1./(2.*z1*z2*groomed_jet.e()*GeVtofm*(1-fastjet::cos_theta(daughter1,daughter2))));
       //std::cout<<"formation time: "<<tf<<std::endl;
      // std::cout<<"groomed_jet.constituents() "<<groomed_jet.constituents()<<std::endl;
       drBranches.push_back(ndrop);
       z.push_back(zg);
       dr12.push_back(deltaR);
       kt.push_back(kT);
       //kappa_.push_back(kappa);
       kappa.push_back(1./getKappa(pt,deltaR,zg));
       //formation time for C/A indexed split
       tfs.push_back(tf);// calculated above
       tfes.push_back(tfe);// calculated above

     }

     // Compute the n-subjettiness ratio for C/A indexed jet
     // N-subjettiness with Unnormalized Measure (in GeV)
     // beta = 1.0:  One-Winner-Take-All kT Axes
     //beta = 2.0:  One-pass E-Scheme kT Axes
     //std::cout<<"groomed_jet "<<groomed_jet<<std::endl;

     double beta = 2;
     fastjet::contrib::NsubjettinessRatio nSub21_beta2(2,1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
     fastjet::contrib::NsubjettinessRatio nSub32_beta2(3,2, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));

     double tau21_beta2 = nSub21_beta2(groomed_jet);
     double tau32_beta2 = nSub32_beta2(groomed_jet);

     tau21.push_back(tau21_beta2);
     tau32.push_back(tau32_beta2);
     
     zg_.push_back(z);
     dr12_.push_back(dr12);
     kt_.push_back(kt);
     drBranches_.push_back(drBranches);
     kappa_.push_back(kappa);
     tfs_.push_back(tfs);
     tfes_.push_back(tfes);
    
     tau21_.push_back(tau21);
     tau32_.push_back(tau32);

    fjOutputs_.push_back(groomed_jet); //put CA reclusterd jet after grooming
    fjDaughters1_.push_back(daughter1);
    fjDaughters2_.push_back(daughter2);
   }

   else {fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(vector<double>());
         dr12_.push_back(vector<double>());
         kt_.push_back(vector<double>());
         drBranches_.push_back(vector<double>());
         kappa_.push_back(vector<double>());
         tau21_.push_back(vector<double>());
         tau32_.push_back(vector<double>());
	 tfs_.push_back(vector<double>());
	 tfes_.push_back(vector<double>());
   }
 }
   return fjOutputs_;
}


#endif
