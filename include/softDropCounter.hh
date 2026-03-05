#ifndef softDropCounter_h
#define softDropCounter_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

#include "jetCollection.hh"
#include "jewelMatcher.hh"

//---------------------------------------------------------------
// Description
// This class runs iterative SoftDrop on a set of jets
// Author: Y. Chen
//---------------------------------------------------------------

class softDropCounter
{
private :
  int fRecursiveAlgo_; //clustering algo to be used for the declustering
  double zcut_;        // ZCut used in SD
  double beta_;        // Beta parameter
  double r0_;          // Jet radius
  double rcut_;        // Termination Criteria

  bool doJewelSub_;                                 // do JEWEL 4MomSub
  std::vector<fastjet::PseudoJet> particlesDummy_;  // dummy particles for 4MomSub

  std::vector<fastjet::PseudoJet>  fjInputs_;    // ungroomed jets
  std::vector<std::vector<double>> droppedBFS_;  // amount of branches dropped before first split
  std::vector<std::vector<double>> droppedTfBFS_;  // amount of branches dropped before first split
  std::vector<std::vector<double>> droppedTfeBFS_;  // amount of branches dropped before first split
  std::vector<std::vector<double>> droppedPts_;
  std::vector<std::vector<double>> droppedKts_;
  std::vector<std::vector<double>> droppedLog1drs_;
  std::vector<std::vector<double>> zgs_;         // all the zg's in the algorithm
  std::vector<std::vector<double>> drs_;         // and the angles in the algorithm
  std::vector<std::vector<double>> pts_;         // pts of each node
  std::vector<std::vector<double>> etas_;        // etas of each node
  std::vector<std::vector<double>> erads_;       // energy of sum of the two branches at each node
  std::vector<std::vector<double>> log1drs_;     // Log(1/angle) in the algorithm
  std::vector<std::vector<double>> logzdrs_;     // Log(z*angle) in the algorithm
  std::vector<std::vector<double>> tfs_;         // formation time tf = soft&collinear
  std::vector<std::vector<double>> tfes_;        // formation time tf = 1 / 2*(z*(1-z)*E*(1-cos(theta)))
  std::vector<std::vector<double>> tfes2_;       // formation time tf = for 2 splits
  std::vector<std::vector<double>> kts_;         // all the kt's in the algorithm
  std::vector<double> tau21s_;      // nsubjettiness
  std::vector<double> tau32s_;      // nsubjettiness
  
public :
  softDropCounter(double z = 0.1, double beta = 0.0, double r0 = 0.4, double rcut = 0.1);
  void setZCut(double c);
  void setBeta(double b);
  void setR0(double r) ;
  void setRCut(double r);
  void setInputJets(const std::vector<fastjet::PseudoJet> &v);

  void setRecursiveAlgo(int r);

  void setJewelSubFlag(bool b);
  void doJewelSub(std::vector<fastjet::PseudoJet> dummies);

  std::vector<std::vector<double>> getDBFSs() const { return droppedBFS_; }
  std::vector<std::vector<double>> getDTfBFSs() const { return droppedTfBFS_; }
  std::vector<std::vector<double>> getDTfeBFSs() const { return droppedTfeBFS_; }
  std::vector<std::vector<double>> getPtBFSs() const { return droppedPts_; }
  std::vector<std::vector<double>> getKtBFSs() const { return droppedKts_; }
  std::vector<std::vector<double>> getLog1DrBFSs() const { return droppedLog1drs_; }
  std::vector<std::vector<double>> getZgs() const { return zgs_; }
  std::vector<std::vector<double>> getDRs() const { return drs_; }
  std::vector<std::vector<double>> getPts() const { return pts_; }
  std::vector<std::vector<double>> getEtas() const { return etas_; }
  std::vector<std::vector<double>> getErads() const { return erads_; }
  std::vector<std::vector<double>> getLog1DRs() const { return log1drs_; }
  std::vector<std::vector<double>> getLogzDRs() const { return logzdrs_; }
  std::vector<std::vector<double>> getTfs() const { return tfs_; }
  std::vector<std::vector<double>> getTfes() const { return tfes_;}
  std::vector<std::vector<double>> getTfes2() const { return tfes2_;}
  std::vector<std::vector<double>> getKts() const { return kts_; }
  std::vector<double> getTau21s() const { return tau21s_; }
  std::vector<double> getTau32s() const { return tau32s_; }

  void run(const jetCollection &c);
  void run(const std::vector<fastjet::PseudoJet> &v);
  void run();
  std::vector<double> calculateNSD(double Kappa, double AngleKappa = 0);
};

softDropCounter::softDropCounter(double z, double beta, double r0, double rcut)
   : zcut_(z), beta_(beta), r0_(r0), rcut_(rcut)
{
  fRecursiveAlgo_ = 0;
  doJewelSub_     = false;
}

void softDropCounter::setZCut(double c)
{
   zcut_ = c;
}

void softDropCounter::setBeta(double b)
{
   beta_ = b;
}

void softDropCounter::setR0(double r)
{
   r0_   = r;
}

void softDropCounter::setRCut(double r)
{
   rcut_ = r;
}

void softDropCounter::setInputJets(const std::vector<fastjet::PseudoJet> &v)
{
   fjInputs_ = v;
}

void softDropCounter::setRecursiveAlgo(int r)
{
  fRecursiveAlgo_ = r;
}

void softDropCounter::setJewelSubFlag(bool b)
{
  doJewelSub_ = b;
}

void softDropCounter::doJewelSub(std::vector<fastjet::PseudoJet> dummies) {
  setJewelSubFlag(true);
  particlesDummy_ = dummies;
}


void softDropCounter::run(const jetCollection &c)
{
   run(c.getJet());
}

void softDropCounter::run(const std::vector<fastjet::PseudoJet> &v)
{
   setInputJets(v);
   run();
}

void softDropCounter::run()
{
  for(fastjet::PseudoJet &jet: fjInputs_)
   {
      if(jet.has_constituents() == false)
      {
	 droppedBFS_.push_back(vector<double>());
	 droppedTfBFS_.push_back(vector<double>());
	 droppedTfeBFS_.push_back(vector<double>());
         droppedPts_.push_back(vector<double>());
         droppedKts_.push_back(vector<double>());
         droppedLog1drs_.push_back(vector<double>());
	 zgs_.push_back(vector<double>());
         drs_.push_back(vector<double>());
         pts_.push_back(vector<double>());
         etas_.push_back(vector<double>());
         erads_.push_back(vector<double>());
         log1drs_.push_back(vector<double>());
         logzdrs_.push_back(vector<double>());
         tfs_.push_back(vector<double>());
         tfes_.push_back(vector<double>());
         tfes2_.push_back(vector<double>());
         kts_.push_back(vector<double>());
         tau21s_.push_back(-1.);
         tau32s_.push_back(-1.);
         continue;
      }

      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      //need flexibility here on how we order the tree before declustering
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

      fastjet::JetDefinition jet_def(jetalgo, fastjet::JetDefinition::max_allowable_R, kpval, static_cast<fastjet::RecombinationScheme>(0));
      fastjet::ClusterSequence cs(particles, jet_def);
      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());

      double hbarc = 0.19732697;
      double GeVtofm = 1./hbarc; //~5.068;

      if(tempJets.size() == 0)
      {
	 droppedBFS_.push_back(vector<double>());
	 droppedTfBFS_.push_back(vector<double>());
	 droppedTfeBFS_.push_back(vector<double>());
	 droppedPts_.push_back(vector<double>());
	 droppedKts_.push_back(vector<double>());
         droppedLog1drs_.push_back(vector<double>());
         zgs_.push_back(vector<double>());
         drs_.push_back(vector<double>());
         pts_.push_back(vector<double>());
         etas_.push_back(vector<double>());
         erads_.push_back(vector<double>());
         log1drs_.push_back(vector<double>());
         logzdrs_.push_back(vector<double>());
         tfs_.push_back(vector<double>());
         tfes_.push_back(vector<double>());
         tfes2_.push_back(vector<double>());
         kts_.push_back(vector<double>());
         tau21s_.push_back(-1.);
	 tau32s_.push_back(-1.);
	 continue;
      }

      fastjet::PseudoJet CurrentJet = tempJets[0];
      fastjet::PseudoJet Part1, Part2;

      std::vector<double> dropped;
      std::vector<double> droppedTf;
      std::vector<double> droppedTfe;
      std::vector<double> droppedPt;
      std::vector<double> droppedKt;
      std::vector<double> droppedLog1dr;

      std::vector<double> z;
      std::vector<double> dr;
      std::vector<double> pt;
      std::vector<double> eta;
      std::vector<double> erad;
      std::vector<double> log1dr;
      std::vector<double> logzdr;
      std::vector<double> tf;
      std::vector<double> tfe;
      std::vector<double> tfe2;
      std::vector<double> kt;
      bool beforeFirstSplit = true;
      bool secondSplit = false;
      bool thirdSplit = false;
      bool fullJet = true;
      double z03 = 0;
      double ptJet = 0;
      fastjet::PseudoJet sj3;//for 1->3 split
      fastjet::PseudoJet FullJet;
      while(CurrentJet.has_parents(Part1, Part2))
      {
         if(CurrentJet.pt2() <= 0)
            break;

         fastjet::PseudoJet sj1;
         fastjet::PseudoJet sj2;
         if(doJewelSub_) {
           sj1 = GetCorrectedJet(Part1,particlesDummy_);
           sj2 = GetCorrectedJet(Part2,particlesDummy_);
         } else  {
           sj1 = Part1;
           sj2 = Part2;
         }
      
         double DeltaR = std::sqrt(sj1.squared_distance(sj2));
         if(DeltaR < rcut_)
            break;

         double PT1 = sj1.pt();
         double PT2 = sj2.pt();
         double zg = -1;
         if(PT1 + PT2 > 0)
           zg = min(PT1, PT2) / (PT1 + PT2);
         else
           break;

         //double z1 = max(sj1.e(),sj2.e())/CurrentJet.e();
         //double z2 = min(sj1.e(),sj2.e())/CurrentJet.e();
         
         if(fullJet == true){
           FullJet = CurrentJet;
	   fullJet = false;
	 }
	 
	 double Threshold = zcut_ * std::pow(DeltaR / r0_, beta_);
         if(zg >= Threshold)   // yay
	 {
           if(fullJet == true){
             FullJet = CurrentJet;
             fullJet = false;
            }
            beforeFirstSplit = false;
            z.push_back(zg);
            dr.push_back(DeltaR);
            pt.push_back(CurrentJet.perp());
            eta.push_back(CurrentJet.eta());
            erad.push_back(sj1.e()+sj2.e());
            log1dr.push_back(log(1./DeltaR));
            logzdr.push_back(log(zg*DeltaR));
            //tf.push_back(2./(zg*(1.-zg)*CurrentJet.perp()*GeVtofm*DeltaR*DeltaR/r0_/r0_));//original
            //tf.push_back(2./(zg*CurrentJet.perp()*GeVtofm*(DeltaR/r0_)*(DeltaR/r0_)));//Johanna soft collinear limit of the tfe (z*E ~ omega, (1-z) ~ 1) cos(theta_12) ~ 1- theta^{2}/2
            //tfe.push_back(1./(2.*z1*z2*CurrentJet.e()*GeVtofm*(1-fastjet::cos_theta(sj1,sj2))));//original
            tf.push_back(1./(2.*zg*(1-zg)*CurrentJet.perp()*GeVtofm*(1-cos(DeltaR/r0_))));//Johanna
            tfe.push_back(hbarc*(sj1.perp() + sj2.perp())/( sj1.perp()*sj2.perp() * DeltaR * DeltaR));//Morgan
            kt.push_back(min(PT1,PT2)*DeltaR);
            
	    if(thirdSplit == true){
	      //only true, once we entered the 2nd split/filled values for the 1->3 formation time
	    }else{
	      if(secondSplit == true){
		double DeltaR13 = 0;
		double DeltaR23 = 0;
		if(PT2<PT1){
		  DeltaR13 = std::sqrt(sj3.squared_distance(sj1));
		  DeltaR23 = std::sqrt(sj3.squared_distance(sj2));
		}
		if(PT1<PT2){
                  DeltaR13 = std::sqrt(sj3.squared_distance(sj2));
                  DeltaR23 = std::sqrt(sj3.squared_distance(sj1));
		}
		double term12 = zg*(1-zg)*(1-z03)*ptJet*(1-cos(DeltaR/r0_));
		double term13 = z03*(1-zg)*(1-z03)*ptJet*(1-cos(DeltaR13/r0_));
		double term23 = z03*zg*(1-z03)*ptJet*(1-cos(DeltaR23/r0_));
		tfe2.push_back(1/(2*GeVtofm*(term12+term13+term23)));
		//std::cout<<"zg: "<<zg<<" z03: "<<z03<<" ptJet: "<<ptJet<<" cos(DeltaR/r0_): "<<cos(DeltaR/r0_)<<" cos(DeltaR13/r0_): "<<cos(DeltaR13/r0_)<<" cos(DeltaR23/r0_): "<<cos(DeltaR23/r0_)<<std::endl;
		//std::cout<<"Term1: "<<term12<<" Term2: "<<term13<< " Term3: "<<term23<<" 1->3 tf: "<<1/(2*GeVtofm*(term12+term13+term23))<<std::endl;
		thirdSplit = true;// to avaiod a nother calculation
	      }
	      z03 = zg;//store value for second split calculation in next iteration
              ptJet = CurrentJet.perp();
	    }
	    secondSplit = true; 
	 }
	 if((beforeFirstSplit == true)){
	    dropped.push_back(DeltaR);
	    droppedTf.push_back(2./(zg*CurrentJet.perp()*GeVtofm*(DeltaR/r0_)*(DeltaR/r0_)));
	    droppedTfe.push_back(1./(2.*zg*(1-zg)*CurrentJet.perp()*GeVtofm*(1-cos(DeltaR/r0_))));
            droppedPt.push_back(CurrentJet.perp());
	    droppedKt.push_back(min(PT1, PT2)*DeltaR);
	    droppedLog1dr.push_back(log(1/DeltaR));
	 }

         if(PT1 > PT2){
            CurrentJet = Part1;
	    sj3 = Part2;
	 }
         if(PT2 > PT1){
            CurrentJet = Part2;
	    sj3 = Part1;
	 }
      }//end while loop
      droppedBFS_.push_back(dropped);
      droppedTfBFS_.push_back(droppedTf);
      droppedTfeBFS_.push_back(droppedTfe);
      droppedLog1drs_.push_back(droppedLog1dr);
      droppedPts_.push_back(droppedPt);
      droppedKts_.push_back(droppedKt);
      zgs_.push_back(z);
      drs_.push_back(dr);
      pts_.push_back(pt);
      etas_.push_back(eta);
      erads_.push_back(erad);
      log1drs_.push_back(log1dr);
      logzdrs_.push_back(logzdr);
      tfs_.push_back(tf);
      tfes_.push_back(tfe);
      tfes2_.push_back(tfe2);
      kts_.push_back(kt);

      //Compute the n-subjettiness ratio for C/A indexed jet
      // N-subjettiness with Unnormalized Measure (in GeV)
      // beta = 1.0:  One-Winner-Take-All kT Axes
      // beta = 2.0:  One-pass E-Scheme kT Axes
      fastjet::PseudoJet FullPart1, FullPart2;
      double tau21_beta2;
      double tau32_beta2;
      if(FullJet.has_structure()){
        double beta = 2;                                                                                 
        fastjet::contrib::NsubjettinessRatio nSub21_beta2(2,1, fastjet::contrib::CA_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));//use CA and kt  --check defualt vs Unnorm.
        fastjet::contrib::NsubjettinessRatio nSub32_beta2(3,2, fastjet::contrib::KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
  
        tau21_beta2 = nSub21_beta2(FullJet);                                                  
        tau32_beta2 = nSub32_beta2(FullJet);                                                  
      }else{
        tau21_beta2 = -99;
	tau32_beta2 = -99;
    
      }
      tau21s_.push_back(tau21_beta2);                                                                         
      tau32s_.push_back(tau32_beta2); 
   }
}

std::vector<double> softDropCounter::calculateNSD(double Kappa, double AngleKappa)
{
   std::vector<double> Result;

   for(int i = 0; i < (int)zgs_.size(); i++)
   {
      double Total = 0;
      for(int j = 0; j < (int)zgs_[i].size(); j++)
         Total = Total + pow(zgs_[i][j], Kappa) * pow(drs_[i][j], AngleKappa);
      Result.push_back(Total);
   }

   return Result;
}

#endif
