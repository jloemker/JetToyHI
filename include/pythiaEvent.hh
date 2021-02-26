#ifndef pythiaEvent_h
#define pythiaEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "Pythia8/Pythia.h"

#include "extraInfo.hh"

//using namespace std;

//---------------------------------------------------------------
// Description
// This class generates a pythia8 event
// Author: M. Verweij
//---------------------------------------------------------------

class pythiaEvent {

private :
  Pythia8::Pythia pythia;
  double pthat_;
  unsigned int tune_;
  double rapMin_;
  double rapMax_;
  bool   partonLevel_;
  bool   vinciaShower_;

  std::vector<fastjet::PseudoJet> partons;

public :
  pythiaEvent(double pthat = 120., unsigned int tune = 14, double rapMin = -3., double rapMax = 3., bool partonLevel = false, bool vinciaShower = false);
  std::vector<fastjet::PseudoJet> createPythiaEvent();
  
  std::vector<fastjet::PseudoJet> getPartonList() const { return partons; }

  void getStat() {pythia.stat();}

};
  
pythiaEvent::pythiaEvent(double pthat, unsigned int tune, double rapMin, double rapMax, bool partonLevel, bool vinciaShower) :
  pthat_(pthat), tune_(tune), rapMin_(rapMin), rapMax_(rapMax), partonLevel_(partonLevel), vinciaShower_(vinciaShower)
{
    
  // Generator. LHC process and output selection. Initialization.
  // tunes: http://home.thep.lu.se/~torbjorn/pythia82html/Tunes.html
  pythia.readString("Beams:eCM = 5002.");
  pythia.readString("HardQCD:all = on");
  pythia.readString(Form("PhaseSpace:pTHatMin = %.1f",pthat_));
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  if(partonLevel_) {
    pythia.readString("HadronLevel:all = off");
  }
  if(vinciaShower_)
    pythia.readString("PartonShowers:Model = 2"); //activate the VINCIA parton shower
  else
      pythia.readString(Form("Tune:pp = %d",tune_));
  
  pythia.init();

}

std::vector<fastjet::PseudoJet> pythiaEvent::createPythiaEvent() {

  pythia.next(); //generate next event
    
  std::vector<fastjet::PseudoJet> particles;
  partons.clear(); //empty list before storing partons of new event

  for (int i = 0; i < pythia.event.size(); ++i) {
    if (pythia.event[i].isFinal()) { //all final state particles
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), 0)); 
      if(p.rap()>rapMin_ && p.rap()<rapMax_)
        particles.push_back(p);
    } else if (pythia.event[i].status()==-23) { //outgoing partons from the hard scattering
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), -1)); 
      partons.push_back(p);
      
      //find the case where the splitting to two separate daughters happens
      int d1 = pythia.event[i].daughter1();
      int d2 = pythia.event[i].daughter2();
      while(d1==d2) {
        d1 = pythia.event[d1].daughter1();
        d2 = pythia.event[d2].daughter2();
      }

      fastjet::PseudoJet pd1(pythia.event[d1].px(),pythia.event[d1].py(),pythia.event[d1].pz(),pythia.event[d1].e());
      pd1.set_user_info(new extraInfo(pythia.event[d1].id(), -2)); 
      partons.push_back(pd1);

      fastjet::PseudoJet pd2(pythia.event[d2].px(),pythia.event[d2].py(),pythia.event[d2].pz(),pythia.event[d2].e());
      pd2.set_user_info(new extraInfo(pythia.event[d2].id(), -2)); 
      partons.push_back(pd2);
      
    }
  }

  // pythia.event.list();

  return particles;
}


#endif
