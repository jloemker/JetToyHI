#include "jetTreeSig.C"
#include "partonTreeSig.C"

void runJetTreeSig() {
  // gROOT->LoadMacro("");
  jetTreeSig oldtree;
  //oldtree.Loop(false, "PYTHIA", "short", "JetToyHIResultTimeClusBkg_13p6_200",0.4);
  oldtree.Loop(false, "JEWEL", "short", "JetToyHIResultTimeClusBkgJewel",0.4);










  /*JetToyHIResultTimeClusBkg_R04_Kt5_100000ev
  TString PYTHIA04Kt5 = "Johanna/Pythia/JetToyHIResultTimeClusBkg_R04_CsPtE09Kt5_100000ev.root";
  TString PYTHIA04NoKt = "Johanna/Pythia/JetToyHIResultTimeClusBkg_R04_CsPtE09NoKt_100000ev.root";
  TString PYTHIA04LP2 = "Johanna/Pythia/JetToyHIResultTimeClusBkg_R04_CsPtE09LP2_100000ev.root";
  TString PYTHIA04LP4 = "Johanna/Pythia/JetToyHIResultTimeClusBkg_R04_CsPtE09LP4_100000ev.root";

        if(reduced==true && generator=="CA"){
          if(!(sigTfsZ.size()>0)) continue;
          if(duration == "supershort" && (sigTfsZ[0] >= 0.5)) continue
          if(duration == "short" && (sigTfsZ[0] > 1)) continue;
          if(duration == "large" && (sigTfsZ[0] < 3)) continue;
          if(duration == "intermed" && ((sigTfsZ[0] < 1) || (sigTfsZ[0] > 3))) continue;
        }
        if(reduced==true && generator=="Tau"){
          if(!(sigTfsZTau.size()>0)) continue;
          if(duration == "supershort" && (sigTfsZ[0] >= 0.5)) continue
          if(duration == "short" && (sigTfsZTau[0] > 1)) continue;
          if(duration == "large" && (sigTfsZTau[0] < 3)) continue;
          if(duration == "intermed" && ((sigTfsZTau[0] < 1) || (sigTfsZTau[0] > 3))) continue;
        
  */
}
