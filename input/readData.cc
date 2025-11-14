#include "TFile.h"
#include "TTree.h"
// add histos per run number 
void read(TString inputFile){
    TFile *fin = new TFile(inputFile);
    fin->ls();
    // define all the column names u want to read in 
    char tracks[10] = "/O2jtrack"; 
    const char *track_name = tracks;
    
    char collision[13] = "/O2collision"; 
    const char *collision_name = collision;
    
    // loop over the directories (DF's) of O2
    for(auto k : *fin->GetListOfKeys()){
        const char *dir_name = k->GetName();// extract name of the current DF
        char buff_trk[100];// add name of DF and track column together
        strncpy(buff_trk, dir_name, sizeof(buff_trk));
        strncat(buff_trk, track_name, sizeof(buff_trk)-10);
        const char *tree_trk = buff_trk;
        TTree* trk = fin->Get<TTree>(tree_trk);// get the tree from the corresponding DF
        if (trk == NULL) continue;// skip the parenFiles
        
        //trk->ls();
        //trk->Print();   
        
        Double_t pt;// declare all variable you want to extract
        Double_t eta;
        Double_t phi;
        Int_t idxColl = 0;
        Int_t idx = 0;// to frind out if we are at the next collision
        trk->SetBranchAddress("fPt",&pt);
        trk->SetBranchAddress("fEta",&eta);
        trk->SetBranchAddress("fPhi",&phi);
        trk->SetBranchAddress("fIndexJCollisions",&idxColl);
        for(int i=0; i<trk->GetEntries(); i++){// loop over the tree (ordered by collisions)
          trk->GetEntry(i);
          if(idxColl > idx ){
            cout<<"================================================================"<<endl;
            idx +=1;
          } 
          cout<<"idxColl "<<idxColl<<endl;
          cout<<"pt: "<<pt;
          cout<<" eta: "<<eta;
          cout<<" phi: "<<phi<<endl;
        }
        
    }
}

void readAndWrite(TString inputFile, TString outputFile){
    TH1D *recoPt = new TH1D("recoPt","recoPt", 200, 0, 200);
    TH1D *recoEta = new TH1D("recoEta","recoEta", 80, -1, 1);
    TH1D *recoPhi = new TH1D("recoPhi","recoPhi", 63, 0, 6.3);

    std::cout<<"Formatting DATA file: "<<inputFile<<std::endl;
    TFile *fin = new TFile(inputFile);
    int eventCount = 0;
    // define all the column names u want to read in 
    char tracks[10] = "/O2jtrack"; 
    const char *track_name = tracks;

    // recreate the new file with vectors
    TFile *fout = TFile::Open(outputFile,"RECREATE");
    if (!fout) { return; }
    std::vector<float> particle_data_pt;
    std::vector<float> particle_data_eta;
    std::vector<float> particle_data_phi;
    std::vector<int> event_number;

    // Create a TTree //per event!
    TTree *tree = new TTree("o2Tracks","Tree with vectors");
    tree->Branch("particle_data_pt",&particle_data_pt);
    tree->Branch("particle_data_eta",&particle_data_eta); 
    tree->Branch("particle_data_phi",&particle_data_phi);
    tree->Branch("event_number",&event_number);

    particle_data_pt.clear();
    particle_data_eta.clear();
    particle_data_phi.clear();
    event_number.clear();

    // loop over the directories (DF's) of O2 - or to compress: o2-aod-merger --input new.txt
    for(auto k : *fin->GetListOfKeys()){
      const char *dir_name = k->GetName();// extract name of the current DF
      char buff_trk[100];// add name of DF and track column together
      strncpy(buff_trk, dir_name, sizeof(buff_trk));
      strncat(buff_trk, track_name, sizeof(buff_trk)-10);
      const char *tree_trk = buff_trk;
      TTree* trk = fin->Get<TTree>(tree_trk);// get the tree from the corresponding DF
      if (trk == NULL) continue;// skip the parenFiles
      Float_t pt;// declare all variable you want to extract
      Float_t eta;
      Float_t phi;
      Int_t idx = 0;// test if we can get all tracks ?
      Int_t idxColl = 0;// to find out if we are at the next collision

      trk->SetBranchAddress("fPt",&pt);
      trk->SetBranchAddress("fEta",&eta);
      trk->SetBranchAddress("fPhi",&phi);
      trk->SetBranchAddress("fIndexJCollisions",&idxColl);

      std::cout<<trk->GetEntries()<<" entries on data level !"<<std::endl;
      for(int i=0; i<trk->GetEntries(); i++){// loop over the tree (ordered by collisions)
        trk->GetEntry(i);

        if((idxColl > idx) || (i==trk->GetEntries())){// fill vectros from previous event || or the last
          event_number.push_back(idxColl);
          tree->Fill();// write vectors to tree
          idx = idxColl;
          eventCount += 1;
          particle_data_pt.clear();// clear vectors for next event
          particle_data_eta.clear();
          particle_data_phi.clear();
        }

        particle_data_pt.push_back(pt);// fill vectors with vectors
        particle_data_eta.push_back(eta);
        particle_data_phi.push_back(phi);

      }// loop over data tree    
    }
    std::cout<<"Number of events (data): "<<eventCount<<std::endl;
    fout->Write();// write tree to file

    TCanvas *reco = new TCanvas("recoQA","recoQA", 1200, 600);
    TPad *pt = new TPad("h1pt","h1pt", 0, 0,1,1);
    reco->Divide(3,0);
    reco->cd(1);
    pt->Draw();
    pt->SetLogy();
    pt->cd();
    recoPt->Scale(1/recoPt->Integral());
    recoPt->SetMarkerStyle(23);
    recoPt->Draw("Esame");
    reco->cd(2);
    recoEta->Scale(1/recoEta->Integral());
    recoEta->SetMarkerStyle(23);
    recoEta->Draw("E");
    reco->cd(3);
    recoPhi->Scale(1/recoPhi->Integral());
    recoPhi->SetMarkerStyle(23);
    recoPhi->Draw("E");
    reco->SaveAs("dataVectorQA.pdf");
    delete fout;// delete pointer / close the file

}


void readAndWriteMC(TString inputFile, TString outputFile){
  // Declare QA histos
  TH1D *recoPt = new TH1D("recoPt","recoPt", 200, 0, 200);
  TH1D *recoEta = new TH1D("recoEta","recoEta", 80, -1, 1);
  TH1D *recoPhi = new TH1D("recoPhi","recoPhi", 63, 0, 6.3);
  TH1D *truthPt = new TH1D("truthPt","truthPt", 200, 0, 200);
  TH1D *truthEta = new TH1D("truthEta","truthEta", 80, -1, 1);
  TH1D *truthPhi = new TH1D("truthPhi","truthPhi", 63, 0, 6.3);

  std::cout<<"Formatting MC file: "<<inputFile<<std::endl;
  TFile *fin = new TFile(inputFile);
  int eventCountTruth = 0;
  int eventCountReco = 0;
  //fin->ls();
  // define all the column names u want to read in 
  char recoTracks[10] = "/O2jtrack"; // reco
  char truthTracks[15] = "/O2jmcparticle"; // truth
  char truthColls[16] = "/O2jmccollision"; // truth
  char runNumber[7]  = "/O2jbc";
  
  const char *reco_name = recoTracks;
  const char *truth_name = truthTracks;
  const char *truth_collname = truthColls;
  const char *run_name = runNumber;

  // recreate the new file with vectors
  TFile *fout = TFile::Open(outputFile,"RECREATE");
  if (!fout) { return; }

  std::vector<float> particle_truth_pt;
  std::vector<float> particle_truth_eta;
  std::vector<float> particle_truth_phi;
  std::vector<int> event_numberMC;

  // Create a TTree //per event!
  TTree *tree = new TTree("o2Tracks","Tree with vectors");
  tree->Branch("particle_truth_pt",&particle_truth_pt);
  tree->Branch("particle_truth_eta",&particle_truth_eta); 
  tree->Branch("particle_truth_phi",&particle_truth_phi);
  tree->Branch("event_numberMC",&event_numberMC);

  particle_truth_pt.clear();
  particle_truth_eta.clear();
  particle_truth_phi.clear();
  event_numberMC.clear();
  // loop over truth
  for(auto k : *fin->GetListOfKeys()){
    const char *dir_name = k->GetName();// extract name of the current DF
    char buff_trkT[100];// add name of DF and track column together
    strncpy(buff_trkT, dir_name, sizeof(buff_trkT));
    strncat(buff_trkT, truth_name, sizeof(buff_trkT)-15);// truth name
    const char *tree_trkT = buff_trkT;
    TTree* trkT = fin->Get<TTree>(tree_trkT);// get the tree from the corresponding DF
    
    char buff_coll[100];
    char buff_bc[100];
    strncpy(buff_coll, dir_name, sizeof(buff_coll));
    strncat(buff_coll, truth_collname, sizeof(buff_coll)-16);
    strncat(buff_bc, run_name, sizeof(buff_bc)-7);
    const char *tree_coll = buff_coll;
    const char *tree_bc = buff_bc;
    TTree* coll = fin->Get<TTree>(tree_coll);
    TTree* bcT = fin->Get<TTree>(tree_bc);

    if (coll == NULL) continue;
    if (trkT == NULL) continue;// skip the parenFiles

    Float_t pt;// declare all variable you want to extract
    Float_t eta;
    Float_t phi;
    Float_t zColl;
    Int_t statusCode = 0; //for empty - 1 for isFinal()
    Int_t idx = 0;
    Int_t idxColl = 0;// to find out if we are at the next collision
    Float_t runNumber = 0;
    trkT->SetBranchAddress("fPt",&pt);
    trkT->SetBranchAddress("fEta",&eta);
    trkT->SetBranchAddress("fPhi",&phi);
    trkT->SetBranchAddress("fIndexJMcCollisions",&idxColl);
    trkT->SetBranchAddress("fHepMCStatusCode", &statusCode);
    coll->SetBranchAddress("fPosZ",&zColl);
    //bcT->SetBranchAddress("fRunNumber", &runNumber);// make a loop oe runnumbers to fill the histos (Th1D's) per run - repeat in data

    bool duplicatePt = false;
    bool duplicateEta = false;
    bool duplicatePhi = false;

    std::cout<<trkT->GetEntries()<<" entries on thruth level !"<<std::endl;
    for(int i=0; i<trkT->GetEntries(); i++){// loop over the tree (ordered by collisions)
      trkT->GetEntry(i);
      if (statusCode != 1) continue;// consider only fs (for now)
      
      truthPt->Fill(pt);
      truthPhi->Fill(phi);
      truthEta->Fill(eta);
      
      duplicatePt = duplicateEta = duplicatePhi = false;
      if((idxColl > idx) || (i==trkT->GetEntries())){// fill vectros from previous event || or the last
        //std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Next Event"<<std::endl;
        if(particle_truth_pt[0]>0) tree->Fill();// write vectors to tree only if there is at least one particle
        idx = idxColl;
        eventCountTruth += 1;
        particle_truth_pt.clear();// clear vectors for next event
        particle_truth_eta.clear();
        particle_truth_phi.clear();
        event_numberMC.clear();
      }


      //coll->GetEntry(idxColl);
      //if(abs(zColl) > 10) continue;// maybe worth to cut of y max as well...
      //if((abs(eta)>1) || (pt < 5) || (pt > 400)) continue;// to save some storage
      
      auto itPt = std::find(particle_truth_pt.begin(), particle_truth_pt.end(), pt);
      auto itEta = std::find(particle_truth_eta.begin(), particle_truth_eta.end(), eta);
      auto itPhi = std::find(particle_truth_phi.begin(), particle_truth_phi.end(), phi);
      // Check if the iterator points to the end of the vector
      if (itPt != particle_truth_pt.end()) duplicatePt = true;
      if (itEta != particle_truth_eta.end()) duplicateEta = true;
      if (itPhi != particle_truth_phi.end()) duplicatePhi = true;
      
      if(((duplicatePt == true) && (duplicateEta == true) && (duplicatePhi == true))){// avoid same particle in event -- requiring && is too tight -- jettoyHI crash
        std::cout << "Duplicated values found! " << std::endl;
        std::cout << "Pt  " << pt << std::endl;
        std::cout << "Eta "<< eta << std::endl;
        std::cout << "Phi "<< phi << std::endl;
        continue;
      }
      
      // fill vector for the event
      particle_truth_pt.push_back(pt);
      particle_truth_eta.push_back(eta);
      particle_truth_phi.push_back(phi);
      event_numberMC.push_back(idxColl);
      //std::cout<<"pt "<<pt<<" eta "<<eta<<" phi "<<phi<<" idxColl: "<<idxColl<<std::endl;
    }
    std::cout<<"End DF Truth "<<std::endl;
    //tree->Fill();// write vectors to tree
  }
  std::cout<<"End loop Truth "<<std::endl;

  std::vector<float> particle_reco_pt;
  std::vector<float> particle_reco_eta;
  std::vector<float> particle_reco_phi;
  std::vector<int> event_number;
  
  tree->Branch("particle_reco_pt",&particle_reco_pt);
  tree->Branch("particle_reco_eta",&particle_reco_eta); 
  tree->Branch("particle_reco_phi",&particle_reco_phi);
  tree->Branch("event_number",&event_number);

  particle_reco_pt.clear();
  particle_reco_eta.clear();
  particle_reco_phi.clear();
  event_number.clear();


  // loop over reco
  for(auto k : *fin->GetListOfKeys()){
    const char *dir_name = k->GetName();// extract name of the current DF
    char buff_trk[100];// add name of DF and track column together
    strncpy(buff_trk, dir_name, sizeof(buff_trk));
    strncat(buff_trk, reco_name, sizeof(buff_trk)-10);//start with reco
    const char *tree_trk = buff_trk;
    TTree* trk = fin->Get<TTree>(tree_trk);// get the tree from the corresponding DF
    if (trk == NULL) continue;// skip the parenFiles
    Float_t pt;// declare all variable you want to extract
    Float_t eta;
    Float_t phi;
    Int_t idx = 0;// test if we can get all tracks ?
    Int_t idxColl = 0.0;// to find out if we are at the next collision

    trk->SetBranchAddress("fPt",&pt);
    trk->SetBranchAddress("fEta",&eta);
    trk->SetBranchAddress("fPhi",&phi);
    trk->SetBranchAddress("fIndexJCollisions",&idxColl);

    std::cout<<trk->GetEntries()<<" entries on reco level !"<<std::endl;
   
    for(int i=0; i<trk->GetEntries(); i++){// loop over the tree (ordered by collisions)
      std::cout<<i<<" entry "<<std::endl;
      
      trk->GetEntry(i);
      
      recoPt->Fill(pt);
      recoPhi->Fill(phi);
      recoEta->Fill(eta);
      std::cout<<" after fill "<<std::endl;
      
      if((idxColl > idx) || (i==trk->GetEntries())){// fill vectros from previous event || or the last
        //std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Next Event"<<std::endl;
        tree->Fill();// write vectors to tree
        idx = idxColl;
        eventCountReco += 1;
        particle_reco_pt.clear();// clear vectors for next event
        particle_reco_eta.clear();
        particle_reco_phi.clear();
      }

      //if((abs(eta)>1.5) || (pt < 1) || (pt > 400)) continue;
      particle_reco_pt.push_back(pt);
      particle_reco_eta.push_back(eta);
      particle_reco_phi.push_back(phi);
      event_number.push_back(idxColl);
    }
    std::cout<<"End DF Reco "<<std::endl;
    //tree->Fill();
  }
  std::cout<<"End Loop Reco "<<std::endl;
  std::cout<<"Number of events (truth): "<<eventCountTruth<<std::endl;
  std::cout<<"Number of events (reco): "<<eventCountReco<<std::endl;

  fout->Write();// write trees to file
  delete fout;// delete pointer / close the file

  //self normalize (& then canvas per runnumber)              -- needs mini fix
  TCanvas *reco = new TCanvas("recoQA","recoQA", 1200, 600);
  TPad *pt = new TPad("h1pt","h1pt", 0, 0,1,1);
  reco->Divide(3,0);
  reco->cd(1);
  pt->Draw();
  pt->SetLogy();
  pt->cd();
  recoPt->Scale(1/recoPt->Integral());
  recoPt->SetMarkerStyle(23);
  recoPt->Draw("Esame");
  reco->cd(2);
  recoEta->Scale(1/recoEta->Integral());
  recoEta->SetMarkerStyle(23);
  recoEta->Draw("E");
  reco->cd(3);
  recoPhi->Scale(1/recoPhi->Integral());
  recoPhi->SetMarkerStyle(23);
  recoPhi->Draw("E");
  reco->SaveAs("recoVectorQA.pdf");

  TCanvas *truth = new TCanvas("truthQA","truthQA", 1200, 600);
  TPad *ptT = new TPad("h1pt","h1pt", 0, 0,1,1);
  truth->Divide(3,0);
  truth->cd(1);
  ptT->Draw();
  ptT->SetLogy();
  ptT->cd();
  truthPt->Scale(1/truthPt->Integral());
  truthPt->SetMarkerStyle(23);
  truthPt->Draw("Esame");
  truth->cd(2);
  truthEta->Scale(1/truthEta->Integral());
  truthEta->SetMarkerStyle(23);
  truthEta->Draw("E");
  truth->cd(3);
  truthPhi->Scale(1/truthPhi->Integral());
  truthPhi->SetMarkerStyle(23);
  truthPhi->Draw("E");
  reco->SaveAs("truthVectorQA.pdf");

  /*
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  +++++++++++++++++++++++++++++++++is this to be expected ??
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++S
  5599 entries on reco level !     -- 10090 events ??
  5551 entries on reco level !
  556818 entries on thruth level ! -- 10092 events ??
  621135 entries on thruth level !

  -- data has 11707 track entries  (one DF)     -- 18989 events (all DF's)
  */
}


void readData(){
    //read("AO2D.root");
    //readAndWrite("TestAO2D.root", "vectorTree_LHC22o_pass7.root");
    
   //readAndWrite("/dcache/alice/jlomker/LHC22o_pass7_minBias/346364", "/dcache/alice/jlomker/LHC22o_pass7_minBias/346364/vectorTree_LHC22o_pass7.root");
   //readAndWriteMC("AO2D_LHC25a2b.root", "vectorTree_LHC25a2b.root");
}
