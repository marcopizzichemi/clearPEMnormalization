#include "CreateTree.hh"
#include <sstream>
#include "TROOT.h"

CreateTree* CreateTree::fInstance = NULL;

using namespace std;

CreateTree::CreateTree(TString name)
{
  gROOT->ProcessLine("#include <vector>"); //this is needed otherwise ROOT will complain about not knowing what a std::vector is...
  
  if(fInstance) 
  {
    return;
  }
  
  this->HITS=true;
  this->ABSORPTIONS=true;
  
  this->fInstance = this;
  this->fname = name;
  this->ftree = new TTree("tree","name");
  
  //global branches
  this->GetTree()->Branch("Seed",&this->Seed,"Seed/L");
  this->GetTree()->Branch("Run",&this->Run,"Run/I");
  this->GetTree()->Branch("Event",&this->Event,"Event/I");
  this->GetTree()->Branch("GammaParity",&this->GammaParity,"GammaParity/I");
  this->GetTree()->Branch("DoEvent",&this->DoEvent,"DoEvent/I");
  
  std::stringstream snames;
  
  //branches for the first gamma
  this->GetTree()->Branch("First_totalEnergyDeposited",&this->First_totalEnergyDeposited,"First_totalEnergyDeposited/F");
  First_pEnergyDeposited = &First_EnergyDeposited;
  First_pCrystal = &First_Crystal;
  First_pPosXEnDep = &First_PosXEnDep;
  First_pPosYEnDep = &First_PosYEnDep;
  First_pPosZEnDep = &First_PosZEnDep;
  
  this->GetTree()->Branch("First_EnergyDeposited","std::vector<float>",&First_pEnergyDeposited);
  this->GetTree()->Branch("First_Crystal","std::vector<int>",&First_pCrystal);
  snames<< "First_PosXEnDep";
  this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&First_pPosXEnDep);
  snames.str("");
  snames << "First_PosYEnDep";
  this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&First_pPosYEnDep);
  snames.str("");
  snames << "First_PosZEnDep";
  this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&First_pPosZEnDep);
  
  //branches for the second gamma
  this->GetTree()->Branch("Second_totalEnergyDeposited",&this->Second_totalEnergyDeposited,"Second_totalEnergyDeposited/F");
  Second_pEnergyDeposited = &Second_EnergyDeposited;
  Second_pCrystal = &Second_Crystal;
  Second_pPosXEnDep = &Second_PosXEnDep;
  Second_pPosYEnDep = &Second_PosYEnDep;
  Second_pPosZEnDep = &Second_PosZEnDep;
  
  this->GetTree()->Branch("Second_EnergyDeposited","std::vector<float>",&Second_pEnergyDeposited);
  this->GetTree()->Branch("Second_Crystal","std::vector<int>",&Second_pCrystal);
  snames.str("");
  snames<< "Second_PosXEnDep";
  this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&Second_pPosXEnDep);
  snames.str("");
  snames << "Second_PosYEnDep";
  this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&Second_pPosYEnDep);
  snames.str("");
  snames << "Second_PosZEnDep";
  this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&Second_pPosZEnDep);
  
  this->Clear();
}

CreateTree::~CreateTree()
{
//   delete DetectorHit;
//   delete CryEnergyDeposited;
//   delete pCryEnergyDeposited;
//   delete PosXEnDep ;
//   delete pPosXEnDep;
//   delete PosYEnDep;
//   delete pPosYEnDep;
//   delete PosZEnDep;
//   delete pPosZEnDep;
}

Bool_t CreateTree::Write()
{
  TString filename = this->GetName();
  filename+=".root";
  TFile* file = new TFile(filename,"RECREATE");
  this->GetTree()->Write();
  file->Write();
  file->Close();
  return true;
}

void CreateTree::Clear()
{
  Run=0;
  Event=0;
  GammaParity=0; //so next one will start from first gamma
  DoEvent = 0;
  
  First_totalEnergyDeposited=0;
  First_EnergyDeposited.clear();
  First_Crystal.clear();
  First_PosXEnDep.clear();
  First_PosYEnDep.clear();
  First_PosZEnDep.clear();
  
  Second_totalEnergyDeposited=0;
  Second_EnergyDeposited.clear();
  Second_Crystal.clear();
  Second_PosXEnDep.clear();
  Second_PosYEnDep.clear();
  Second_PosZEnDep.clear();
  
  
}
