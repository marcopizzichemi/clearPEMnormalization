#include "CreateTree.hh"
#include <sstream>
#include "TROOT.h"

CreateTree* CreateTree::fInstance = NULL;

using namespace std;

CreateTree::CreateTree(TString name/*, int x, int y, int z, int k*/)
//   : nCrystalsX(x),
//     nCrystalsY(y),
//     nDetectorsX(z), 
//     nDetectorsY(k)
{
  
  //G4cout << "DETECTORS "  << " " << nDetectorsX << " " << nDetectorsY << G4endl; ;
  
  //create the vectors
//   DetectorHit         = new Short_t [nDetectorsX*nDetectorsY];
//   CryEnergyDeposited  = new std::vector<float> [nCrystalsX*nCrystalsY];
//   pCryEnergyDeposited = new std::vector<float>* [nCrystalsX*nCrystalsY];
//   PosXEnDep           = new std::vector<float> [nCrystalsX*nCrystalsY];
//   pPosXEnDep          = new std::vector<float>* [nCrystalsX*nCrystalsY];
//   PosYEnDep           = new std::vector<float> [nCrystalsX*nCrystalsY];
//   pPosYEnDep          = new std::vector<float>* [nCrystalsX*nCrystalsY];
//   PosZEnDep           = new std::vector<float> [nCrystalsX*nCrystalsY];
//   pPosZEnDep          = new std::vector<float>* [nCrystalsX*nCrystalsY];
  
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
  
  
  this->GetTree()->Branch("Seed",&this->Seed,"Seed/L");
  this->GetTree()->Branch("Run",&this->Run,"Run/I");
  this->GetTree()->Branch("Event",&this->Event,"Event/I");
  this->GetTree()->Branch("EventTag",&this->EventTag,"EventTag/L");
  this->GetTree()->Branch("totalEnergyDeposited",&this->totalEnergyDeposited,"totalEnergyDeposited/F");
  this->GetTree()->Branch("NumOptPhotons",&this->NumOptPhotons,"NumOptPhotons/I");
  this->GetTree()->Branch("NumCherenkovPhotons",&this->NumCherenkovPhotons,"NumCherenkovPhotons/I");
  this->GetTree()->Branch("FirstIsCoincidenceCandidate",&this->FirstIsCoincidenceCandidate,"FirstIsCoincidenceCandidate/I");
  
  pEnergyDeposited = &EnergyDeposited;
  this->GetTree()->Branch("EnergyDeposited","std::vector<float>",&pEnergyDeposited);
//   this->GetTree()->Branch("DH1EnergyDeposited",&this->DH1EnergyDeposited,"DH1EnergyDeposited/F");
  
  pCrystal = &Crystal;
  this->GetTree()->Branch("Crystal","std::vector<int>",&pCrystal);
  
  
  
//   
//   for(int i = 0; i < nCrystalsX*nCrystalsY ; i++) 
//   {
//     pCryEnergyDeposited[i] = &CryEnergyDeposited[i];
    pPosXEnDep = &PosXEnDep;
    pPosYEnDep = &PosYEnDep;
    pPosZEnDep = &PosZEnDep;
//   }
//   
//   for (int i = 0 ; i < nCrystalsX*nCrystalsY ; i++) 
//   {
    std::stringstream snames;
//     snames << "cry" << i;
//     this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&pCryEnergyDeposited[i]);
//     snames.str("");
    snames<< "PosXEnDep";
    this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&pPosXEnDep);
    snames.str("");
    snames << "PosYEnDep";
    this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&pPosYEnDep);
    snames.str("");
    snames << "PosZEnDep";
    this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&pPosZEnDep);
    
    
    
//   }
//   
//   for (int i = 0 ; i < nDetectorsX*nDetectorsY ; i++) 
//   {
//     std::stringstream snames,stypes;
//     snames << "detector" << i;
//     stypes << "detector" << i << "/S"; 
//     this->GetTree()->Branch(snames.str().c_str(),&DetectorHit[i],stypes.str().c_str());
//   }
//   
//   // 	pTransmissionX = &TransmissionX;
//   // 	pTransmissionY = &TransmissionY;
//   // 	pTransmissionZ = &TransmissionZ;
//   // 	
//   // 	this->GetTree()->Branch("TransmissionX","std::vector<float>",&pTransmissionX);
//   // 	this->GetTree()->Branch("TransmissionY","std::vector<float>",&pTransmissionY);
//   // 	this->GetTree()->Branch("TransmissionZ","std::vector<float>",&pTransmissionZ);
//   // 	
//   pPositionX = &PositionX;
//   pPositionY = &PositionY;
//   
//   this->GetTree()->Branch("PositionX","std::vector<float>",&pPositionX);
//   this->GetTree()->Branch("PositionY","std::vector<float>",&pPositionY);
//   
//   pGlobalTime = &GlobalTime;
//   
//   this->GetTree()->Branch("GlobalTime","std::vector<float>",&pGlobalTime);
  
  
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
  NumOptPhotons=0;
  EventTag = 0;
  NumCherenkovPhotons=0;
  totalEnergyDeposited=0;
  FirstIsCoincidenceCandidate=0;
  EnergyDeposited.clear();
  Crystal.clear();
  
//   for(int i = 0; i < nCrystalsX*nCrystalsY ; i++)
//   {
//     CryEnergyDeposited[i].clear();
  PosXEnDep.clear();
  PosYEnDep.clear();
  PosZEnDep.clear();
  
  
   
    
//   }
//   
//   for (int i = 0 ; i < nDetectorsX*nDetectorsY ; i++)//
//   {
//     DetectorHit[i] = 0;
//   }
//   
//   // 	TransmissionX.clear();
//   // 	TransmissionY.clear();
//   // 	TransmissionZ.clear();
//   
//   PositionX.clear();
//   PositionY.clear();
//   
//   GlobalTime.clear();
  
}
