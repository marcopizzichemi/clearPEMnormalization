#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "globals.hh"

class CreateTree
{
  private:

  TTree*              ftree;
  TString             fname;

  Bool_t              HITS;
  Bool_t              ABSORPTIONS;

  static const Int_t  MaxNum = 2000000;
  static const Int_t  MaxNumPro = 100;
  
  int nCrystalsX, nCrystalsY, nDetectorsX, nDetectorsY;

  public:

  CreateTree(TString name/*,int x, int y, int z, int k*/);
  ~CreateTree();

  
  
  //void SetModuleElements(int x, int y,int z,int k){ nCrystalsX = x; nCrystalsY = y ; nDetectorsX = z ; nDetectorsY = k; };
  
  
  
  
  TTree*              GetTree() const { return ftree; };
  TString             GetName() const { return fname;};
  Int_t               Fill() { return this->GetTree()->Fill(); };
  Bool_t              Write();
  void                Clear();
  static CreateTree*  Instance() { return fInstance; };
  static CreateTree*  fInstance;

  Bool_t              Hits() const { return this->HITS; };
  Bool_t              Absorptions() const { return this->ABSORPTIONS; };

  long int            Seed;
  Int_t               Run;
  Int_t               Event;
  Int_t               NumOptPhotons;
  Int_t               NumCherenkovPhotons;
  Float_t             totalEnergyDeposited;

  Float_t             enPerCrystal[6144]; 
  
  long int            EventTag;
  
  Short_t*            DetectorHit;          
//   int            DetectorHit[16];
  
  
//   std::vector<float>*  CryEnergyDeposited;   
//   std::vector<float>** pCryEnergyDeposited;
//   
  std::vector<float>  EnergyDeposited;
  std::vector<float>*  pEnergyDeposited;
  
  std::vector<int>    Crystal;
  std::vector<int>*    pCrystal;
  
  
  
  std::vector<float>  PosXEnDep; 
  std::vector<float>* pPosXEnDep;
  std::vector<float>  PosYEnDep; 
  std::vector<float>* pPosYEnDep;
  std::vector<float>  PosZEnDep; 
  std::vector<float>* pPosZEnDep;
  
  
//   std::vector<float> PositionX;
//   std::vector<float> *pPositionX;
//   std::vector<float> PositionY;
//   std::vector<float> *pPositionY;
//   
//   std::vector<float> GlobalTime;
//   std::vector<float> *pGlobalTime;
  
  
 
};
