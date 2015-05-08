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
  
  TTree*              GetTree() const { return ftree; };
  TString             GetName() const { return fname;};
  Int_t               Fill() { return this->GetTree()->Fill(); };
  Bool_t              Write();
  void                Clear();
  static CreateTree*  Instance() { return fInstance; };
  static CreateTree*  fInstance;

  Bool_t              Hits() const { return this->HITS; };
  Bool_t              Absorptions() const { return this->ABSORPTIONS; };

  //global variables
  long int            Seed;
  Int_t               Run;
  Int_t               Event;
  Int_t               GammaParity; // 0 if it's First gamma, 1 if it's Second
  
  //variables for first photon 
  Float_t             First_totalEnergyDeposited;
  std::vector<float>  First_EnergyDeposited;
  std::vector<float>* First_pEnergyDeposited;
  std::vector<int>    First_Crystal;
  std::vector<int>*   First_pCrystal;
  std::vector<float>  First_PosXEnDep; 
  std::vector<float>* First_pPosXEnDep;
  std::vector<float>  First_PosYEnDep; 
  std::vector<float>* First_pPosYEnDep;
  std::vector<float>  First_PosZEnDep; 
  std::vector<float>* First_pPosZEnDep;
  
  
  //variables for first photon 
  Float_t             Second_totalEnergyDeposited;
  std::vector<float>  Second_EnergyDeposited;
  std::vector<float>* Second_pEnergyDeposited;
  std::vector<int>    Second_Crystal;
  std::vector<int>*   Second_pCrystal;
  std::vector<float>  Second_PosXEnDep; 
  std::vector<float>* Second_pPosXEnDep;
  std::vector<float>  Second_PosYEnDep; 
  std::vector<float>* Second_pPosYEnDep;
  std::vector<float>  Second_PosZEnDep; 
  std::vector<float>* Second_pPosZEnDep;
  
  
 
};
