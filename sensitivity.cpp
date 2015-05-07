// compile with

// g++ -o ../analysis-bin/sensitivity sensitivity.cpp `root-config --cflags --glibs`

// so you need to have a ../analysis-bin/ folder ready (this way the binary is not in the folder watched by git)

// this code transforms the output of the normalization "simple" simulation into the elm2 format 
// that is usually produced by the detector (and that should be fed to norm_total_gen)
// it will produce only a listmode file with true coincidences

// the idea is to assign x and y coordinates as the coordinates of the crystal identified, while z is the 
// weighted average of the real z information given by the simulation data (weighted on energy deposited)
// being there no detector in this simulation, we can't take the DOI from anywhere else...



#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TF2.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TGraph2D.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TChain.h"
#include <vector>
#include <algorithm>

#include <stdlib.h> 
#include <stdio.h> 
#include <unistd.h>
#include <cmath> 

//declares the struct of events
struct EventFormat 
{
  double ts;				// time of the event, in seconds. if i'm not mistaken is the absolute machine time in seconds 
  u_int8_t random;                      // set probably at the level of acq, says if the event is random (when it's value is 1)
  float d;                              // distance between the heads, fixed to the head distance value... 
  float yozRot;		                // rotation of the heads 
  float x1;                             // x coordinate of the event, for detector 1		
  float y1;                             // y coordinate of the event, for detector 1
  float z1;                             // z coordinate of the event, for detector 1
  float e1;                             // energy deposited in detector 1
  u_int8_t n1;                          // some check done probably at the level of acq, it's 1 if the event is ok, otherwise is not ok
  float x2;                             // x coordinate of the event, for detector 2
  float y2;                             // y coordinate of the event, for detector 2
  float z2;                             // z coordinate of the event, for detector 2
  float e2;                             // energy deposited in detector 2
  u_int8_t n2;                          // some check done probably at the level of acq, it's 1 if the event is ok, otherwise is not ok
  float dt;                             // delta time between event in detector 1 and event in detector 2 
} __attribute__((__packed__));     



int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>");
  
  TChain *tree =  new TChain("tree");
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  
  std::string oFile = "t.dat";
  std::ofstream offs;
  offs.open (oFile.c_str(), std::ofstream::out);
  
  long int Seed;
  int Run;
  int Event;
  long int EventTag;
  float totalEnergyDeposited;
  int NumOptPhotons;
  int NumCherenkovPhotons;
  
  std::vector<float> *pEn = 0;
  std::vector<int> *pCry = 0;
  std::vector<float> *pPosXEnDep = 0;
  std::vector<float> *pPosYEnDep = 0;
  std::vector<float> *pPosZEnDep = 0;
  
//   TBranch *bvpx = 0;
  
  
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("EventTag",&EventTag);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);
  
  tree->SetBranchAddress("EnergyDeposited",&pEn);
  tree->SetBranchAddress("Crystal",&pCry);
  tree->SetBranchAddress("PosXEnDep",&pPosXEnDep);
  tree->SetBranchAddress("PosYEnDep",&pPosYEnDep);
  tree->SetBranchAddress("PosZEnDep",&pPosZEnDep);
  
  TH2F *projectionSlice = new TH2F("projection","projection", 64, -90, 90, 48, -80, 80);
  
  
  TH2F *efficiency_dh0 = new TH2F("efficiency_dh0","efficiency_dh0", 64, 0, 64, 48, 0, 48);
  TH2F *efficiency_dh1 = new TH2F("efficiency_dh1","efficiency_dh1", 64, 0, 64, 48, 0, 48);
  
  float enDH0, enDH1;
  float avgX,avgX1,avgY,avgY1,avgZ,avgZ1;
  float sX0,sX1,sY0,sY1;
  
  
  
  
  
  
  long int counter = 0;
  long int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  
  float x1,y1,z1,x2,y2,z2;
  float e1,e2;
  
  long int results[6144];
  
  for(int i = 0; i < 6144; i++)
  {
    results[i] = 0;
  }
  
//   TTree* exitTree = new TTree("tree","tree");
//   
//   exitTree->Branch("e1",&e1,"e1/F"); 	
//   exitTree->Branch("x1",&x1,"x1/F");
//   exitTree->Branch("y1",&y1,"y1/F");
//   exitTree->Branch("z1",&z1,"z1/F");
//   exitTree->Branch("e2",&e2,"e2/F"); 	
//   exitTree->Branch("x2",&x2,"x2/F");
//   exitTree->Branch("y2",&y2,"y2/F");
//   exitTree->Branch("z2",&z2,"z2/F");
  float maxX,maxY;
  float minX,minY;
  
  //fake coincidence sorting
  for(int i = 0; i < nEntries ; i=i+2)
  { 
    counter=counter+2;
//     std::cout << "get "<< std::endl;
    tree->GetEvent(i);
    
    std::vector<int> cryIDs;
    std::vector<double> energy;
    std::vector<double> doi;
    
    //check which crystals have energy deposition
    for(int j = 0;  j < pCry->size() ; j++) 
    {
      
      bool found = false;
      for(int k = 0 ; k < cryIDs.size() ; k++) //check if this crystal ID was already added to the cryIDs vector
      {
	if(pCry->at(j) == cryIDs.at(k))
	{ 
	  found = true;
	}
      }
      if(!found) //if not, add it
	cryIDs.push_back(pCry->at(j));
    }
    
    //for each crystal ID, sum the energy deposited into it
    //calculate the weighted average z position at the same time for each crystal ID
    for(int j = 0;  j < cryIDs.size() ; j++) 
    {
      energy.push_back(0);
      doi.push_back(0);
      for(int k = 0 ; k < pCry->size() ; k++)
      {
	if(pCry->at(k) == cryIDs.at(j)) //if the energy was deposited in that crystal
	{
	  energy[j] += pEn->at(k);
	  doi[j] += pEn->at(k)*pPosZEnDep->at(k);
	}
      }
    }
    
    for(int j = 0;  j < energy.size() ; j++) 
    {
      doi[j] = doi[j]/energy[j];
    }
    
    
    
    
    
    
//     for(int j = 0;  j < pCry->size() ; j++) 
//     {
//       std::cout <<  pCry->at(j) << "\t";
//     }
//     std::cout << std::endl;
//     
//     for(int j = 0;  j < pEn->size() ; j++) 
//     {
//       std::cout <<  pEn->at(j) << "\t";
//     }
//     std::cout << std::endl;
//     
//     for(int j = 0;  j < cryIDs.size() ; j++) 
//     {
//       std::cout <<  cryIDs.at(j) << "\t";
//     }
//     std::cout << std::endl;
//     
//     for(int j = 0;  j < energy.size() ; j++) 
//     {
//       std::cout <<  energy.at(j) << "\t";
//     }
//     std::cout << std::endl;
//     
    int cry1;
    bool coinc1 = false; //flag "i'm a candidate for coincidence"
    float energy1;
    float doi1;
    
    for(int j = 0 ; j < cryIDs.size() ; j++)
    {
      if(energy.at(j) > 0.5) //cannot really happen for 2 crystals, given a source of 511kev and gammas been taken one by one...
      {
	cry1 = cryIDs.at(j);
	energy1 = (float) energy.at(j);
	doi1 = (float) doi.at(j);
	coinc1 = true;
      }
    }
    
    
    //second event
    tree->GetEvent(i+1);
    
    
    cryIDs.clear();
    energy.clear();
    doi.clear()
    
    //check which crystals have energy deposition, sum the energy deposited
    for(int j = 0;  j < pCry->size() ; j++) 
    {
      
      bool found = false;
      for(int k = 0 ; k < cryIDs.size() ; k++)
      {
	if(pCry->at(j) == cryIDs.at(k))
	{ 
	  found = true;
	}
      }
      if(!found)
	cryIDs.push_back(pCry->at(j));
    }
    
    for(int j = 0;  j < cryIDs.size() ; j++) 
    {
      energy.push_back(0);
      doi.push_back(0);
      for(int k = 0 ; k < pCry->size() ; k++)
      {
	if(pCry->at(k) == cryIDs.at(j))
	{
	  energy[j] += pEn->at(k);
	  doi[j] += pEn->at(k)*pPosZEnDep->at(k);
	}
      }
    }
    
    for(int j = 0;  j < energy.size() ; j++) 
    {
      doi[j] = doi[j]/energy[j];
    }
    
    
//     for(int j = 0;  j < pCry->size() ; j++) 
//     {
//       std::cout <<  pCry->at(j) << "\t";
//     }
//     std::cout << std::endl;
//     
//     for(int j = 0;  j < pEn->size() ; j++) 
//     {
//       std::cout <<  pEn->at(j) << "\t";
//     }
//     std::cout << std::endl;
//     
//     for(int j = 0;  j < cryIDs.size() ; j++) 
//     {
//       std::cout <<  cryIDs.at(j) << "\t";
//     }
//     std::cout << std::endl;
//     
//     for(int j = 0;  j < energy.size() ; j++) 
//     {
//       std::cout <<  energy.at(j) << "\t";
//     }
//     std::cout << std::endl;
    
    int cry2;
    bool coinc2 = false;
    float energy2;
    float doi2;
    
    for(int j = 0 ; j < cryIDs.size() ; j++)
    {
      if(energy.at(j) > 0.5)
      {
	cry2 = cryIDs.at(j);
	energy2 = (float) energy.at(j);
	doi2 = (float) doi.at(j);
	coinc2 = true;
      }
    }
    
    
    //now if both the flags of coincidence are on
    if(coinc1 == true && coinc2 == true)
    {
      results[cry1]++;
      results[cry2]++;
      
      //offs << cry1 << "\t" << cry2 << std::endl;
      if(cry1 < 6144/2)  //then we are in the first head
      {
	int x,y;
	x = cry1 / 48;
	y = cry1 % 48;
	
// 	offs /*<< "cry1 dh0 "*/ << x << "\t" << y << std::endl;
	efficiency_dh0->Fill(x,y);
      }
      else //we are in the second head
      {
	cry1 = cry1 - (6144/2);
// 	std::cout << cry1  << std::endl;
	int x,y;
	x = cry1 / 48;
	y = cry1 % 48;
// 	offs /*<< "cry1 dh1 "*/ << x << "\t" << y << std::endl;
	efficiency_dh1->Fill(x,y);
      }
      
      if(cry2 < 6144/2)
      {
	int x,y;
	x = cry2 / 48;
	y = cry2 % 48;
// 	offs << /*"cry2 dh0 " <<*/ x << "\t" << y << std::endl;
	efficiency_dh0->Fill(x,y);
      }
      else
      {
	cry2 = cry2 - (6144/2);
// 	std::cout << cry2  << std::endl;
	int x,y;
	x = cry2 / 48;
	y = cry2 % 48;
// 	offs << /*"cry2 dh1 " <<*/ x << "\t" << y << std::endl;
	efficiency_dh1->Fill(x,y);
      }
    }
    
//     avgX = avgY = avgZ = 0;
// //     enDep = 0;
//     maxX = maxY = -1e6;
//     minX = minY = 1e6;
// //     std::cout << Event << "\t" << EventTag << "\t" << totalEnergyDeposited;
// //     if(totalEnergyDeposited != 0)
// //       std::cout << "\t" <<  pEn->size()<< std::endl;
// //     else
// //       std::cout << std::endl;
//     if(totalEnergyDeposited != 0)
//     {
//       
//       
//       //std::cout << "1 "<< std::endl;
//       //std::cout << totalEnergyDeposited << std::endl;
//       //std::cout << pEn->size() << std::endl;
//       for (int j = 0 ; j < pEn->size() ; j++)
//       {
// // 	std::cout << "1 "<< std::endl;
// 	avgX += pPosXEnDep->at(j) * pEn->at(j);
// 	avgY += pPosYEnDep->at(j) * pEn->at(j);
// 	avgZ += pPosZEnDep->at(j) * pEn->at(j);
// 	
// 	if(pPosXEnDep->at(j) > maxX)
// 	  maxX = pPosXEnDep->at(j);
// 	if(pPosXEnDep->at(j) < minX)
// 	  minX = pPosXEnDep->at(j);
// 	
// 	if(pPosYEnDep->at(j) > maxY)
// 	  maxY = pPosYEnDep->at(j);
// 	if(pPosYEnDep->at(j) < minY)
// 	  minY = pPosYEnDep->at(j);
// 	
//       }
//       avgX = avgX / totalEnergyDeposited;
//       avgY = avgY / totalEnergyDeposited;
//       avgZ = avgZ / totalEnergyDeposited;
//     }
//     else
//       continue;
// //     std::cout << maxX << "\t" << minX << "\t"<< maxY << "\t"<< minY << std::endl;
//     
//     if( (abs(maxX - minX) > 2.0) | (abs (maxY -minY) > 2.0) ) //gamma not contained in a "crystal"
//       continue;
// //     
//     x1 = avgX;
//     y1 = avgY;
//     z1 = avgZ;
//     
//     if(z1 < 100 && z1 > -100) //insane event
//       continue;
//     
//     e1 = totalEnergyDeposited;
    
    
//     avgX = avgY = avgZ = 0;
//     maxX = maxY = -1e6;
//     minX = minY = 1e6;
//     
//     if(totalEnergyDeposited != 0)
//     {
//       //std::cout << "2 "<< std::endl;
//       for (int j = 0 ; j < pEn->size() ; j++)
//       {
// 	avgX += pPosXEnDep->at(j) * pEn->at(j);
// 	avgY += pPosYEnDep->at(j) * pEn->at(j);
// 	avgZ += pPosZEnDep->at(j) * pEn->at(j);
// 	
// 	if(pPosXEnDep->at(j) > maxX)
// 	  maxX = pPosXEnDep->at(j);
// 	if(pPosXEnDep->at(j) < minX)
// 	  minX = pPosXEnDep->at(j);
// 	
// 	if(pPosYEnDep->at(j) > maxY)
// 	  maxY = pPosYEnDep->at(j);
// 	if(pPosYEnDep->at(j) < minY)
// 	  minY = pPosYEnDep->at(j);
// 	
//       }
//       avgX = avgX / totalEnergyDeposited;
//       avgY = avgY / totalEnergyDeposited;
//       avgZ = avgZ / totalEnergyDeposited;
//     }
//     else
//       continue;
// //     std::cout << maxX << "\t" << minX << "\t"<< maxY << "\t"<< minY << std::endl;
//     if( (abs(maxX - minX) > 2.0) | (abs (maxY -minY) > 2.0) ) //gamma not contained in a "crystal"
//       continue;
//     
//     x2 = avgX;
//     y2 = avgY;
//     z2 = avgZ;
//     
//     if(z2 < 100 && z2 > -100) //insane event
//       continue;
//     
//     e2 = totalEnergyDeposited;
    
//   exitTree->Fill();
    
    
    
//     std::cout << "det0" << std::endl;
//     //on detector 0
//     if(pDH0->size() != 0)
//     {
//       std::cout << "size 0 = " <<pDH0->size() <<  std::endl;
//       for (int j = 0 ; j < pDH0->size() ; j++)
//       {
// 	avgX0 += pDH0PosXEnDep->at(j) / pDH0PosXEnDep->size();
// 	avgY0 += pDH0PosYEnDep->at(j) / pDH0PosYEnDep->size();
// 	avgZ0 += pDH0PosZEnDep->at(j) / pDH0PosZEnDep->size();
// 	enDH0 += pDH0->at(j);
//       }
//     }
//     else
//     {
//       enDH0 = 0;
//     }
//     
//     std::cout << "det1" << std::endl;
//     // on detector 1
//     if(pDH1->size() != 0)
//     {
//       for (int j = 0 ; j < pDH1->size() ; j++)
//       {
// 	avgX1 += pDH1PosXEnDep->at(j) / pDH1PosXEnDep->size();
// 	avgY1 += pDH1PosYEnDep->at(j) / pDH1PosYEnDep->size();
// 	avgZ1 += pDH1PosZEnDep->at(j) / pDH1PosZEnDep->size();
// 	enDH1 += pDH1->at(j);
//       }
//     }
//     else
//     {
//       enDH1 = 0;
//     }
//     
//     if( (e1 > 0.5) && (e2 > 0.5) )
//     {
//       //calculates the coordinate u,v,w (distance on x,y,z for events in detector 1 and 2)
//       float u = x1 - x2;
//       float v = y1 - y2;
//       float w = z1 - z2;
//       
//       //some opearations for calculating the intersection of LORs with z=0
//       float s = (0 - z1)/w;
//       float mx = x1 + s*u;
//       float my = y1 + s*v;
//       
//       projectionSlice->Fill(mx, my);
//     }
    
    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
    
    
  }
  std::cout << std::endl;
 
  
//   std::string outFile = "test.dat";
//   std::ofstream ofs;
//   ofs.open (outFile.c_str(), std::ofstream::out);
//   
//   for(int i = 0; i < 6144; i++)
//   {
//     if(((i % 64) == 0) && (i!=0) )
//       ofs << std::endl;
//     ofs << results[i] << "\t";
//     
//   }
//   ofs.close();
  
  
  std::string effFile = "eff_" + std::string(argv[1]);
  TFile* fOut = new TFile(effFile.c_str(),"recreate");
  efficiency_dh0->Write();
  efficiency_dh1->Write();
  fOut->Close();
  
//   std::string outFile = "coincidence_" + std::string(argv[1]);
//   TFile* fOut = new TFile(outFile.c_str(),"recreate");
//   projectionSlice->Write();  
//   fOut->Close();
//   
//   projectionSlice->RebinX(2);
//   projectionSlice->RebinY(2);
//   
//   
//   TCanvas *Canvas = new TCanvas("","",1200,800);
//   std::string outplot = "plot_" + std::string(argv[1]) + ".gif";
//   projectionSlice->Draw("COLZ");
//   Canvas->Print(outplot.c_str());
  
  return 0;
}