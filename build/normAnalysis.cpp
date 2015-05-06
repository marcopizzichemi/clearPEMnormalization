// g++ -o normAnalysis normAnalysis.cpp `root-config --cflags --glibs`


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


int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>");
  
  TChain *tree =  new TChain("tree");
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  
  long int Seed;
  int Run;
  int Event;
  long int EventTag;
  float totalEnergyDeposited;
  int NumOptPhotons;
  int NumCherenkovPhotons;
  
  std::vector<float> *pEn = 0;  
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
  tree->SetBranchAddress("PosXEnDep",&pPosXEnDep);
  tree->SetBranchAddress("PosYEnDep",&pPosYEnDep);
  tree->SetBranchAddress("PosZEnDep",&pPosZEnDep);
  
  TH2F *projectionSlice = new TH2F("projection","projection", 48, -80, 80, 64, -90, 90);
  
  float enDH0, enDH1;
  float avgX,avgX1,avgY,avgY1,avgZ,avgZ1;
  float sX0,sX1,sY0,sY1;
  
  long int counter = 0;
  long int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  
  float x1,y1,z1,x2,y2,z2;
  float e1,e2;
  
  
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
  for(int i = 0; i < nEntries ; i++)
  { 
    counter++;
//     std::cout << "get "<< std::endl;
    tree->GetEvent(i);
    avgX = avgY = avgZ = 0;
//     enDep = 0;
    maxX = maxY = -1e6;
    minX = minY = 1e6;
//     std::cout << Event << "\t" << EventTag << "\t" << totalEnergyDeposited;
//     if(totalEnergyDeposited != 0)
//       std::cout << "\t" <<  pEn->size()<< std::endl;
//     else
//       std::cout << std::endl;
    if(totalEnergyDeposited != 0)
    {
      
      
      //std::cout << "1 "<< std::endl;
      //std::cout << totalEnergyDeposited << std::endl;
      //std::cout << pEn->size() << std::endl;
      for (int j = 0 ; j < pEn->size() ; j++)
      {
// 	std::cout << "1 "<< std::endl;
	avgX += pPosXEnDep->at(j) * pEn->at(j);
	avgY += pPosYEnDep->at(j) * pEn->at(j);
	avgZ += pPosZEnDep->at(j) * pEn->at(j);
	
	if(pPosXEnDep->at(j) > maxX)
	  maxX = pPosXEnDep->at(j);
	if(pPosXEnDep->at(j) < minX)
	  minX = pPosXEnDep->at(j);
	
	if(pPosYEnDep->at(j) > maxY)
	  maxY = pPosYEnDep->at(j);
	if(pPosYEnDep->at(j) < minY)
	  minY = pPosYEnDep->at(j);
	
      }
      avgX = avgX / totalEnergyDeposited;
      avgY = avgY / totalEnergyDeposited;
      avgZ = avgZ / totalEnergyDeposited;
    }
    else
      continue;
//     std::cout << maxX << "\t" << minX << "\t"<< maxY << "\t"<< minY << std::endl;
    
    if( (abs(maxX - minX) > 2.0) | (abs (maxY -minY) > 2.0) ) //gamma not contained in a "crystal"
      continue;
//     
    x1 = avgX;
    y1 = avgY;
    z1 = avgZ;
    
    if(z1 < 100 && z1 > -100) //insane event
      continue;
    
    e1 = totalEnergyDeposited;
    
    tree->GetEvent(i+1);
    avgX = avgY = avgZ = 0;
    maxX = maxY = -1e6;
    minX = minY = 1e6;
    
    if(totalEnergyDeposited != 0)
    {
      //std::cout << "2 "<< std::endl;
      for (int j = 0 ; j < pEn->size() ; j++)
      {
	avgX += pPosXEnDep->at(j) * pEn->at(j);
	avgY += pPosYEnDep->at(j) * pEn->at(j);
	avgZ += pPosZEnDep->at(j) * pEn->at(j);
	
	if(pPosXEnDep->at(j) > maxX)
	  maxX = pPosXEnDep->at(j);
	if(pPosXEnDep->at(j) < minX)
	  minX = pPosXEnDep->at(j);
	
	if(pPosYEnDep->at(j) > maxY)
	  maxY = pPosYEnDep->at(j);
	if(pPosYEnDep->at(j) < minY)
	  minY = pPosYEnDep->at(j);
	
      }
      avgX = avgX / totalEnergyDeposited;
      avgY = avgY / totalEnergyDeposited;
      avgZ = avgZ / totalEnergyDeposited;
    }
    else
      continue;
//     std::cout << maxX << "\t" << minX << "\t"<< maxY << "\t"<< minY << std::endl;
    if( (abs(maxX - minX) > 2.0) | (abs (maxY -minY) > 2.0) ) //gamma not contained in a "crystal"
      continue;
    
    x2 = avgX;
    y2 = avgY;
    z2 = avgZ;
    
    if(z2 < 100 && z2 > -100) //insane event
      continue;
    
    e2 = totalEnergyDeposited;
    
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
    if( (e1 > 0.5) && (e2 > 0.5) )
    {
      //calculates the coordinate u,v,w (distance on x,y,z for events in detector 1 and 2)
      float u = x1 - x2;
      float v = y1 - y2;
      float w = z1 - z2;
      
      //some opearations for calculating the intersection of LORs with z=0
      float s = (0 - z1)/w;
      float mx = x1 + s*u;
      float my = y1 + s*v;
      
      projectionSlice->Fill(mx, my);
    }
    
    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
    
    
  }
  
  std::string outFile = "coincidence_" + std::string(argv[1]);
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  projectionSlice->Write();  
  fOut->Close();
  
  projectionSlice->RebinX(2);
  projectionSlice->RebinY(2);
  
  
  TCanvas *Canvas = new TCanvas("","",1200,800);
  std::string outplot = "plot_" + std::string(argv[1]) + ".gif";
  projectionSlice->Draw("COLZ");
  Canvas->Print(outplot.c_str());
  
  return 0;
}