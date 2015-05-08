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
  gROOT->ProcessLine("#include <vector>"); //needed to be able to use std::vectors into TTree
  
  //store the first input file name
  std::string filename = argv[1];
  
  //open the input file(s) using a TChain
  TChain *tree =  new TChain("tree");
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  
  //output file in text format - FIXME clean this
//   std::string oFile = "t.dat";
//   std::ofstream offs;
//   offs.open (oFile.c_str(), std::ofstream::out);
  
  //variables to read the TTree file
  // N.B. The structure of the input file is
  // first part some data relative to the entire run (i.e. to what happened for this SINGLE gamma)
  long int Seed;
  int Run;
  int Event;
  long int EventTag;
  float totalEnergyDeposited;
  int NumOptPhotons;
  int NumCherenkovPhotons;
  // second part, some vectors carrying information related to the children of the original gamma
  // if the gamma has deposited energy, this could have happened in one shot or in many different steps
  // so each characteristic is recorded as an entry in the relative std::vector below
  // because of the way TTree are done, these info have to be retreived via pointers
  std::vector<float> *pEn = 0;              // energy deposited in this step
  std::vector<int>   *pCry = 0;             // crystal where the energy was deposited (i.e. the crystal ID, from 0 to 6143, in both the detectors)
  std::vector<float> *pPosXEnDep = 0;       // global x position of this energy deposition
  std::vector<float> *pPosYEnDep = 0;       // global y position of this energy deposition
  std::vector<float> *pPosZEnDep = 0;       // global z position of this energy deposition
  // therefore, if for example the energy of this gamma was deposited in, say, two crystals (ex. 34 and 46), and in 3 different position of the first 
  // crystal and 2 positions in the second, this event will have the above mentioned vectors that would look like follows
  
  // pEn[5]        = 0.3  , 0.05  , 0.1  , 0.05  , 0.011 
  // pCry[5]       = 34   , 34    , 34   , 45    , 45 
  // pPosXEnDep[5] = 4.5  , 4.52  , 4.51 , 6.2   , 6.21
  // pPosYEnDep[5] = -2.3 , -2.32 , -2.3 , -2.2  , -2.21
  // pPosZEnDep[5] = -120 , -122  , -127 , -122  , -123
  
  // so first energy deposition was 0.3 Mev, in cristal 34, in x = 4.5 , y = -2.3, z = -120 and so on.
  // then next event could have energy deposited in any other number of steps, so the length of the std::vectors will be different
  // (which is why we use std::vectors in the first place)
  
  //set branch addresses
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
  
  //some histograms that will be plotted at the same time...
  TH2F *projectionSlice = new TH2F("projection","projection", 64, -90, 90, 48, -80, 80);
  TH2F *efficiency_dh0 = new TH2F("efficiency_dh0","efficiency_dh0", 64, 0, 64, 48, 0, 48);
  TH2F *efficiency_dh1 = new TH2F("efficiency_dh1","efficiency_dh1", 64, 0, 64, 48, 0, 48);
  
  long int counter = 0;
  long int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  //float x1,y1,z1,x2,y2,z2;
  //float e1,e2;
  
  //binary output: create the name based on the first input file name
  std::string outFileName;
  outFileName = filename.substr(0,filename.length() -5 );
  outFileName += ".elm2";
  ofstream output_file(outFileName.c_str(), std::ios::binary);
  
  
  u_int8_t randomSet = 0;
  u_int8_t n1Set = 1;
  u_int8_t n2Set = 1;
  Float_t distance = 209.059998; //real heads distance, in mm //FIXME hardcoded for now
  Float_t yozRot = 0.0; //head rotation, in radiands          //again hardcoded but this is for normalization, so..
  Float_t absTime = 0.0; //t_start of this take, from the beginning of the "exam"
  
  //   if(argc > 2)
  //     distance = atof(argv[2]);
  //   if(argc > 3)
  //     yozRot = atof(argv[3]);
  //   if(argc > 4)
  //     absTime = atof(argv[4]);
  //   
  
//   long int results[6144];
//   
//   for(int i = 0; i < 6144; i++)
//   {
//     results[i] = 0;
//   }
//   
  //   float maxX,maxY;
  //   float minX,minY;
  
  //fake coincidence sorting
  for(int i = 0; i < nEntries /* TEMPORARY LIMITATION FOR CHECKING THE PROGRAM  100000 */; i=i+2) //read the event 2 by 2, as they are coupled like this (gamma and back to back gamma)
  { 
    counter=counter+2;
    
    EventFormat fe; // prepare the output stucture 
    //global values
    fe.ts = absTime + i*(1e-9); //not really important here
    fe.random = randomSet;
    fe.d = distance;
    fe.yozRot = yozRot;
    fe.n1 = n1Set;
    fe.n2 = n2Set;
    fe.dt = 1e-9; // set it to 1ns, so they will always be considered trues
    
    //FIRST GAMMA SHOT BY THE SIMULATION 
    //get the event data saved in the input Root TTree
    tree->GetEvent(i); 
    //prepare some useful vector
    std::vector<int> cryIDs;              //vector of the IDs of crystals where, in this event (= because of this gamma), there has been energy deposition
    std::vector<float> energy;            //vector of the energies deposited in each crystal
    std::vector<float> avgX,avgY,avgZ;    //vector of the average energy deposition position
    
    //check which crystals have energy deposition
    for(int j = 0;  j < pCry->size() ; j++) //run on the length of the crystal ID std::vector, which is the same length of all the other std::vectors (see above)
    {
      // now the vector pCry has as many elements as were the energy deposition "events", but most likely means it contains repetitions of some IDs 
      // whenever energy was deposited more than once per crystal
      // here we want to produce another std::vector of the cry IDs without repetitions
      bool found = false;
      for(int k = 0 ; k < cryIDs.size() ; k++) //check if this crystal ID was already added to the cryIDs vector (the first time, j=0, of course the program will not enter this cycle)
      {
	if(pCry->at(j) == cryIDs.at(k))
	{ 
	  found = true;
	}
      }
      if(!found) //if not, add it (and the first time, j=0, this will always be the case)
	cryIDs.push_back(pCry->at(j));
    }
    
    //for each crystal ID found, sum the energy deposited into that crystal
    //calculate the weighted average x,y,z positions at the same time for each crystal ID
    for(int j = 0;  j < cryIDs.size() ; j++) //cycle on the IDs found before
    {
      //initialize the vectors
      energy.push_back(0);
      avgX.push_back(0);
      avgY.push_back(0);
      avgZ.push_back(0);
      //now for this ID, run on all the entries of the input std::vectors (which means on each deposition event)...
      for(int k = 0 ; k < pCry->size() ; k++)
      {
	if(pCry->at(k) == cryIDs.at(j)) //...and if the energy was deposited in this crystal...
	{
	  energy[j] += pEn->at(k);                      // ... sum it to the total energy deposited in this crystal...
	  avgX[j]   += pEn->at(k)*pPosXEnDep->at(k);    // ... update the weighted average position
	  avgY[j]   += pEn->at(k)*pPosYEnDep->at(k);    // ... update the weighted average position
	  avgZ[j]   += pEn->at(k)*pPosZEnDep->at(k);    // ... update the weighted average position
	}
      }
    }
    
    // now the weighted averages still need to be divided by the total energy in order to have the correct x,y,z values
    for(int j = 0;  j < cryIDs.size() ; j++) 
    {
      avgX[j] = avgX[j] / energy[j];
      avgY[j] = avgY[j] / energy[j];
      avgZ[j] = avgZ[j] / energy[j];
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
    
    //now prepare the event for being written in the elm2 file and to enter the histograms
    int cry1;
    bool coinc1 = false; //flag "i'm a candidate for coincidence"
    float energy1;
    
    for(int j = 0 ; j < cryIDs.size() ; j++)
    {
      if(energy.at(j) > 0.5) //cannot really happen for 2 crystals, given a source of 511kev and gammas been taken one by one...
      {
	cry1 = cryIDs.at(j);
	energy1 = energy.at(j);
	coinc1 = true;
	fe.x1 = avgX.at(j);
	fe.y1 = avgY.at(j);
	fe.z1 = avgZ.at(j);
	fe.e1 = energy.at(j) * 1000.0;  //energies in the ClearPEM software are in KeV, in the simulation in MeV
      }
    }
    
    
    //SECOND EVENT 
    tree->GetEvent(i+1);
    
    if(coinc1 == true) // if the first event was not a good candidate, it makes no sense to check the second in this condition
    {
      //clear the vectors
      cryIDs.clear();
      energy.clear();
      avgX.clear();
      avgY.clear();
      avgZ.clear();
      
      //then do exactly as before...
      
      //check which crystals have energy deposition
      for(int j = 0;  j < pCry->size() ; j++) //run on the length of the crystal ID std::vector, which is the same length of all the other std::vectors (see above)
      {
	// now the vector pCry has as many elements as were the energy deposition "events", but most likely means it contains repetitions of some IDs 
	// whenever energy was deposited more than once per crystal
	// here we want to produce another std::vector of the cry IDs without repetitions
	bool found = false;
	for(int k = 0 ; k < cryIDs.size() ; k++) //check if this crystal ID was already added to the cryIDs vector (the first time, j=0, of course the program will not enter this cycle)
	{
	  if(pCry->at(j) == cryIDs.at(k))
	  { 
	    found = true;
	  }
	}
	if(!found) //if not, add it (and the first time, j=0, this will always be the case)
	  cryIDs.push_back(pCry->at(j));
      }
      
      //for each crystal ID found, sum the energy deposited into that crystal
      //calculate the weighted average x,y,z positions at the same time for each crystal ID
      for(int j = 0;  j < cryIDs.size() ; j++) //cycle on the IDs found before
      {
	//initialize the vectors
	energy.push_back(0);
	avgX.push_back(0);
	avgY.push_back(0);
	avgZ.push_back(0);
	//now for this ID, run on all the entries of the input std::vectors (which means on each deposition event)...
	for(int k = 0 ; k < pCry->size() ; k++)
	{
	  if(pCry->at(k) == cryIDs.at(j)) //...and if the energy was deposited in this crystal...
	  {
	    energy[j] += pEn->at(k);                      // ... sum it to the total energy deposited in this crystal...
	    avgX[j]   += pEn->at(k)*pPosXEnDep->at(k);    // ... update the weighted average position
	    avgY[j]   += pEn->at(k)*pPosYEnDep->at(k);    // ... update the weighted average position
	    avgZ[j]   += pEn->at(k)*pPosZEnDep->at(k);    // ... update the weighted average position
	  }
	}
      }
      
      // now the weighted averages still need to be divided by the total energy in order to have the correct x,y,z values
      for(int j = 0;  j < cryIDs.size() ; j++) 
      {
	avgX[j] = avgX[j] / energy[j];
	avgY[j] = avgY[j] / energy[j];
	avgZ[j] = avgZ[j] / energy[j];
      }
      
      //       //check which crystals have energy deposition, sum the energy deposited
      //       for(int j = 0;  j < pCry->size() ; j++) 
      //       {
      // 	
      // 	bool found = false;
      // 	for(int k = 0 ; k < cryIDs.size() ; k++)
      // 	{
      // 	  if(pCry->at(j) == cryIDs.at(k))
      // 	  { 
      // 	    found = true;
      // 	  }
      // 	}
      // 	if(!found)
      // 	  cryIDs.push_back(pCry->at(j));
      //       }
      //       
      //       for(int j = 0;  j < cryIDs.size() ; j++) 
      //       {
      // 	energy.push_back(0);
      // 	doi.push_back(0);
      // 	for(int k = 0 ; k < pCry->size() ; k++)
      // 	{
      // 	  if(pCry->at(k) == cryIDs.at(j))
      // 	  {
      // 	    energy[j] += pEn->at(k);
      // 	    doi[j] += pEn->at(k)*pPosZEnDep->at(k);
      // 	  }
      // 	}
      //       }
      //       
      //       for(int j = 0;  j < energy.size() ; j++) 
      //       {
      // 	doi[j] = doi[j]/energy[j];
      //       }
      //       
      
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
      
      for(int j = 0 ; j < cryIDs.size() ; j++)
      {
	if(energy.at(j) > 0.5) //cannot really happen for 2 crystals, given a source of 511kev and gammas been taken one by one...
	{
	  cry2 = cryIDs.at(j);
	  energy2 = energy.at(j);
	  coinc2 = true;
	  fe.x2 = avgX.at(j);
	  fe.y2 = avgY.at(j);
	  fe.z2 = avgZ.at(j);
	  fe.e2 = energy.at(j) * 1000.0;  //energies in the ClearPEM software are in KeV, in the simulation in MeV
	}
      }
      
      
      //now if both the flags of coincidence are true
      if(/*coinc1 == true &&*/ coinc2 == true)
      {
// 	results[cry1]++;
// 	results[cry2]++;
	
	// first the hitogram part, not so interesting...
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
	
	// now if this is a coincidence candidate and the events are in different detectors 
	// (we still check this but it's probably not necessary)
	if(fe.z1*fe.z2 < 0) // i.e. the two z coordinates have different signs
        {
	  output_file.write((char*)&fe,sizeof(fe)); //store the event in the elm2 file
	}
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
  output_file.close();
  
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