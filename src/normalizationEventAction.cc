//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: normalizationEventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
// 
/// \file normalizationEventAction.cc
/// \brief Implementation of the normalizationEventAction class

#include "normalizationEventAction.hh"
#include "normalizationRunAction.hh"
// #include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include "CreateTree.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

normalizationEventAction::normalizationEventAction(normalizationRunAction* runAction)
 : G4UserEventAction(),
      fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

normalizationEventAction::~normalizationEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void normalizationEventAction::BeginOfEventAction(const G4Event* event)
{  

  //print event number
  G4int eventID = event->GetEventID();
  G4cout << "//***********************************************//" << G4endl;
  G4cout << "//              Event Number " << eventID		<< G4endl;
  G4cout << "//***********************************************//" << G4endl;

  
  Int_t run = CreateTree::Instance()->Run;
  
  long int eventTag = CreateTree::Instance()->EventTag;
  long int seed = CreateTree::Instance()->Seed;
  CreateTree::Instance()->Clear();
  CreateTree::Instance()->Seed = seed;
  CreateTree::Instance()->Run = run;
  CreateTree::Instance()->Event = event->GetEventID();
  if(eventID % 2 == 0)
  {
    if(eventID != 0)
      CreateTree::Instance()->EventTag = eventTag+1;
    else 
      CreateTree::Instance()->EventTag = 0;
  }
  else 
  {
    CreateTree::Instance()->EventTag = eventTag;
  }
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void normalizationEventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //
  CreateTree::Instance()->Fill();
  //CreateTree::Instance()->Clear();
  
  G4cout << "Total Energy deposited in this event = "<< CreateTree::Instance()->totalEnergyDeposited << " MeV" << G4endl;
  
  
  
  //BEGIN of debug output - not a good idea, works only for 64 crystals - 16 detectors, otherwise gives seg fault. all info is anyway in the root file...
//   //get total number of crystals
//   int totCryNum = sizeof(CreateTree::Instance()->pCryEnergyDeposited)/sizeof(CreateTree::Instance()->pCryEnergyDeposited[0]);
//   //get total number of detectors
//   int totDetNum =  sizeof(CreateTree::Instance()->DetectorHit)/sizeof(CreateTree::Instance()->DetectorHit[0]);
//   
//   G4cout << "sizeof crystals "<<  sizeof(CreateTree::Instance()->pCryEnergyDeposited)/sizeof(CreateTree::Instance()->pCryEnergyDeposited[0]) << G4endl;
//   G4cout << "sizeof detectors "<<  sizeof(CreateTree::Instance()->DetectorHit)/sizeof(CreateTree::Instance()->DetectorHit[0]) << G4endl;
//   
//   for (int i = 0 ; i < 64 ; i++) //FIXME generilize to NxM crystals
//   {
//     if(CreateTree::Instance()->pCryEnergyDeposited[i]->size() != 0)
//     {
//       G4cout << "Crystal" << i;
//       for(int j = 0 ; j < CreateTree::Instance()->pCryEnergyDeposited[i]->size() ; j++)
//       {
// 	G4cout << "\t" << CreateTree::Instance()->pCryEnergyDeposited[i]->at(j);
//       }
//       G4cout << G4endl;
//       G4cout << "PosX ";
//       for(int j = 0 ; j < CreateTree::Instance()->pPosXEnDep[i]->size() ; j++)
//       {
// 	G4cout << "\t" << CreateTree::Instance()->pPosXEnDep[i]->at(j);
//       }
//       G4cout << G4endl;
//       G4cout << "PosY ";
//       for(int j = 0 ; j < CreateTree::Instance()->pPosYEnDep[i]->size() ; j++)
//       {
// 	G4cout << "\t" << CreateTree::Instance()->pPosYEnDep[i]->at(j);
//       }
//       G4cout << G4endl;
//       G4cout << "PosZ ";
//       for(int j = 0 ; j < CreateTree::Instance()->pPosZEnDep[i]->size() ; j++)
//       {
// 	G4cout << "\t" << CreateTree::Instance()->pPosZEnDep[i]->at(j);
//       }
//       G4cout << G4endl;
//     }
//   }
//   //output Detector hits
//   for(int i = 0 ; i < 16 ; i++)//FIXME generilize to NxM detectors
//   {
//     G4cout << "Detector " << i << " = " << CreateTree::Instance()->DetectorHit[i] << G4endl;
//   }
  //END of debug output
  
} 
