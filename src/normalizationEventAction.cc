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

  //store the info that would be destroied by the clear command, so you can use it immediately later
  Int_t run = CreateTree::Instance()->Run;
  long int eventTag = CreateTree::Instance()->EventTag;
  long int seed = CreateTree::Instance()->Seed;
  int FirstIsCoincidenceCandidate = CreateTree::Instance()->FirstIsCoincidenceCandidate;
  
  //clear the create tree
  CreateTree::Instance()->Clear();
  
  //these remained stored in the same way
  CreateTree::Instance()->Seed = seed;
  CreateTree::Instance()->Run = run;
  //store the event number
  CreateTree::Instance()->Event = event->GetEventID();
  //store the event tag (pairs)
  if(eventID % 2 == 0)  //this are the "first" events in the pairs
  {
    if(eventID != 0)   // if it's not the very first event
      CreateTree::Instance()->EventTag = eventTag+1; // increase the event tag +1
    else 
      CreateTree::Instance()->EventTag = 0; // otherwise start from 0
  }
  else 
  {
    CreateTree::Instance()->EventTag = eventTag; // so if it's a "second" event of the pair, keep the same event tag
  }
  
  //now in the same way, store the crucial info if the first was or not a candidate
  if(eventID % 2 == 0)  //this are the "first" events in the pairs
  {
    // so before entering the event, they are candidate by definition
    CreateTree::Instance()->FirstIsCoincidenceCandidate = 0;
  }
  else //for the second ones, they keep the FirstIsCoincidenceCandidate to 0 only if the one before was actually a candidate
  {
    CreateTree::Instance()->FirstIsCoincidenceCandidate = FirstIsCoincidenceCandidate;
  }
  //in this way, we will generate always the first gamma, but the second will be generated only if the first was a candidate

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void normalizationEventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //
  if(CreateTree::Instance()->totalEnergyDeposited > 0.5)
    CreateTree::Instance()->FirstIsCoincidenceCandidate = 0;
  else
    CreateTree::Instance()->FirstIsCoincidenceCandidate = 1;
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
