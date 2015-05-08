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
  //G4cout << "//***********************************************//" << G4endl;
  G4cout << "\r";
  G4cout << "Event Number " << eventID;//		<< G4endl;
  //G4cout << "//***********************************************//" << G4endl;

  //store the info that would be destroied by the clear command, so you can use it immediately later
  Int_t run = CreateTree::Instance()->Run;
  //long int eventTag = CreateTree::Instance()->EventTag;
  long int seed = CreateTree::Instance()->Seed;
  //int FirstIsCoincidenceCandidate = CreateTree::Instance()->FirstIsCoincidenceCandidate;
  
  //clear the create tree
  //CreateTree::Instance()->Clear(); // <----------------
  
  //these remained stored in the same way
  CreateTree::Instance()->Seed = seed;
  CreateTree::Instance()->Run = run;
  //store the event number
  CreateTree::Instance()->Event = event->GetEventID();
  //store the event tag (pairs)
  if(eventID % 2 == 0)  //this are the "first" events in the pairs
  {
    CreateTree::Instance()->GammaParity = 0;
    CreateTree::Instance()->DoEvent = 0; //always do the simulation for the first events
  }
  else 
  {
    CreateTree::Instance()->GammaParity = 1;
    // leave the DoEvent flag set as it has been set by the end of previous event
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void normalizationEventAction::EndOfEventAction(const G4Event* event)
{
  //now at the end of each event check this
  // if it's a first gamma do nothing
  // if it's a second gamma, and both energies were above 0.5, then save the event in the root file //TODO implement skip second if first is not candidate?
  // otherwise do nothing and clear
  
  if(CreateTree::Instance()->GammaParity == 0) //it's first
  {
    if(CreateTree::Instance()->First_totalEnergyDeposited < 0.5) // so if the First gamma is not candidate
    {
      CreateTree::Instance()->DoEvent = 1; //tell to the primary generator not to generate the gamma in the next event
    }
  }
  else  //it's second
  {
    if( (CreateTree::Instance()->First_totalEnergyDeposited > 0.5) && (CreateTree::Instance()->Second_totalEnergyDeposited > 0.5) ) // both gammas are candidate
    {
      CreateTree::Instance()->Fill();
    }
    CreateTree::Instance()->Clear();
  }
  
  
//   if(CreateTree::Instance()->totalEnergyDeposited > 0.5)
//     CreateTree::Instance()->FirstIsCoincidenceCandidate = 0;
//   else
//     CreateTree::Instance()->FirstIsCoincidenceCandidate = 1;
//   CreateTree::Instance()->Fill();
  //CreateTree::Instance()->Clear();
  
  //G4cout << "Total Energy deposited in this event = "<< CreateTree::Instance()->totalEnergyDeposited << " MeV" << G4endl;
  
  
} 
