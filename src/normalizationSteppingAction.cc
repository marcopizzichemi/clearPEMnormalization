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
// $Id: normalizationSteppingAction.cc 71007 2013-06-09 16:14:59Z maire $
//
/// \file normalizationSteppingAction.cc
/// \brief Implementation of the normalizationSteppingAction class

#include "normalizationSteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
// #include <sstream>
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "normalizationEventAction.hh"
// #include "B4Analysis.hh"
#include "CreateTree.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


normalizationSteppingAction::normalizationSteppingAction(
                      ConfigFile& config)
: G4UserSteppingAction()
    //fEventAction(eventAction)
{ 
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;
  fOpticalPhotonsStopped = 0;
  fGammaStopped = 0;
  fElectronStopped = 0;
  
//   quantumEff = config.read<double>("quantumEff");
  
  
  //ch[] = {0};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

normalizationSteppingAction::~normalizationSteppingAction()
{ ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void normalizationSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int eventNumber = G4RunManager::GetRunManager()->
                                              GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber) {
     G4cout << " Number of Scintillation Photons in previous event: "
            << fScintillationCounter << G4endl;
     G4cout << " Number of Cerenkov Photons in previous event: "
            << fCerenkovCounter << G4endl;
     fEventNumber = eventNumber;
     fScintillationCounter = 0;
     fCerenkovCounter = 0;
  }
  
  
	
  G4Track* track = step->GetTrack();
  G4String ParticleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();  
  G4String materialName = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName();
  
  G4VPhysicalVolume* thePrePV  = step->GetPreStepPoint() ->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = step->GetPostStepPoint()->GetPhysicalVolume();
                     
  G4String PreMaterialName ;
  G4String PostMaterialName ;
  
  if(thePrePV) PreMaterialName = thePrePV->GetLogicalVolume()->GetMaterial()->GetName();
  if(thePostPV) PostMaterialName = thePostPV->GetLogicalVolume()->GetMaterial()->GetName();
			       
			       
  //getting information from the particles
  
  //for opticalphoton, we want to know where they stopped 
  if (ParticleName == "opticalphoton") //if it's an opticalphoton
  {
    
    //BEGIN of debug part (flood maps not symmetric issue)
//     G4OpBoundaryProcessStatus boundaryStatus=Undefined;
//     static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;
//     
//     //find the boundary process only once
//     if(!boundary)
//     {
//       G4ProcessManager* pm = track->GetDefinition()->GetProcessManager();
//       G4int nprocesses = pm->GetProcessListLength();
//       G4ProcessVector* pv = pm->GetProcessList();
//       G4int i;
//       for( i=0;i<nprocesses;i++)
//       {
// 	if((*pv)[i]->GetProcessName()=="OpBoundary")
// 	{
// 	  boundary = (G4OpBoundaryProcess*)(*pv)[i];
// 	  break;
// 	}
//       }
//     }
//     boundaryStatus = boundary->GetStatus();
//     
//     if((materialName == "AirThinLayer") && (boundaryStatus == Transmission))
//     {
//       G4ThreeVector TransmissionPosition = track->GetStep()->GetPostStepPoint()->GetPosition(); //get the position vector
//       CreateTree::Instance()->TransmissionX.push_back(TransmissionPosition.getX());
//       CreateTree::Instance()->TransmissionY.push_back(TransmissionPosition.getY());
//       CreateTree::Instance()->TransmissionZ.push_back(TransmissionPosition.getZ());   
//     }
//     
    //END of debug part (flood maps not symmetric issue)
    
    
    //CreateTree::Instance()->NumOptPhotonsAbsorbed++;
//     if(track->GetTrackStatus()==fStopAndKill) //if it just died here
//     {
//       if(materialName == "SilicioMPPC")
//       {
// 	G4String detectorName = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
// 	int numb;
// 	std::istringstream ( detectorName ) >> numb;
	
	//take into account quantum efficiency
// 	G4double rand = G4UniformRand();
	
// 	if(rand < quantumEff) //absorb the optical photon, save the relevant data
// 	{
// 	  CreateTree::Instance()->DetectorHit[numb]++;
	  
// 	  G4ThreeVector OnDetectorPosition = track->GetStep()->GetPreStepPoint()->GetPosition(); //get the position vector
// 	  G4double globalTime = track->GetGlobalTime();
	  
// 	  CreateTree::Instance()->PositionX.push_back(OnDetectorPosition.getX());
//           CreateTree::Instance()->PositionY.push_back(OnDetectorPosition.getY());
// 	  CreateTree::Instance()->GlobalTime.push_back(globalTime/CLHEP::ns);
// 	}

//       } 
//     }
  }
  else //there's only gammas and electrons. for them, we want to know where they left energy (opticalphoton don't really leave any)
  {
    G4double edep;
    edep = step->GetTotalEnergyDeposit()/CLHEP::MeV; //check if there was an energy deposition
    if(edep != 0) //if there was en dep, save the data
    {
      if(materialName == "LYSO") //if the gamma is interacting with the detector, it does a huge number of cherenkov
      {
	//add total energy deposited
	CreateTree::Instance()->totalEnergyDeposited += edep;
	//take crystal name
	G4String crystalName = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
	//convert it to int
	
// 	if(crystalName == "DH0")
// 	{
// 	  CreateTree::Instance()->DH0EnergyDeposited += edep;
// 	  //add position of deposited energy
// 	  G4ThreeVector positionVector = track->GetStep()->GetPreStepPoint()->GetPosition(); //get the position vector
// 	  //call the methods to fill the position vectors
// 	  CreateTree::Instance()->DH0PosXEnDep.push_back(positionVector.getX());
// 	  CreateTree::Instance()->DH0PosYEnDep.push_back(positionVector.getY());
// 	  CreateTree::Instance()->DH0PosZEnDep.push_back(positionVector.getZ());
// 	}
// 	else if (crystalName == "DH1")
// 	{
// 	  CreateTree::Instance()->DH1EnergyDeposited += edep;
// 	  //add position of deposited energy
// 	  G4ThreeVector positionVector = track->GetStep()->GetPreStepPoint()->GetPosition(); //get the position vector
// 	  //call the methods to fill the position vectors
// 	  CreateTree::Instance()->DH1PosXEnDep.push_back(positionVector.getX());
// 	  CreateTree::Instance()->DH1PosYEnDep.push_back(positionVector.getY());
// 	  CreateTree::Instance()->DH1PosZEnDep.push_back(positionVector.getZ());
// 	}
	
// 	if((crystalName == "DH0") | (crystalName == "DH1"))
// 	{
	CreateTree::Instance()->EnergyDeposited.push_back(edep);
	G4ThreeVector positionVector = track->GetStep()->GetPreStepPoint()->GetPosition(); //get the position vector
	//call the methods to fill the position vectors
	CreateTree::Instance()->PosXEnDep.push_back(positionVector.getX());
	CreateTree::Instance()->PosYEnDep.push_back(positionVector.getY());
	CreateTree::Instance()->PosZEnDep.push_back(positionVector.getZ());
// 	}
	
	
	int numb;
	std::istringstream ( crystalName ) >> numb;
	//add energy deposited in this step to this crystal
	CreateTree::Instance()->Crystal.push_back(numb);
	
      }
    }
  }
  
  if (ParticleName == "opticalphoton") return;
  const std::vector<const G4Track*>* secondaries =
                                            step->GetSecondaryInCurrentStep();
  if (secondaries->size()>0) {
     for(unsigned int i=0; i<secondaries->size(); ++i) {
        if (secondaries->at(i)->GetParentID()>0) {
           if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
               == G4OpticalPhoton::OpticalPhotonDefinition()){
              if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
               == "Scintillation")
	      {
		CreateTree::Instance()->NumOptPhotons++;
		fScintillationCounter++;
	      }
              if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
               == "Cerenkov")
	      {
		fCerenkovCounter++;
		CreateTree::Instance()->NumCherenkovPhotons++;
	      }
           }
        }
     }
  }
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
