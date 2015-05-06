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
// $Id: normalizationActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file normalizationActionInitialization.cc
/// \brief Implementation of the normalizationActionInitialization class

#include "normalizationActionInitialization.hh"
#include "normalizationPrimaryGeneratorAction.hh"
#include "normalizationDetectorConstruction.hh"
#include "normalizationRunAction.hh"
#include "normalizationEventAction.hh"
#include "normalizationSteppingAction.hh"
#include "normalizationStackingAction.hh"
#include "normalizationSteppingVerbose.hh"
#include "ConfigFile.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

normalizationActionInitialization::normalizationActionInitialization(ConfigFile& config)
 : G4VUserActionInitialization(),
   fConfig(config)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

normalizationActionInitialization::~normalizationActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void normalizationActionInitialization::BuildForMaster() const
{
  SetUserAction(new normalizationRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void normalizationActionInitialization::Build() const
{
  
  
  SetUserAction(new normalizationPrimaryGeneratorAction(fConfig));
  normalizationRunAction* runAction = new normalizationRunAction();
  SetUserAction(runAction);
  normalizationEventAction* eventAction = new normalizationEventAction(runAction);
  SetUserAction(eventAction);
  SetUserAction(new normalizationSteppingAction(fConfig));
  SetUserAction(new normalizationStackingAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VSteppingVerbose*
               normalizationActionInitialization::InitializeSteppingVerbose() const
{
  return new normalizationSteppingVerbose();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
