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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef normalizationPrimaryGeneratorAction_h
#define normalizationPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "ConfigFile.h"

class G4ParticleGun;
class G4Event;
class normalizationPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class normalizationPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    normalizationPrimaryGeneratorAction(ConfigFile& config);
    virtual ~normalizationPrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

//     void SetOptPhotonPolar();
//     void SetOptPhotonPolar(G4double);

  private:
    G4ParticleGun* fParticleGun;
    normalizationPrimaryGeneratorMessenger* fGunMessenger;
    G4double sourcex;  
    G4double sourcey;
    G4double sourcez;
    
    G4double preSourcex;
    G4double preSourcey;
    G4double preSourcez;
    
    G4double distance; 
    G4double energy;   
    G4double direction;
    
    G4double theta;
    G4double phi ;
    
    G4double firstTheta;
    G4double firstPhi;
    
    G4double esrThickness;
    G4double crystalx;
    G4int ncrystalx;
    
    G4bool isFirst;
    
    
    G4double phantomx;
    G4double phantomy;
    G4double phantomz;
  
    G4double posphantomx;
    G4double posphantomy;
    G4double posphantomz;
    
    
        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*normalizationPrimaryGeneratorAction_h*/
