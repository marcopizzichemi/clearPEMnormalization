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

#include "normalizationPrimaryGeneratorAction.hh"
// #include "normalizationPrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "CreateTree.hh"

#include <iostream>
#include <fstream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

normalizationPrimaryGeneratorAction::normalizationPrimaryGeneratorAction(ConfigFile& config)
 : G4VUserPrimaryGeneratorAction(), 
   fParticleGun(0)
//    ,fConfig(config)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  //create a messenger for this class
  //fGunMessenger = new normalizationPrimaryGeneratorMessenger(this);

  //default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0*ns);
  
  //source 
  //x position
//   sourcex = config.read<double>("sourcex");
//   G4cout << "Gamma source x position [mm]: " << sourcex << G4endl;
//   //y position of source
//   sourcey = config.read<double>("sourcey");
//   G4cout << "Gamma source y position [mm]: " << sourcey << G4endl;
//   //x position of source
//   distance = config.read<double>("distance");
//   G4cout << "Distance of source from back ESR [mm]: " << distance << G4endl;
//   
//   //retrieve dimensions of the other elements and calculate the position of the source
//   G4double crystalz = config.read<double>("crystalz");
//   G4double greaseBack = config.read<double>("greaseBack");
//   G4double glassBack = config.read<double>("glassBack");
//   G4double airBack = config.read<double>("airBack");
//   //multiply for units [mm]
//   crystalz   = crystalz*mm;
//   greaseBack = greaseBack*mm;
//   glassBack  = glassBack*mm;
//   airBack    = airBack*mm;
//   G4double fakeAir = 0.1*mm;  //fixed to 0.1, it has no physical meaning in our simulation
//   G4double sourcez = -(distance + (crystalz/2.0) + greaseBack + glassBack + airBack + fakeAir);
//   G4cout << "Gamma source z position [mm]: " << sourcez << G4endl;
//   //energy of gammas
//   energy = config.read<double>("energy");
//   G4cout << "Energy of gamma source [KeV]: " << energy << G4endl;
//   //gamma direction
//   direction = config.read<int>("direction");
//   if(direction == 0)
//   {
//     G4cout << "Gammas shoot parallel to the z axis" << G4endl;
//   }
//   else
//   {
//     G4cout << "Gammas shoot randomly towards the matrix" << G4endl;
//   }
//   //find the x and y of matrix
//   crystalx = config.read<double>("crystalx");
//   //G4double crystaly = config.read<double>("crystaly");
//   ncrystalx = config.read<int>("ncrystalx");
//   //G4int ncrystaly = config.read<int>("ncrystaly");
//   esrThickness = config.read<double>("esrThickness");
  
  phantomx = config.read<double>("phantomx");
  phantomy = config.read<double>("phantomy");
  phantomz = config.read<double>("phantomz");
  
//   G4cout << "Phantom x length [mm]: " << phantomx << G4endl;
//   G4cout << "Phantom y length [mm]: " << phantomy << G4endl;
//   G4cout << "Phantom z length [mm]: " << phantomz << G4endl;
  
  posphantomx = config.read<double>("posphantomx");
  posphantomy = config.read<double>("posphantomy");
  posphantomz = config.read<double>("posphantomz");

//   G4cout << "Phantom x pos [mm]: " << posphantomx << G4endl;
//   G4cout << "Phantom y pos [mm]: " << posphantomy << G4endl;
//   G4cout << "Phantom z pos [mm]: " << posphantomz << G4endl;

  
  fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0));
  
  theta = 0;
  phi = 0;
    
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta) ));
  fParticleGun->SetParticleEnergy(511*keV);
  isFirst = true;
  firstTheta = theta;
  firstPhi = phi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

normalizationPrimaryGeneratorAction::~normalizationPrimaryGeneratorAction()
{
  delete fParticleGun;
//   delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void normalizationPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
     //cylindrical source constants
     double alpha = G4UniformRand()*2*CLHEP::pi;
     double ranradius = sqrt(G4UniformRand());
     double radius = ranradius*104.54;
  
  
//   G4double halfXatrix = ((crystalx + esrThickness) * ncrystalx) / 2.0;
//   G4double angleLimit = atan(halfXatrix / distance);
  
//   theta = (G4UniformRand() * 2.0*angleLimit) - angleLimit;
//   phi = G4UniformRand() * CLHEP::pi;
  
  if(isFirst)
  {
    //G4cout << "First" << G4endl;
    //theta = G4UniformRand() * 2.0 * CLHEP::pi;
   
    phi = G4UniformRand() * CLHEP::pi;
    double nrand = G4UniformRand()*2.0-1.0;
    theta = acos (nrand);
    
     
    //cylindrical source rotated
    sourcey = radius*sin(alpha);
    sourcez = radius*cos(alpha);
    sourcex = (G4UniformRand() * phantomx) - ( phantomx/2.0) + posphantomx;
        
    //planar source
     //sourcex = (G4UniformRand() * phantomx) - ( phantomx/2.0) + posphantomx;
     //sourcey = (G4UniformRand() * phantomy) - ( phantomy/2.0) + posphantomy;
     //sourcez = (G4UniformRand() * phantomz) - ( phantomz/2.0) + posphantomz;
       
    //G4cout << "sourcex = " << sourcex << G4endl;
    //G4cout << "sourcey = " << sourcey << G4endl;
    //G4cout << "sourcez = " << sourcez << G4endl;
    //G4cout << "Theta = " << theta << G4endl;
    //G4cout << "Phi = " << phi << G4endl;
    
    //octagonal source
    //sourcex = (G4UniformRand() * phantomx) - ( phantomx/2.0) + posphantomx;
    //sourcey = INFINITY;
    //double distanza = 209.08;
    //while( (sourcey > (sourcez-147.82)) )
    //{
    //sourcez = (G4UniformRand() * (distanza)) - (distanza/2.0);
    //sourcey = (G4UniformRand() * (distanza)) - (distanza/2.0);
    //}
    
    preSourcex = sourcex;
    preSourcey = sourcey;
    preSourcez = sourcez;
    
    firstTheta = theta ;
    firstPhi = phi;
    isFirst = false;
    }
    else
    {
//     G4cout << "Second" << G4endl;
    theta = firstTheta + CLHEP::pi;
    phi = firstPhi;
    
    sourcex = preSourcex;
    sourcey = preSourcey;
    sourcez = preSourcez;
    
//     G4cout << "Theta = " << theta << G4endl;
//     G4cout << "Phi = " << phi << G4endl;
    
    isFirst = true;
  }
  
  //G4cout << "Theta = " << theta << G4endl;
  //G4cout << "Phi = " << phi << G4endl;
  
//   if(direction == 0)
//   {
//     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
//   }
//   else
//   {
  fParticleGun->SetParticlePosition(G4ThreeVector(sourcex,sourcey,sourcez));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta) ));
  
  
//   }
  if(CreateTree::Instance()->DoEvent == 0) //generate the gamma only if the flag is 0, which will always be the case for First gammas, while for 
    //Second gammas depends on the energy deposited by the First (it makes no sense to simulate the second gamma if anyway this pair would not do a coincidence)
  {
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  
  ofstream myfile;
  myfile.open ("SourcePosition.txt",ios::app);
  myfile << sourcex << " " << sourcey << " " << sourcez << " " << firstTheta << " " << firstPhi << std::endl;
  myfile.close();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void normalizationPrimaryGeneratorAction::SetOptPhotonPolar()
// {
//  G4double angle = G4UniformRand() * 360.0* deg;
//  SetOptPhotonPolar(angle);
// }
// 
// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// 
// void normalizationPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
// {
//  if (fParticleGun->GetParticleDefinition()->GetParticleName()!="opticalphoton")
//    {
//      G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
//                "the particleGun is not an opticalphoton" << G4endl;
//      return;
//    }
// 
//  G4ThreeVector normal (1., 0., 0.);
//  G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
//  G4ThreeVector product = normal.cross(kphoton);
//  G4double modul2       = product*product;
//  
//  G4ThreeVector e_perpend (0., 0., 1.);
//  if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
//  G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
//  
//  G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
//  fParticleGun->SetParticlePolarization(polar);
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
