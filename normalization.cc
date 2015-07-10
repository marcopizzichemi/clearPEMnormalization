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
//
//     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "normalizationPhysicsList.hh"

#include "normalizationDetectorConstruction.hh"
#include "normalizationActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

// #include "ConfigFile.hh"
#include "ConfigFile.h"
#include "CreateTree.hh"

#include <TNamed.h>
#include <TString.h>
#include <time.h>

#include "TRandom1.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " normalization [-m macro ] [-u UIsession] [-t nThreads] [-r seed] "
           << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

long int CreateSeed();


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

  if(argc<2) {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " normalization <config_file> [-m macro ] [-r seed] [-o output]" << G4endl;
    G4cerr << "   note: -m -r -o options are optional" << G4endl;
    exit(0);
  }
  
  //----------------------------------------------------------//
  //  Parse the config file                                   //
  //----------------------------------------------------------//
  
  std::cout<<"\n"<<std::endl;
  std::cout<<"###########################################################"<<std::endl;  
  std::cout<<"#                                                         #"<<std::endl;
  std::cout<<"#           Clear PEM normalization simulation            #"<<std::endl;  
  std::cout<<"#                                                         #"<<std::endl;  
  std::cout<<"###########################################################"<<std::endl;
  std::cout<<"\n\n"<<std::endl;
  std::cout<<"=====>   C O N F I G U R A T I O N   <====\n"<<std::endl;
   
  //config file
  G4cout << "Configuration file: '" << argv[1] << "'"<< G4endl;
  ConfigFile config(argv[1]);
  
  //output root file
  std::string filename;
  config.readInto(filename,"output");
  //G4cout << "Writing data to file '" << filename << ".root '..."<< G4endl;

  //run macro
  G4String macro;
  config.readInto(macro,"macro");

  //seed for random generator
  G4long myseed = config.read<long int>("seed");
  
  //type of session
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  
  //if passed by command line, override the filename, macro and myseed (the others are not important for us)
  for (int i = 0; i < argc ; i++)
  {
    if(std::string(argv[i]) == std::string("-m")) macro   = argv[i+1];
    else if(std::string(argv[i]) == std::string("-r")) myseed  = atoi(argv[i+1]);
    else if(std::string(argv[i]) == std::string("-o")) filename = argv[i+1];
  }
  
  if(myseed==-1) 
  {
     G4cout<<"Creating random seed..."<<G4endl;   
     myseed=CreateSeed();
  }
  
  G4cout << "Reading macro '" << macro << "' ..." << G4endl;
  G4cout<< "Random seed : "<< myseed <<G4endl;
  G4cout << "Output file '" << filename << ".root '..."<< G4endl;
  
  
  G4double headx = config.read<double>("headx");
  G4double heady = config.read<double>("heady");
  G4double headz = config.read<double>("headz");
  
  G4cout << "Heads x length [mm]: " << headx << G4endl;
  G4cout << "Heads y length [mm]: " << heady << G4endl;
  G4cout << "Heads z length [mm]: " << headz << G4endl;

  G4double posdh0x = config.read<double>("posdh0x"); 
  G4double posdh0y = config.read<double>("posdh0y");
  G4double posdh0z = config.read<double>("posdh0z");
  
  G4cout << "DH0 x pos [mm]: " << posdh0x << G4endl;
  G4cout << "DH0 y pos [mm]: " << posdh0y << G4endl;
  G4cout << "DH0 z pos [mm]: " << posdh0z << G4endl;
  
  G4double posdh1x = config.read<double>("posdh1x");
  G4double posdh1y = config.read<double>("posdh1y");
  G4double posdh1z = config.read<double>("posdh1z");
  
  G4cout << "DH1 x pos [mm]: " << posdh1x << G4endl;
  G4cout << "DH1 y pos [mm]: " << posdh1y << G4endl;
  G4cout << "DH1 z pos [mm]: " << posdh1z << G4endl;

  G4double phantomx = config.read<double>("phantomx");
  G4double phantomy = config.read<double>("phantomy");
  G4double phantomz = config.read<double>("phantomz");
  
  G4cout << "Phantom x length [mm]: " << phantomx << G4endl;
  G4cout << "Phantom y length [mm]: " << phantomy << G4endl;
  G4cout << "Phantom z length [mm]: " << phantomz << G4endl;

  G4double posphantomx = config.read<double>("posphantomx");
  G4double posphantomy = config.read<double>("posphantomy");
  G4double posphantomz = config.read<double>("posphantomz");

  G4cout << "Phantom x pos [mm]: " << posphantomx << G4endl;
  G4cout << "Phantom y pos [mm]: " << posphantomy << G4endl;
  G4cout << "Phantom z pos [mm]: " << posphantomz << G4endl;
  
  //read crystal dimensions
//   G4double crystalx = config.read<double>("crystalx");
//   G4double crystaly = config.read<double>("crystaly");
//   G4double crystalz = config.read<double>("crystalz");
//   G4cout << "Crystal x length [mm]: " << crystalx << G4endl;
//   G4cout << "Crystal y length [mm]: " << crystaly << G4endl;
//   G4cout << "Crystal z length [mm]: " << crystalz << G4endl;
//   
//   G4int ncrystalx = config.read<int>("ncrystalx");
//   G4int ncrystaly = config.read<int>("ncrystaly");
//   G4cout << "Number crystal in x direction: " << ncrystalx << G4endl;
//   G4cout << "Number crystal in y direction: " << ncrystaly << G4endl;
//   
//   //esr thickness
//   G4double esrThickness = config.read<double>("esrThickness");
//   G4cout << "ESR thickness [mm]: " << esrThickness << G4endl;
// 
//   //esr lateral and back existance
//   G4bool lateralEsr = (bool) config.read<int>("lateralEsr");
//   G4bool backEsr = (bool) config.read<int>("backEsr");
//   if(lateralEsr)
//     G4cout << "Crystal(s) with lateral ESR " << G4endl;
//   else
//     G4cout << "No lateral ESR " << G4endl;
//   if(lateralEsr)
//     G4cout << "Crystal(s) with back ESR " << G4endl;
//   else
//     G4cout << "No back ESR " << G4endl;
//   
//   //grease front from matrix to glass
//   G4double greaseFront1 = config.read<double>("greaseFront1");
//   G4cout << "Grease layer between matrix and front glass [mm]: " << greaseFront1 << G4endl;
//   //glass front thickness
//   G4double glassFront = config.read<double>("glassFront");
//   G4cout << "Front glass light guide thickness [mm]: " << glassFront << G4endl;
//   //grease front from glass to mppc
//   G4double greaseFront2 = config.read<double>("greaseFront2");
//   G4cout << "Grease layer between front glass and mppc [mm]: " << greaseFront2 << G4endl;
//   //epoxy
//   G4double epoxy = config.read<double>("epoxy");
//   G4cout << "Epoxy layer of mppc [mm]: " << epoxy << G4endl;
//   //mppc dimensions
//   G4double mppcx = config.read<double>("mppcx");
//   G4double mppcy = config.read<double>("mppcy");
//   G4double mppcz = config.read<double>("mppcz");
//   G4double mppcGap = config.read<double>("mppcGap");
//   G4cout << "Mppc x dimension [mm]: " << mppcx << G4endl;
//   G4cout << "Mppc y dimension [mm]: " << mppcy << G4endl;
//   G4cout << "Mppc z dimension [mm]: " << mppcz << G4endl;
//   G4cout << "Gap between MPPC detectors [mm]: " << mppcGap << G4endl;
//   G4int nmppcx = config.read<int>("nmppcx");
//   G4int nmppcy = config.read<int>("nmppcy");
//   G4cout << "Number mppc in x direction: " << nmppcx << G4endl;
//   G4cout << "Number mppcs in y direction: " << nmppcy << G4endl;
//   
//   //grease back from matrix to glass
//   G4double greaseBack = config.read<double>("greaseBack");
//   G4cout << "Grease layer between matrix and back glass [mm]: " << greaseBack << G4endl;
//   //glass back
//   G4double glassBack = config.read<double>("glassBack");
//   G4cout << "Back glass light guide thickness [mm]: " << glassBack << G4endl;
//   //grease back from matrix to glass
//   G4double airBack = config.read<double>("airBack");
//   G4cout << "Air layer between back glass and esr [mm]: " << airBack << G4endl;
//   
//   //material etc
//   G4double lightyield = config.read<double>("lightyield");
//   G4cout << "Scintillator Light Yield [Ph/MeV]: " << lightyield << G4endl;
//   // energy resolution
//   G4double resolutionScale = config.read<double>("resolutionScale");
//   
//   //lateral depolishing
//   G4bool latdepolished = (bool) config.read<int>("latdepolished");
//   G4double latsurfaceroughness = config.read<double>("latsurfaceroughness");
//   G4double latsigmaalpha = config.read<double>("latsigmaalpha");
//   G4double latspecularlobe = config.read<double>("latspecularlobe");
//   G4double latspecularspike = config.read<double>("latspecularspike");
//   G4double latbackscattering = config.read<double>("latbackscattering");
//   G4cout << "Crystal surface finish: "; 
//   if(latdepolished)
//   {
//     G4cout << "Lateral depolishing:" << G4endl;
//     G4cout << "Sigma Alpha [radiants]: " << latsigmaalpha << G4endl;
//     G4cout << "Specular Lobe:  " << latspecularlobe << G4endl;
//     G4cout << "Specular Spike: " << latspecularspike << G4endl;
//     G4cout << "Back Scattering: " << latbackscattering << G4endl;
//     if(latsurfaceroughness != 0)
//       G4cout << "Surface roughness [nm]: " << latsurfaceroughness << G4endl; 
//   }
//   else
//   {
//     G4cout << " Lateral: perfect polished" << G4endl;
//   }
//   
//   //front and back depolishing (for real polishing state)
//   G4bool realdepolished = (bool) config.read<int>("realdepolished");
//   G4double realsurfaceroughness = config.read<double>("realsurfaceroughness");
//   G4double realsigmaalpha = config.read<double>("realsigmaalpha");
//   G4double realspecularlobe = config.read<double>("realspecularlobe");
//   G4double realspecularspike = config.read<double>("realspecularspike");
//   G4double realbackscattering = config.read<double>("realbackscattering");
//   G4cout << "Front and Back crystal surfaces finish: "; 
//   if(realdepolished)
//   {
//     G4cout << "Front and Back depolishing:" << G4endl;
//     G4cout << "Sigma Alpha [radiants]: " << realsigmaalpha << G4endl;
//     G4cout << "Specular Lobe:  " << realspecularlobe << G4endl;
//     G4cout << "Specular Spike: " << realspecularspike << G4endl;
//     G4cout << "Back Scattering: " << realbackscattering << G4endl;
//     if(realsurfaceroughness != 0)
//       G4cout << "Surface roughness [nm]: " << realsurfaceroughness << G4endl; 
//   }
//   else
//   {
//     G4cout << " Front and Back: perfect polished" << G4endl;
//   }
//   
//   
//   G4double esrTransmittance = config.read<double>("esrTransmittance");
//   if(esrTransmittance == -1)
//   {
//     G4cout << "Using esr transmittance values from measurement" << G4endl;
//   }
//   else
//   {
//     G4cout << "Esr transmittance fixed to " << esrTransmittance << " for all wavelenghts" << G4endl;
//   }
//   
//   
//   // quantum eff of the detectors
//   G4double quantumEff = config.read<double>("quantumEff");
//   
//   G4cout << "Resolution Scale: " << resolutionScale << G4endl; 
//   G4cout << "here" << G4endl;
//   G4cout << "Quantum efficiency: " << quantumEff << G4endl; 
//   
//   //distance of source from back of the module
//   G4double distance = config.read<double>("distance");
  
  
  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Seed the random number generator manually
  G4Random::setTheSeed(myseed);

  
  //create output ttree
//   CreateTree* mytree = new CreateTree("StandardTree",ncrystalx,ncrystaly,nmppcx,nmppcy,true,true);
  
  CreateTree* mytree = new CreateTree("StandardTree"/*,ncrystalx,ncrystaly,nmppcx,nmppcy*/);
  //((CreateTree*)mytree)->SetModuleElements(ncrystalx,ncrystaly,nmppcx,nmppcy);
  //mytree->SetModuleElements(ncrystalx,ncrystaly,nmppcx,nmppcy);
  
  //save the seed in the ttree
  CreateTree::Instance()->Seed = myseed;
  
  //create the detector construction
  normalizationDetectorConstruction* detector = new normalizationDetectorConstruction();
  
  //set the parameters of detector
  
  
  ((normalizationDetectorConstruction*)detector)->SetHeadsDimensions(headx,heady,headz);
  ((normalizationDetectorConstruction*)detector)->SetDH0Position(posdh0x,posdh0y,posdh0z);
  ((normalizationDetectorConstruction*)detector)->SetDH1Position(posdh1x,posdh1y,posdh1z);
  ((normalizationDetectorConstruction*)detector)->SetPhantomDimensions(phantomx,phantomy,phantomz);
  ((normalizationDetectorConstruction*)detector)->SetPhantomPosition(posphantomx,posphantomy,posphantomz);
  
  
//   ((normalizationDetectorConstruction*)detector)->SetCrystalDimensions(crystalx,crystaly,crystalz);
//   ((normalizationDetectorConstruction*)detector)->SetNumberOfCrystals(ncrystalx,ncrystaly);
//   ((normalizationDetectorConstruction*)detector)->SetThinAirThickness(esrThickness);
//   ((normalizationDetectorConstruction*)detector)->SetLateralEsr(lateralEsr);
//   ((normalizationDetectorConstruction*)detector)->SetBackEsr(backEsr);
//   ((normalizationDetectorConstruction*)detector)->SetGreaseFrontOne(greaseFront1);
//   ((normalizationDetectorConstruction*)detector)->SetGreaseFrontTwo(greaseFront2);
//   ((normalizationDetectorConstruction*)detector)->SetGlassFront(glassFront);
//   ((normalizationDetectorConstruction*)detector)->SetEpoxy(epoxy);
//   ((normalizationDetectorConstruction*)detector)->SetMppcDimensions(mppcx,mppcy,mppcz);
//   ((normalizationDetectorConstruction*)detector)->SetMppcGap(mppcGap);
//   ((normalizationDetectorConstruction*)detector)->SetNumberOfMPPC(nmppcx,nmppcy);
//   ((normalizationDetectorConstruction*)detector)->SetGreaseBack(greaseBack);
//   ((normalizationDetectorConstruction*)detector)->SetGlassBack(glassBack);
//   ((normalizationDetectorConstruction*)detector)->SetAirBack(airBack);
//   ((normalizationDetectorConstruction*)detector)->SetLightYield(lightyield);
//   
//   ((normalizationDetectorConstruction*)detector)->SetLateralDepolished(latdepolished);
//   ((normalizationDetectorConstruction*)detector)->SetLateralSurfaceRoughness(latsurfaceroughness);
//   ((normalizationDetectorConstruction*)detector)->SetLateralSurfaceSigmaAlpha(latsigmaalpha);
//   ((normalizationDetectorConstruction*)detector)->SetLateralSurfaceSpecularLobe(latspecularlobe);
//   ((normalizationDetectorConstruction*)detector)->SetLateralSurfaceSpecularSpike(latspecularspike);
//   ((normalizationDetectorConstruction*)detector)->SetLateralSurfaceBackScattering(latbackscattering);
//   
//   ((normalizationDetectorConstruction*)detector)->SetRealDepolished(realdepolished);
//   ((normalizationDetectorConstruction*)detector)->SetRealSurfaceRoughness(realsurfaceroughness);
//   ((normalizationDetectorConstruction*)detector)->SetRealSurfaceSigmaAlpha(realsigmaalpha);
//   ((normalizationDetectorConstruction*)detector)->SetRealSurfaceSpecularLobe(realspecularlobe);
//   ((normalizationDetectorConstruction*)detector)->SetRealSurfaceSpecularSpike(realspecularspike);
//   ((normalizationDetectorConstruction*)detector)->SetRealSurfaceBackScattering(realbackscattering);
//   
//   
//   ((normalizationDetectorConstruction*)detector)->SetEsrTransmittance(esrTransmittance);
//   ((normalizationDetectorConstruction*)detector)->SetSourceDistance(distance);
//   ((normalizationDetectorConstruction*)detector)->SetResolutionScale(resolutionScale);
  
  // Set mandatory initialization classes
  //
  // Detector construction
  runManager-> SetUserInitialization(detector);
  // Physics list
  runManager-> SetUserInitialization(new normalizationPhysicsList);
  // User action initialization
  runManager->SetUserInitialization(new normalizationActionInitialization(config));
  // Initialize G4 kernel
  //
  runManager->Initialize();

#ifdef G4VIS_USE
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer(); 
   
  if ( macro.size() ) {
     // Batch mode
     G4String command = "/control/execute ";
     UImanager->ApplyCommand(command+macro);
  }
  else // Define UI session for interactive mode
  {
#ifdef G4UI_USE
     G4UIExecutive * ui = new G4UIExecutive(argc,argv,session);
#ifdef G4VIS_USE
     UImanager->ApplyCommand("/control/execute vis.mac");
#else
     UImanager->ApplyCommand("/control/execute normalization.in");
#endif
     if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");
     ui->SessionStart();
     delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  
  
  
#ifdef G4VIS_USE
  delete visManager;
    
#endif
  delete runManager;
  
  //output root file
  std::stringstream outfileName;
  outfileName << filename << ".root";

  TFile* outfile = new TFile(TString(outfileName.str().c_str()),"RECREATE");
  outfile->cd();
  G4cout<<"Writing tree to file "<< outfileName.str() <<" ..."<<G4endl;
  mytree->GetTree()->Write();
  outfile->Write();
  outfile->Close();
  
  
  return 0;
}


long int CreateSeed()
{
  TRandom1 rangen;
  long int ss = time(0);
  std::cout<<"Time : "<<ss<<std::endl;
  ss+=getpid();
  std::cout<<"PID  : "<<getpid()<<std::endl;
  //
  FILE * fp = fopen ("/proc/uptime", "r");
  int uptime,upsecs;
  if (fp != NULL)
    {
      char buf[BUFSIZ];
      int res;
      char *b = fgets (buf, BUFSIZ, fp);
      if (b == buf)
        {
          /* The following sscanf must use the C locale.  */
          setlocale (LC_NUMERIC, "C");
          res = sscanf (buf, "%i", &upsecs);
          setlocale (LC_NUMERIC, "");
          if (res == 1)
            uptime = (time_t) upsecs;
        }

      fclose (fp);
    }

  std::cout<<"Uptime: "<<upsecs<<std::endl;
  ss+=upsecs;
  //
  std::cout<<"Seed for srand: "<<ss<<std::endl;
  srand(ss);
  rangen.SetSeed(rand());
  long int seed = round(1000000*rangen.Uniform());
  return seed;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
