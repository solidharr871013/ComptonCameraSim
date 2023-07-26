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
/// \file ComptonRunAction.cc
/// \brief Implementation of the ComptonRunAction class

#include "ComptonRunAction.hh"
#include "ComptonPrimaryGeneratorAction.hh"
#include "ComptonDetectorConstruction.hh"
#include "ComptonEventAction.hh"
#include "ComptonAnalysis.hh"
// #include "ComptonRun.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonRunAction::ComptonRunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.),
  f3LayerCoincidence(0),
  fTotalNumber(0),
  fOrdered3(0),
  fOrdered2(0),
  f2LayerCoincidence(0)
{ 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2); 
  accumulableManager->RegisterAccumulable(f3LayerCoincidence);
  accumulableManager->RegisterAccumulable(fTotalNumber);
  accumulableManager->RegisterAccumulable(fOrdered3);
  accumulableManager->RegisterAccumulable(fOrdered2);
  accumulableManager->RegisterAccumulable(f2LayerCoincidence);

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->SetVerboseLevel(1);
  analysis->SetFileName("ComptonHisto");
  analysis->CreateH1("3layerEdep1","edep1",2000,0,2*MeV);//ID=0
  analysis->CreateH1("3layerEdep2","edep2",2000,0,2*MeV); //ID=1
  analysis->CreateH1("3layerEdep3","edep3",2000,0,2*MeV); //ID=2
  analysis->CreateH1("2layerEdepScat","edepScat",2000,0,2*MeV); //ID=3
  analysis->CreateH1("2layerEdepAbs","edepAbs",2000,0,2*MeV); //ID=4

  /**************************
  analysis->CreateH2("scatVsAbs","ScatVsAbs",2000,0,2*MeV,
                                             2000,0,2*MeV);//id=0
  analysis->CreateH2("Lay1VsLay2","1Vs2",2000,0,2*MeV,
                                            2000,0,2*MeV);//id=1
  analysis->CreateH2("Lay2VsLay3","2Vs3",2000,0,2*MeV,
                                            2000,0,2*MeV);//id=2
  ********************************/

  /********************** 3layer coincidence, id=0*************/
  analysis->CreateNtuple("3layerEventRecorder","Event Info");

  analysis->CreateNtupleIColumn("Event_No.");//0

  analysis->CreateNtupleSColumn("1_particle");//1
  analysis->CreateNtupleDColumn("1_eDep");//
  analysis->CreateNtupleDColumn("1_xPos");//
  analysis->CreateNtupleDColumn("1_ypos");//
  analysis->CreateNtupleDColumn("1_zPos");//
  analysis->CreateNtupleDColumn("1_time");//
  analysis->CreateNtupleSColumn("1_process");//

  analysis->CreateNtupleSColumn("2_particle");//8
  analysis->CreateNtupleDColumn("2_eDep");//
  analysis->CreateNtupleDColumn("2_xPos");//
  analysis->CreateNtupleDColumn("2_ypos");//
  analysis->CreateNtupleDColumn("2_zPos");//
  analysis->CreateNtupleDColumn("2_time");//
  analysis->CreateNtupleSColumn("2_process");//

  analysis->CreateNtupleSColumn("3_particle");//15
  analysis->CreateNtupleDColumn("3_eDep");//
  analysis->CreateNtupleDColumn("3_xPos");//
  analysis->CreateNtupleDColumn("3_ypos");//
  analysis->CreateNtupleDColumn("3_zPos");//
  analysis->CreateNtupleDColumn("3_time");//
  analysis->CreateNtupleSColumn("3_process");//
  analysis->FinishNtuple();
  
  /********************2layer coincidence, id =1*******************/
  analysis->CreateNtuple("2layerEventRecorder","Event Info");

  analysis->CreateNtupleIColumn("Event No.");//0

  analysis->CreateNtupleSColumn("Scat_particle");//1
  analysis->CreateNtupleDColumn("Scat_eDep");//
  analysis->CreateNtupleDColumn("Scat_xPos");//
  analysis->CreateNtupleDColumn("Scat_ypos");//
  analysis->CreateNtupleDColumn("Scat_zPos");//
  analysis->CreateNtupleDColumn("Scat_time");//
  analysis->CreateNtupleSColumn("Scat_process");//

  analysis->CreateNtupleSColumn("Abs_particle");//8
  analysis->CreateNtupleDColumn("Abs_eDep");//
  analysis->CreateNtupleDColumn("Abs_xPos");//
  analysis->CreateNtupleDColumn("Abs_ypos");//
  analysis->CreateNtupleDColumn("Abs_zPos");//
  analysis->CreateNtupleDColumn("Abs_time");//
  analysis->CreateNtupleSColumn("Abs_process");//
  analysis->FinishNtuple();
  
  /************************* run recorder, id=2***********/
  
  analysis->CreateNtuple("run_recorder","run info");
  analysis->CreateNtupleDColumn("2layerCoincidence_number");
  analysis->CreateNtupleDColumn("ordered_2layerCoincidence");
  analysis->CreateNtupleDColumn("3layerCoincidence_number");
  analysis->CreateNtupleDColumn("ordered_3layerCoincidence");
  analysis->CreateNtupleDColumn("totalEventNumber");
  analysis->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonRunAction::~ComptonRunAction()
{delete G4AnalysisManager::Instance();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonRunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  G4double edepTest = fEdep.GetValue();
  G4cout << "the test energy value is " << edepTest/keV << G4endl;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  // cast int to double type
  G4double NO3LayerCoincidence = f3LayerCoincidence.GetValue();
  G4double NOTotalEvent = fTotalNumber.GetValue();
  G4double NoOrdered3 = fOrdered3.GetValue();
  G4double NoOrdered2 = fOrdered2.GetValue();
  G4double No2LayerCoincidence = f2LayerCoincidence.GetValue();

  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const ComptonDetectorConstruction* detectorConstruction
   = static_cast<const ComptonDetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetLayer2LogicalVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  G4double Compton3EventRatio = NO3LayerCoincidence/NOTotalEvent,
           Compton2EventRatio = No2LayerCoincidence/NOTotalEvent;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const ComptonPrimaryGeneratorAction* generatorAction
   = static_cast<const ComptonPrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  G4double thicknessLayer1 = detectorConstruction->GetLayer1Thickness();


  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The thickness of layer1 is " << G4BestUnit(thicknessLayer1, "Length") << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : " 
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "The number of total event is " << NOTotalEvent << G4endl
     << "The number of 3LayerCoincidence event is " << NO3LayerCoincidence << G4endl
     << "The number of ordered 3LayerCoincident event is " << NoOrdered3 << G4endl
     << "The number of 2LayerCoincidence event is " << No2LayerCoincidence << G4endl
     << "The number of ordered 2LayerCoincident event is " << NoOrdered2 << G4endl
     << "The 3 compton event ratio is " << std::setw(3) << Compton3EventRatio << G4endl
     << "The 2 compton event ratio is " << std::setw(3) << Compton2EventRatio << G4endl
     << "The total compton event ratio is " << std::setw(3)
     << Compton2EventRatio+Compton3EventRatio
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  
  analysis->FillNtupleDColumn(2,0,No2LayerCoincidence);
  analysis->FillNtupleDColumn(2,1,NoOrdered2);
  analysis->FillNtupleDColumn(2,2,NO3LayerCoincidence);
  analysis->FillNtupleDColumn(2,3,NoOrdered3);
  analysis->FillNtupleDColumn(2,4,NOTotalEvent);
  analysis->AddNtupleRow(2);
  
  analysis->Write();
  analysis->CloseFile();

/****************Start txt format output***************

  std::ofstream SingleComptonEfficiency;


    SingleComptonEfficiency.open(fOutput + "SingleComptonEfficiency.txt", std::ios_base::app);
    if (SingleComptonEfficiency.is_open()){

      SingleComptonEfficiency << "The thickness is " << G4BestUnit(thicknessLayer1, "Length") << "; "
                              << "the efficiency is" << std::setw(4) << ComptonEventRatio << ";"
                              << std::endl;
      }
    else SingleComptonEfficiency << "Unable to open file" << std::endl;

****************End txt format output*********************/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonRunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

void ComptonRunAction::Add2LayerCoincidence(G4int number){
    f2LayerCoincidence += number;
}

void ComptonRunAction::Add3LayerCoincidence(G4int number){
    f3LayerCoincidence += number;
}

void ComptonRunAction::AddOrdered2(G4int number){
    fOrdered2 += number;
}

void ComptonRunAction::AddOrdered3(G4int number){
    fOrdered3 += number;
}

void ComptonRunAction::AddTotalNumber(G4int totalnumber){
    fTotalNumber += totalnumber;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

