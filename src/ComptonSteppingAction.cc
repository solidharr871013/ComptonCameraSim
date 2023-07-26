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
/// \file ComptonSteppingAction.cc
/// \brief Implementation of the ComptonSteppingAction class

#include "ComptonSteppingAction.hh"
#include "ComptonEventAction.hh"
#include "ComptonDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "random"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonSteppingAction::ComptonSteppingAction(ComptonEventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume1(nullptr),
  fScoringVolume2(nullptr),
  fScoringVolume3(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonSteppingAction::~ComptonSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonSteppingAction::UserSteppingAction(const G4Step* step)
{
  if ( !fScoringVolume1 || !fScoringVolume2 || !fScoringVolume3  ) {
   const ComptonDetectorConstruction* detectorConstruction
      = static_cast<const ComptonDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume1 = detectorConstruction->GetLayer1LogicalVolume();
    fScoringVolume2 = detectorConstruction->GetLayer2LogicalVolume();
    fScoringVolume3 = detectorConstruction->GetLayer3LogicalVolume();

  }


  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();



  // check if we are in scoring volume
  if ( volume != fScoringVolume1
       && volume != fScoringVolume2
       && volume != fScoringVolume3 ) return;

  // get information of the step
  G4StepPoint* preStep = step->GetPreStepPoint();
  G4StepPoint* postStep = step->GetPostStepPoint();

  G4Track* track = step->GetTrack();

  //G4LogicalVolume* PresentVolume = track->GetTouchable()->GetVolume()->GetLogicalVolume();
  G4String PresentVolumeName = track->GetTouchable()->GetVolume()->GetName();

 // G4LogicalVolume* NextVolume = track->GetNextTouchable()->GetVolume()->GetLogicalVolume();
  G4String NextVolumeName = track->GetNextTouchable()->GetVolume()->GetName();


  G4String PreStepVolumeName = preStep->GetTouchable()->GetVolume()->GetName();
  G4String PostStepVolumeName = postStep->GetTouchable()->GetVolume()->GetName();

 /* if using the prestepPoint without this code, it will get Segmentation Fault(because it
  might point to a empty pointer).  */
  G4String ProcessName = "";
  if(postStep!= nullptr)
  {
    const G4VProcess* proc = postStep->GetProcessDefinedStep();
    if(proc != nullptr)
    {
      ProcessName = proc->GetProcessName();
    }
  }
  //if(track->GetParentID()>0){
    //  G4String createProcess = track->GetCreatorProcess()->GetProcessName();
      //G4cout << "Creating process " << createProcess << G4endl;
  //}



  if(volume == fScoringVolume2){
  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();

  fEventAction->AddEdep(edepStep);
 // G4cout << "the energy deposites in " << volume->GetName() << G4endl;

}

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

