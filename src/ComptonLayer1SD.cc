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
/// \file ComptonLayer1SD.cc
/// \brief Implementation of the ComptonLayer1 class

#include "ComptonLayer1SD.hh"
#include "ComptonLayer1Hit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "G4EventManager.hh"
#include "G4SteppingManager.hh"

#include "random"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonLayer1SD::ComptonLayer1SD(G4String name)
: G4VSensitiveDetector(name), 
  fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("LayerColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonLayer1SD::~ComptonLayer1SD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonLayer1SD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection 
    = new ComptonLayer1HitsCollection(SensitiveDetectorName,collectionName[0]);

  if (fHCID<0) { 
     fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);

 // G4cout << "the hit collection ID of Layer1 is " << fHCID << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ComptonLayer1SD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto track = step->GetTrack();
  
  auto preStepPoint = step->GetPreStepPoint();
  auto postStepPoint = step->GetPostStepPoint();

  auto touchable = preStepPoint->GetTouchable();
  //auto motherPhysical = touchable->GetVolume(1); // mother
  auto copyNo = touchable->GetVolume()->GetCopyNo();

  auto VolumeName = touchable->GetVolume()->GetLogicalVolume()->GetName();

  auto worldPos = postStepPoint->GetPosition();
  auto ScatteringPos
    = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

  auto proc = step->GetPostStepPoint()->GetProcessDefinedStep();
  if(!proc) { return  true; }

  auto procName = proc->GetProcessName();
  if(!procName){ return true;}

  auto particleName = track->GetParticleDefinition()->GetParticleName();

  G4double Edep = step->GetTotalEnergyDeposit();

  auto parentID = track->GetParentID();
  
  auto hit = new ComptonLayer1Hit;
  hit->SetCellID(copyNo);
  hit->SetScatteringPos(worldPos);
  hit->AddEdep(Edep);
  hit->SetTime(preStepPoint->GetGlobalTime());
  hit->SetProcessName(procName);
  hit->SetParticleName(particleName);
  hit->SetVolumeName(VolumeName);
  hit->SetParentID(parentID);
  
  fHitsCollection->insert(hit);
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
