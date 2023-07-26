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
/// \file ComptonLayer1Hit.cc
/// \brief Implementation of the ComptonLayer1Hit class

#include "ComptonLayer1Hit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<ComptonLayer1Hit>* ComptonLayer1HitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonLayer1Hit::ComptonLayer1Hit()
: G4VHit(), 
   fEdep(0),
   fTime(0),
   fScatteringPos(0),
   fProcessName(""),
   fVolumeName(""),
   fParticleName(""),
   fCellID(-1),
   fParentID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//ComptonLayer1Hit::ComptonLayer1Hit(G4int layerID)
//: G4VHit(),
 // fLayerID(layerID), fTime(0.), fLocalPos(0), fWorldPos(0), fProcessName("")
//{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonLayer1Hit::~ComptonLayer1Hit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonLayer1Hit::ComptonLayer1Hit(const ComptonLayer1Hit &right)
: G4VHit(),
  fEdep(right.fEdep),
  fTime(right.fTime),
  fScatteringPos(right.fScatteringPos),
  fProcessName(right.fProcessName),
  fVolumeName(right.fVolumeName),
  fParticleName(right.fParticleName),
  fCellID(right.fCellID),
  fParentID(right.fParentID)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const ComptonLayer1Hit& ComptonLayer1Hit::operator=(const ComptonLayer1Hit &right)
{
  fEdep = right.fEdep;
  fTime = right.fTime;
  fScatteringPos = right.fScatteringPos;
  fProcessName = right.fProcessName;
  fVolumeName = right.fVolumeName;
  fParticleName = right.fParticleName;
  fCellID = right.fCellID;
  fParentID = right.fParentID;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ComptonLayer1Hit::operator==(const ComptonLayer1Hit &/*right*/) const
{
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonLayer1Hit::Draw()
{
  auto visManager = G4VVisManager::GetConcreteInstance();
  if (! visManager) return;

  G4Circle circle(fScatteringPos);
  circle.SetScreenSize(2);
  circle.SetFillStyle(G4Circle::filled);
  G4Colour colour(0.4,0.4,0.);
  G4VisAttributes attribs(colour);
  circle.SetVisAttributes(attribs);
  visManager->Draw(circle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* ComptonLayer1Hit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("ComptonLayer1Hit",isNew);

  if (isNew) {
      (*store)["HitType"] 
        = G4AttDef("HitType","Hit Type","Physics","","G4String");
      
      (*store)["ID"]
        = G4AttDef("ID","ID","Physics","","G4int");

      (*store)["PID"]
        = G4AttDef("PID","ParentID","Physics","","G4int");

      (*store)["EnergyDeposition"]
        = G4AttDef("Edep","EnergyDeposition","Physics","G4BestUnit","G4double");
      
      (*store)["Time"] 
        = G4AttDef("Time","Time","Physics","G4BestUnit","G4double");
      
      (*store)["Pos"] 
        = G4AttDef("Pos", "Position", "Physics","G4BestUnit","G4ThreeVector");

      (*store)["Process"]
              =G4AttDef("Process", "ProcessName", "Physics" ,"", "G4String");

      (*store)["Particle"]
              =G4AttDef("Particle", "ParitcleName", "Physics" ,"", "G4String");

      (*store)["Vol"]
              =G4AttDef("Vol", "VolumeName", "Physics", "", "G4String");
  }
  
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* ComptonLayer1Hit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;
  
  values
    ->push_back(G4AttValue("HitType","ComptonLayer1Hit",""));
  values
    ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fCellID),""));
  values
    ->push_back(G4AttValue("PID",G4UIcommand::ConvertToString(fParentID),""));
  values
    ->push_back(G4AttValue("Edep",G4BestUnit(fEdep,"Energy"),""));
  values
    ->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
  values
    ->push_back(G4AttValue("Pos",G4BestUnit(fScatteringPos,"Length"),""));
  values
    ->push_back(G4AttValue("Process", fProcessName,""));
  values
    ->push_back(G4AttValue("Particle", fParticleName,""));
  values
    ->push_back(G4AttValue("Vol", fVolumeName, ""));
  
  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonLayer1Hit::Print()
{
  G4cout << "] : "
  << " (nsec) --- local (x,y) " << fScatteringPos.x()
  << ", " << fScatteringPos.y() << G4endl
  << "process name " << fProcessName << G4endl
  << "Volume name " << fVolumeName << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
