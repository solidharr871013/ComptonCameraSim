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
/// \file ComptonLayer1Hit.hh
/// \brief Definition of the ComptonLayer1Hit class

#ifndef ComptonLayer1Hit_h
#define ComptonLayer1Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Drift chamber hit
///
/// It records:
/// - the layer ID
/// - the particle time
/// - the particle local and global positions

class ComptonLayer1Hit : public G4VHit
{
  public:
    ComptonLayer1Hit();
    //ComptonLayer1Hit(G4int layerID);
    ComptonLayer1Hit(const ComptonLayer1Hit &right);
    virtual ~ComptonLayer1Hit();

    const ComptonLayer1Hit& operator=(const ComptonLayer1Hit &right);
    G4bool operator==(const ComptonLayer1Hit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();

    void SetCellID(G4int z) { fCellID = z; }
    G4int GetCellID() const { return fCellID; }

    void SetEdep(G4double edep) { fEdep = edep; }
    void AddEdep(G4double edep){fEdep += edep;}
    G4double GetEdep() const { return fEdep; }

    void SetScatteringPos(G4ThreeVector xyz) { fScatteringPos = xyz; }
    G4ThreeVector GetScatteringPos() const { return fScatteringPos; }

    void SetProcessName(G4String procName) { fProcessName = procName; }
    G4String GetProcessName() const { return fProcessName; }

    void SetParticleName(G4String pName) { fParticleName = pName; }
    G4String GetParticleName() const { return fParticleName; }

    void SetVolumeName(G4String volName) { fVolumeName = volName; }
    G4String GetVolumeName() const { return fVolumeName; }

    void SetTime(G4double time){fTime = time;}
    G4double GetTime() const {return fTime;}

    void SetParentID(G4int PID){fParentID = PID;}
    G4int GetParentID() const {return  fParentID;}
    
  private:
    G4double fEdep;
    G4double fTime;
    G4ThreeVector fScatteringPos;
    G4String fProcessName;
    G4String fVolumeName;
    G4String fParticleName;
    G4int fCellID;
    G4int fParentID;

};

using ComptonLayer1HitsCollection = G4THitsCollection<ComptonLayer1Hit>;

extern G4ThreadLocal G4Allocator<ComptonLayer1Hit>* ComptonLayer1HitAllocator;

inline void* ComptonLayer1Hit::operator new(size_t)
{
  if (!ComptonLayer1HitAllocator) {
       ComptonLayer1HitAllocator = new G4Allocator<ComptonLayer1Hit>;
  }
  return (void*)ComptonLayer1HitAllocator->MallocSingle();
}

inline void ComptonLayer1Hit::operator delete(void* aHit)
{
  ComptonLayer1HitAllocator->FreeSingle((ComptonLayer1Hit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
