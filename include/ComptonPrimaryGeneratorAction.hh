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
/// \file ComptonPrimaryGeneratorAction.hh
/// \brief Definition of the ComptonPrimaryGeneratorAction class

#ifndef ComptonPrimaryGeneratorAction_h
#define ComptonPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4Box;
class G4GenericMessenger;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

class ComptonPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ComptonPrimaryGeneratorAction();
    virtual ~ComptonPrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);

    void setPosX(G4double x){xpos=x;}
    void setPosY(G4double y){ypos=y;}
    void setPosZ(G4double z){zpos=z;}

    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

    G4double getGunXPosition() const {return xpos;}
    G4double getGunYPosition() const {return ypos;}
    G4double getGunZPosition() const {return zpos;}

  private:
    void DefineCommands();
    G4GenericMessenger *fMessenger;

    G4ParticleGun*  fParticleGun;
    G4double xpos, ypos, zpos;
    // pointer a to G4 gun class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
