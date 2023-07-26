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
/// \file ComptonEventAction.hh
/// \brief Definition of the ComptonEventAction class

#ifndef ComptonEventAction_h
#define ComptonEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class ComptonRunAction;
class G4GenericMessenger;


/// Event action class
///

class ComptonEventAction : public G4UserEventAction
{
  public:
    ComptonEventAction(ComptonRunAction* runAction);
    virtual ~ComptonEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdep(G4double edep) { fEdep += edep; }
    void Add3LayerCoincidence(G4int number) { f3LayerCoincidence += number; }
    void AddTotalNumber(G4int totalnumber) {fTotalNumber += totalnumber;}

    void Set2TotalEHigh(G4double val);
    void Set2TotalELow(G4double val);
    void Set2ScatterEHigh(G4double val);
    void Set2ScatterELow(G4double val);

    G4double AddFluction(G4double val);

    G4int GetTotalNumber() { return fTotalNumber; }

  private:
    void DefineCommands();
    G4GenericMessenger *fMessenger;

    ComptonRunAction* fRunAction;
    G4double     fEdep;
    G4int        f3LayerCoincidence;
    G4int        fTotalNumber;
    G4int fOrdered3, fOrdered2;
    G4int f2LayerCoincidence;

    G4int fLayer1ID;
    G4int fLayer2ID;
    G4int fLayer3ID;

    G4double f2TotalE_low, f2TotalE_high, f2ScatterE_low, f2ScatterE_high;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
