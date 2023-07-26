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
/// \file ComptonRunAction.hh
/// \brief Definition of the ComptonRunAction class

#ifndef ComptonRunAction_h
#define ComptonRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

//txt format output
#include "iostream"

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class ComptonRunAction : public G4UserRunAction
{
  public:
    ComptonRunAction();
    virtual ~ComptonRunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    void AddEdep (G4double edep);
    void Add3LayerCoincidence (G4int number);
    void AddTotalNumber (G4int totalnumber);
    void AddOrdered3 (G4int number);
    void AddOrdered2(G4int numner);
    void Add2LayerCoincidence (G4int number);


    //txt format output
    void OutputFile(std::string outputFileName) {fOutput = outputFileName;}

  private:
    G4Accumulable<G4double> fEdep;
    G4Accumulable<G4double> fEdep2;
    G4Accumulable<G4int> f3LayerCoincidence;
    G4Accumulable<G4int> fTotalNumber;
    G4Accumulable<G4int> fOrdered3, fOrdered2;
    G4Accumulable<G4int> f2LayerCoincidence;

    //txt format output
    std::string fOutput;
};

#endif

