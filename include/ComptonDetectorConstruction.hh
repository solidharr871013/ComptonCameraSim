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
/// \file ComptonDetectorConstruction.hh
/// \brief Definition of the ComptonDetectorConstruction class

#ifndef ComptonDetectorConstruction_h
#define ComptonDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;

class G4GenericMessenger;
/// Detector construction class to define materials and geometry.

class ComptonDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ComptonDetectorConstruction();
    virtual ~ComptonDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    void SetLayer1Thickness(G4double val);
    G4double GetLayer1Thickness() const { return fLayer1Thickness; }

    void SetLayer2Thickness(G4double val);
    G4double GetLayer2Thickness() const { return fLayer2Thickness; }
    
    void SetLayer3Thickness(G4double val);
    G4double GetLayer3Thickness() const {return fLayer3Thickness;}

    void SetLayer1ZPosition(G4double val);

    void SetLayer2ZPosition(G4double val);

    void SetLayer3ZPosition(G4double val);

    void SetCell1Size(G4double val);
    void SetCell2Size(G4double val);
    void SetCell3Size(G4double val);

    void Cell1Placement(G4int cellNO, G4double Celllen, G4String cellname, G4double gap,
                       G4LogicalVolume* celllog, G4LogicalVolume* layerlog);

    void Cell2Placement(G4int cellNO, G4double Celllen, G4String cellname, G4double gap,
                       G4LogicalVolume* celllog, G4LogicalVolume* layerlog);

    void Cell3Placement(G4int cellNO, G4double Celllen, G4String cellname, G4double gap,
                       G4LogicalVolume* celllog, G4LogicalVolume* layerlog);

    void Layer1Constructor(G4int cellNO, G4double Celllen,
                           G4double gap, G4double zpos,
                           G4LogicalVolume* layermotherlog);
    void Layer2Constructor(G4int cellNO, G4double Celllen,
                           G4double gap, G4double zpos,
                           G4LogicalVolume* layermotherlog);
    void Layer3Constructor(G4int cellNO, G4double Celllen,
                           G4double gap, G4double zpos,
                           G4LogicalVolume* layermotherlog);


    G4LogicalVolume* GetLayer1LogicalVolume() const { return fCell1Log; }
    G4LogicalVolume* GetLayer2LogicalVolume() const { return fCell2Log; }
    G4LogicalVolume* GetLayer3LogicalVolume() const {return fCell3Log;}

    void ConstructMaterial();

  private:
    void DefineCommands();
    G4GenericMessenger* fMessenger;

    G4Box* fSolidLayer1, *fSolidLayer2, *fSolidLayer3;

    G4Box* fSolidCell1, *fSolidCell2, *fSolidCell3;

    G4VPhysicalVolume* fLayer1Phys, *fLayer2Phys, *fLayer3Phys;

    std::vector<G4VPhysicalVolume*> fCell1Phys, fCell2Phys, fCell3Phys;

    G4LogicalVolume*  fLayer1Log, *fLayer2Log, *fLayer3Log;

    G4LogicalVolume* fCell1Log, *fCell2Log, *fCell3Log;

    G4LogicalVolume* fLogicalWorld;

    G4double fLayer1Thickness, fLayer2Thickness, fLayer3Thickness,
             fCell1Side, fCell2Side, fCell3Side;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

