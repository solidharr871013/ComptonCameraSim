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
/// \file ComptonDetectorConstruction.cc
/// \brief Implementation of the ComptonDetectorConstruction class

#include "ComptonDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4GenericMessenger.hh"

#include "ComptonLayer1SD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonDetectorConstruction::ComptonDetectorConstruction()
: G4VUserDetectorConstruction(),
  fMessenger(nullptr),
  fLayer1Log(nullptr),
  fLayer2Log(nullptr),
  fLayer3Log(nullptr),
  fLayer1Thickness(15.*mm),
  fLayer2Thickness(10.*mm),
  fLayer3Thickness(10.*mm),
  fCell1Side(2.*mm),
  fCell2Side(2.*mm),
  fCell3Side(2.*mm)
{
    DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonDetectorConstruction::~ComptonDetectorConstruction()
{
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ComptonDetectorConstruction::Construct()
{
   //
   // Construct the material to be used
   //
   ConstructMaterial();
   G4Material* world_mat = G4Material::GetMaterial("G4_AIR");
   G4Material* Siboard_mat = G4Material::GetMaterial("Siboard");
   //G4Material* cell_mat = G4Material::GetMaterial("Gd3Ca5Al5O12");
   //G4Material* cell_mat = G4Material::GetMaterial("Lu2SiO5");
   //G4Material* layer_mat = G4Material::GetMaterial("BaSO4");

  //
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 7.5*m;
  G4double world_sizeZ  = 7.5*m;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  fLogicalWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(nullptr,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      fLogicalWorld,            //its logical volume
                      "World",               //its name
                      nullptr,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  G4double gap = 0.02*mm;
  G4int layer1_cell_num = 16, layer2_cell_num = 16, layer3_cell_num =16;
  G4double layer1ZPosition = -18*mm, layer2ZPosition = 0., layer3ZPosition = 18.0*mm;

  Layer1Constructor(layer1_cell_num,fCell1Side,gap,layer1ZPosition,fLogicalWorld);
  Layer2Constructor(layer2_cell_num,fCell2Side,gap,layer2ZPosition,fLogicalWorld);
  Layer3Constructor(layer3_cell_num,fCell3Side,gap,layer3ZPosition,fLogicalWorld);

  Cell1Placement(layer1_cell_num,fCell1Side,"cell1",gap,fCell1Log,fLayer1Log);
  Cell2Placement(layer2_cell_num,fCell2Side,"cell2",gap,fCell2Log,fLayer2Log);
  Cell3Placement(layer3_cell_num,fCell3Side,"cell3",gap,fCell3Log,fLayer3Log);

  G4double SiBoard_side = 5.*cm, SiBoard_thick = 4.*mm;

  G4Box* SiboardSolid = new G4Box("Siboard", 0.5*SiBoard_side,
                                             0.5*SiBoard_side,
                                             0.5*SiBoard_thick);

  G4LogicalVolume* SiboardLog = new G4LogicalVolume(SiboardSolid,
                                                    Siboard_mat,
                                                    "Siboard");

  G4ThreeVector boardPos1 = G4ThreeVector(0,0,25*mm),
                boardPos2 = G4ThreeVector(0,0,11*mm),
                boardPos3 = G4ThreeVector(0,0,7*mm),
                boardPos4 = G4ThreeVector(0,0,-7*mm),
                boardPos5 = G4ThreeVector(0,0,-11*mm),
                boardPos6 = G4ThreeVector(0,0,-25*mm);


  new G4PVPlacement(nullptr,
                    boardPos1,
                    SiboardLog,
                    "board1",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);

  new G4PVPlacement(nullptr,
                    boardPos2,
                    SiboardLog,
                    "board1",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);

  new G4PVPlacement(nullptr,
                    boardPos3,
                    SiboardLog,
                    "board1",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);

  new G4PVPlacement(nullptr,
                    boardPos4,
                    SiboardLog,
                    "board1",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);

  new G4PVPlacement(nullptr,
                    boardPos5,
                    SiboardLog,
                    "board1",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);

  new G4PVPlacement(nullptr,
                    boardPos6,
                    SiboardLog,
                    "board1",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);




  fLayer1Log->SetVisAttributes(G4VisAttributes::GetInvisible());
  fLayer2Log->SetVisAttributes(G4VisAttributes::GetInvisible());
  fLayer3Log->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4Color yellow(1.0,1.0,0.0);
  G4VisAttributes* cellVisAttributes = new G4VisAttributes(yellow);
  fCell1Log->SetVisAttributes(cellVisAttributes);
  fCell2Log->SetVisAttributes(cellVisAttributes);
  fCell3Log->SetVisAttributes(cellVisAttributes);

  //
  //always return the physical World
  //
  return physWorld;
}

void ComptonDetectorConstruction::ConstructMaterial(){
    // Get NistManager
    G4NistManager* man = G4NistManager::Instance();

    man->FindOrBuildMaterial("G4_AIR");

    // GAGG material
    G4Element*  O = man->FindOrBuildElement( "O", true);
    G4Element*  S = man->FindOrBuildElement( "S", true);
    G4Element* Si = man->FindOrBuildElement("Si", true);
    G4Element* Lu = man->FindOrBuildElement("Lu", true);
    G4Element* Gd = man->FindOrBuildElement("Gd", true);
    G4Element* Ca = man->FindOrBuildElement("Ca", true);
    G4Element* Ba = man->FindOrBuildElement("Ba", true);
    G4Element* Al = man->FindOrBuildElement("Al", true);
    G4Element* Zn = man->FindOrBuildElement("Zn", true);
    G4Element* Te = man->FindOrBuildElement("Te", true);
    G4Element* Cd = man->FindOrBuildElement("Cd", true);

    G4Material* LSO = new G4Material("Lu2SiO5",7.1*g/cm3,3);
    LSO->AddElement(Lu,2);
    LSO->AddElement(Si,1);
    LSO->AddElement(O,5);

    G4Material* GAGG = new G4Material("Gd3CaAl4O12", 6.63*g/cm3,4);
    GAGG->AddElement(Gd,3);
    GAGG->AddElement(Ca,1);
    GAGG->AddElement(Al,4);
    GAGG->AddElement(O,12);

    G4Material* BaSO4 = new G4Material("BaSO4", 4.49*g/cm3,3);
    BaSO4->AddElement(Ba,1);
    BaSO4->AddElement(S,1);
    BaSO4->AddElement(O,4);

    G4Material* CZT = new G4Material("CZT", 5.78*g/cm3, 3);
    CZT->AddElement(Cd,9);
    CZT->AddElement(Zn,1);
    CZT->AddElement(Te,10);


    G4Material* Siboard = new G4Material("Siboard", 2.33*g/cm3,1);
    Siboard->AddElement(Si,1);

     man->FindOrBuildMaterial("G4_SODIUM_IODIDE");

}

void ComptonDetectorConstruction::ConstructSDandField(){

    auto sdManger = G4SDManager::GetSDMpointer();
    G4String SDname;

    auto Layer1 = new ComptonLayer1SD(SDname = "/Layer1");
    sdManger->AddNewDetector(Layer1);
    fCell1Log->SetSensitiveDetector(Layer1);

    auto Layer2 = new ComptonLayer1SD(SDname = "/Layer2");
    sdManger->AddNewDetector(Layer2);
    fCell2Log->SetSensitiveDetector(Layer2);

    auto Layer3 = new ComptonLayer1SD(SDname = "/Layer3");
    sdManger->AddNewDetector(Layer3);
    fCell3Log->SetSensitiveDetector(Layer3);


}

void ComptonDetectorConstruction::SetLayer1Thickness(G4double val){
    fLayer1Thickness = val;
    //fSolidLayer1->SetZHalfLength(0.5*fLayer1Thickness);
    //fSolidCell1->SetZHalfLength(0.5*fLayer1Thickness);

    G4RunManager::GetRunManager()->ReinitializeGeometry();
    // tell G4RunManager that we change the geometry
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


void ComptonDetectorConstruction::SetLayer2Thickness(G4double val){
    fLayer2Thickness = val;
    //fSolidLayer2->SetZHalfLength(0.5*fLayer2Thickness);
    //fSolidCell2->SetZHalfLength(0.5*fLayer2Thickness);

    G4RunManager::GetRunManager()->ReinitializeGeometry();
    // tell G4RunManager that we change the geometry    
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void ComptonDetectorConstruction::SetLayer3Thickness(G4double val){
    fLayer3Thickness = val;
    //fSolidLayer3->SetZHalfLength(0.5*fLayer3Thickness);
    //fSolidCell3->SetZHalfLength(0.5*fLayer3Thickness);

    G4RunManager::GetRunManager()->ReinitializeGeometry();
    // tell G4RunManager that we change the geometry
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

/************** Setting the Layer position ***************/
void ComptonDetectorConstruction::SetLayer1ZPosition(G4double val){

    fLayer1Phys->SetTranslation(G4ThreeVector(0,0,val));

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void ComptonDetectorConstruction::SetLayer2ZPosition(G4double val){
    fLayer2Phys->SetTranslation(G4ThreeVector(0,0,val));

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void ComptonDetectorConstruction::SetLayer3ZPosition(G4double val){
    fLayer3Phys->SetTranslation(G4ThreeVector(0,0,val));

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

/****************** Setting the Cell size **********/
void ComptonDetectorConstruction::SetCell1Size(G4double val){
    fCell1Side = val;
    G4int cellNO = 16;
    G4double gap = 0.*mm;

    //change the cell size
    fSolidCell1->SetXHalfLength(0.5*fCell1Side);
    fSolidCell1->SetYHalfLength(0.5*fCell1Side);

    //new layer size
    fSolidLayer1->SetXHalfLength(0.5*((cellNO+1)*gap+cellNO*fCell1Side));
    fSolidLayer1->SetYHalfLength(0.5*((cellNO+1)*gap+cellNO*fCell1Side));

    //replace the cells, fisrt place back, then place them to the new places
    G4double adjust_displacement_x = (cellNO*fCell1Side+(cellNO+1)*gap)*0.5,
             adjust_displacement_y = (cellNO*fCell1Side+(cellNO+1)*gap)*0.5;
    G4int cellItr = 255, rowNo = 0, columnNo = 0;

    for(G4int i=0; i!=cellItr+1; ++i){

        rowNo = i%cellNO;
        columnNo = i/cellNO;

        G4ThreeVector newTrans = G4ThreeVector((rowNo+0.5)*fCell1Side+(rowNo+1)*gap
                                               -adjust_displacement_x,
                                                ((columnNo+0.5)*fCell1Side+(columnNo+1)*gap
                                               -adjust_displacement_y),
                                                0);

        fCell1Phys[i]->SetTranslation(newTrans);

    }
    for(G4int i=0; i!=cellItr+1; ++i){
        fCell1Phys[i]->CheckOverlaps();
    }

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void ComptonDetectorConstruction::SetCell2Size(G4double val){
    fCell2Side = val;
    G4int cellNO = 16;
    G4double gap = 0.*mm;

    //change the cell size
    fSolidCell2->SetXHalfLength(0.5*fCell2Side);
    fSolidCell2->SetYHalfLength(0.5*fCell2Side);

    //new layer size
    fSolidLayer2->SetXHalfLength(0.5*((cellNO+1)*gap+cellNO*fCell2Side));
    fSolidLayer2->SetYHalfLength(0.5*((cellNO+1)*gap+cellNO*fCell2Side));

    //replace the cells, fisrt place back, then place them to the new places
    G4double adjust_displacement_x = (cellNO*fCell2Side+(cellNO+1)*gap)*0.5,
             adjust_displacement_y = (cellNO*fCell2Side+(cellNO+1)*gap)*0.5;
    G4int cellItr = 255, rowNo = 0, columnNo = 0;

    for(G4int i=0; i!=cellItr+1; ++i){

        rowNo = i%cellNO;
        columnNo = i/cellNO;

        G4ThreeVector newTrans = G4ThreeVector((rowNo+0.5)*fCell2Side+(rowNo+1)*gap
                                               -adjust_displacement_x,
                                                ((columnNo+0.5)*fCell2Side+(columnNo+1)*gap
                                               -adjust_displacement_y),
                                                0);

        fCell2Phys[i]->SetTranslation(newTrans);

    }
    for(G4int i=0; i!=cellItr+1; ++i){
        fCell2Phys[i]->CheckOverlaps();
    }

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void ComptonDetectorConstruction::SetCell3Size(G4double val){
    fCell3Side = val;
    G4int cellNO = 16;
    G4double gap = 0.*mm;

    //change the cell size
    fSolidCell3->SetXHalfLength(0.5*fCell3Side);
    fSolidCell3->SetYHalfLength(0.5*fCell3Side);

    //new layer size
    fSolidLayer3->SetXHalfLength(0.5*((cellNO+1)*gap+cellNO*fCell3Side));
    fSolidLayer3->SetYHalfLength(0.5*((cellNO+1)*gap+cellNO*fCell3Side));

    //replace the cells, fisrt place back, then place them to the new places
    G4double adjust_displacement_x = (cellNO*fCell3Side+(cellNO+1)*gap)*0.5,
             adjust_displacement_y = (cellNO*fCell3Side+(cellNO+1)*gap)*0.5;
    G4int cellItr = 255, rowNo = 0, columnNo = 0;

    for(G4int i=0; i!=cellItr+1; ++i){

        rowNo = i%cellNO;
        columnNo = i/cellNO;

        G4ThreeVector newTrans = G4ThreeVector((rowNo+0.5)*fCell3Side+(rowNo+1)*gap
                                               -adjust_displacement_x,
                                                ((columnNo+0.5)*fCell3Side+(columnNo+1)*gap
                                               -adjust_displacement_y),
                                                0);

        fCell3Phys[i]->SetTranslation(newTrans);

    }
    for(G4int i=0; i!=cellItr+1; ++i){
        fCell3Phys[i]->CheckOverlaps();
    }

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

/******************** setting the physical volume of the cells ************/
void ComptonDetectorConstruction::Cell1Placement(G4int cellNO, G4double Celllen, G4String cellname,
                                                 G4double gap, G4LogicalVolume *celllog,
                                                 G4LogicalVolume *layerlog){
    G4double adjust_displacement_x = (cellNO*Celllen+(cellNO+1)*gap)*0.5,
             adjust_displacement_y = (cellNO*Celllen+(cellNO+1)*gap)*0.5;
    G4int layer_NOCopy = 0;

    for (G4int i=0; i!=cellNO; ++i) {
        for (G4int j=0; j!=cellNO; ++j) {


           G4VPhysicalVolume* cellphys = new G4PVPlacement(nullptr,
                              G4ThreeVector((i+0.5)*Celllen+(i+1)*gap
                                           -adjust_displacement_x,
                                            ((j+0.5)*Celllen+(j+1)*gap
                                           -adjust_displacement_y),
                                            0),
                              celllog,
                              cellname,
                              layerlog,
                              false,
                              layer_NOCopy,
                              true);
           fCell1Phys.push_back(cellphys);
              ++layer_NOCopy;
        }

    }

}

void ComptonDetectorConstruction::Cell2Placement(G4int cellNO, G4double Celllen, G4String cellname,
                                                 G4double gap, G4LogicalVolume *celllog,
                                                 G4LogicalVolume *layerlog
                                                 ){
    G4double adjust_displacement_x = (cellNO*Celllen+(cellNO+1)*gap)*0.5,
             adjust_displacement_y = (cellNO*Celllen+(cellNO+1)*gap)*0.5;
    G4int layer_NOCopy = 0;

    for (G4int i=0; i!=cellNO; ++i) {
        for (G4int j=0; j!=cellNO; ++j) {


           G4VPhysicalVolume* cellphys = new G4PVPlacement(nullptr,
                              G4ThreeVector((i+0.5)*Celllen+(i+1)*gap
                                           -adjust_displacement_x,
                                            ((j+0.5)*Celllen+(j+1)*gap
                                           -adjust_displacement_y),
                                            0),
                              celllog,
                              cellname,
                              layerlog,
                              false,
                              layer_NOCopy,
                              true);
           fCell2Phys.push_back(cellphys);
              ++layer_NOCopy;
        }

    }

}

void ComptonDetectorConstruction::Cell3Placement(G4int cellNO, G4double Celllen, G4String cellname,
                                                 G4double gap, G4LogicalVolume *celllog,
                                                 G4LogicalVolume *layerlog
                                                 ){
    G4double adjust_displacement_x = (cellNO*Celllen+(cellNO+1)*gap)*0.5,
             adjust_displacement_y = (cellNO*Celllen+(cellNO+1)*gap)*0.5;
    G4int layer_NOCopy = 0;

    for (G4int i=0; i!=cellNO; ++i) {
        for (G4int j=0; j!=cellNO; ++j) {


           G4VPhysicalVolume* cellphys = new G4PVPlacement(nullptr,
                              G4ThreeVector((i+0.5)*Celllen+(i+1)*gap
                                           -adjust_displacement_x,
                                            ((j+0.5)*Celllen+(j+1)*gap
                                           -adjust_displacement_y),
                                            0),
                              celllog,
                              cellname,
                              layerlog,
                              false,
                              layer_NOCopy,
                              true);
           fCell3Phys.push_back(cellphys);
              ++layer_NOCopy;
        }

    }

}

/******************** constructing the Layer and cells ***********/
void ComptonDetectorConstruction::Layer1Constructor(G4int cellNO, G4double Celllen, G4double gap,
                                                    G4double zpos, G4LogicalVolume *layermotherlog){
    //G4Material* cell_mat = G4Material::GetMaterial("Gd3CaAl4O12");
    //G4Material* cell_mat = G4Material::GetMaterial("Lu2SiO5");
    G4Material* cell_mat = G4Material::GetMaterial("CZT");
    G4Material* layer_mat = G4Material::GetMaterial("BaSO4");

    //
    // cell3
    //
    G4double cell1_x = Celllen, cell1_y = Celllen, cell1_z = fLayer1Thickness;
    fSolidCell1 = new G4Box("cell1",0.5*cell1_x,0.5*cell1_y,0.5*cell1_z);
    fCell1Log = new G4LogicalVolume(fSolidCell1,cell_mat,"cell1Log");

    //
    // layer3
    //

    G4ThreeVector pos1 = G4ThreeVector(0, 0, zpos);
    G4double layer1_dx = (cellNO+1)*gap+cellNO*Celllen,
             layer1_dy = (cellNO+1)*gap+cellNO*Celllen,
             layer1_dz  = fLayer1Thickness;


    fSolidLayer1 =
      new G4Box("layer1",                      //its name
                0.5*layer1_dx, 0.5*layer1_dy, 0.5*layer1_dz); //its size

    fLayer1Log =
      new G4LogicalVolume(fSolidLayer1,         //its solid
                          layer_mat,          //its material
                          "layer1Log");           //its name

    fLayer1Phys =  new G4PVPlacement(nullptr,                       //no rotation
                      pos1,                    //at position
                      fLayer1Log,             //its logical volume
                      "layer1",                //its name
                      layermotherlog,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      true);          //overlaps checking
}

void ComptonDetectorConstruction::Layer2Constructor(G4int cellNO, G4double Celllen, G4double gap, G4double zpos, G4LogicalVolume *layermotherlog){
    //G4Material* cell_mat = G4Material::GetMaterial("Gd3CaAl4O12");
    G4Material* cell_mat = G4Material::GetMaterial("Lu2SiO5");
    G4Material* layer_mat = G4Material::GetMaterial("BaSO4");

    //
    // cell3
    //
    G4double cell2_x = Celllen, cell2_y = Celllen, cell2_z = fLayer2Thickness;
    fSolidCell2 = new G4Box("cell1",0.5*cell2_x,0.5*cell2_y,0.5*cell2_z);
    fCell2Log = new G4LogicalVolume(fSolidCell2,cell_mat,"cell2Log");

    //
    // layer3
    //

    G4ThreeVector pos2 = G4ThreeVector(0, 0, zpos);
    G4double layer2_dx = (cellNO+1)*gap+cellNO*Celllen,
             layer2_dy = (cellNO+1)*gap+cellNO*Celllen,
             layer2_dz  = fLayer2Thickness;


    fSolidLayer2 =
      new G4Box("layer2",                      //its name
                0.5*layer2_dx, 0.5*layer2_dy, 0.5*layer2_dz); //its size

    fLayer2Log =
      new G4LogicalVolume(fSolidLayer2,         //its solid
                          layer_mat,          //its material
                          "layer2Log");           //its name

    fLayer2Phys =  new G4PVPlacement(nullptr,                       //no rotation
                      pos2,                    //at position
                      fLayer2Log,             //its logical volume
                      "layer2",                //its name
                      layermotherlog,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      true);          //overlaps checking
}

void ComptonDetectorConstruction::Layer3Constructor(G4int cellNO, G4double Celllen, G4double gap, G4double zpos, G4LogicalVolume *layermotherlog){
    //G4Material* cell_mat = G4Material::GetMaterial("Gd3CaAl4O12");
    G4Material* cell_mat = G4Material::GetMaterial("Lu2SiO5");
    G4Material* layer_mat = G4Material::GetMaterial("BaSO4");

    //
    // cell3
    //
    G4double cell3_x = Celllen, cell3_y = Celllen, cell3_z = fLayer1Thickness;
    fSolidCell3 = new G4Box("cell3",0.5*cell3_x,0.5*cell3_y,0.5*cell3_z);
    fCell3Log = new G4LogicalVolume(fSolidCell3,cell_mat,"cell3Log");

    //
    // layer3
    //

    G4ThreeVector pos3 = G4ThreeVector(0, 0, zpos);
    G4double layer3_dx = (cellNO+1)*gap+cellNO*Celllen,
             layer3_dy = (cellNO+1)*gap+cellNO*Celllen,
             layer3_dz  = fLayer1Thickness;


    fSolidLayer3 =
      new G4Box("layer1",                      //its name
                0.5*layer3_dx, 0.5*layer3_dy, 0.5*layer3_dz); //its size

    fLayer3Log =
      new G4LogicalVolume(fSolidLayer3,         //its solid
                          layer_mat,          //its material
                          "layer3Log");           //its name

    fLayer3Phys =  new G4PVPlacement(nullptr,                       //no rotation
                      pos3,                    //at position
                      fLayer3Log,             //its logical volume
                      "layer3",                //its name
                      layermotherlog,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      true);          //overlaps checking
}

void ComptonDetectorConstruction::DefineCommands(){

    fMessenger = new G4GenericMessenger(this,
                                        "/ComptonCamera/detector/",
                                        "Detector control");

    //  command setting the thickness of Layer1
    auto& Layer1ThicknessCmd
      = fMessenger->DeclareMethodWithUnit("Layer1Thickness","mm",
                                  &ComptonDetectorConstruction::SetLayer1Thickness,
                                  "Set thickness of the Layer1.");
    Layer1ThicknessCmd.SetParameterName("Thickness", true);
    Layer1ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    Layer1ThicknessCmd.SetDefaultValue("8.");

    // command setting the thickness of layer2
    auto& Layer2ThicknessCmd
      = fMessenger->DeclareMethodWithUnit("Layer2Thickness","mm",
                                  &ComptonDetectorConstruction::SetLayer2Thickness,
                                  "Set thickness of the Layer2.");
    Layer2ThicknessCmd.SetParameterName("Thickness", true);
    Layer2ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    Layer2ThicknessCmd.SetDefaultValue("8.");

    // command setting the thickness of layer3
    auto& Layer3ThicknessCmd
      = fMessenger->DeclareMethodWithUnit("Layer3Thickness","mm",
                                  &ComptonDetectorConstruction::SetLayer3Thickness,
                                  "Set thickness of the Layer3.");
    Layer3ThicknessCmd.SetParameterName("Thickness", true);
    Layer3ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    Layer3ThicknessCmd.SetDefaultValue("8.");

    auto& Cell1SizeCmd
      = fMessenger->DeclareMethodWithUnit("Cell1Size","mm",
                                  &ComptonDetectorConstruction::SetCell1Size,
                                  "Set size of the Cell1.");
    Cell1SizeCmd.SetParameterName("cellsize", true);
    Cell1SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
    Cell1SizeCmd.SetDefaultValue("2.");

    auto& Cell2SizeCmd
      = fMessenger->DeclareMethodWithUnit("Cell2Size","mm",
                                  &ComptonDetectorConstruction::SetCell2Size,
                                  "Set size of the Cell2.");
    Cell2SizeCmd.SetParameterName("cellsize", true);
    Cell2SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
    Cell2SizeCmd.SetDefaultValue("2.");

    auto& Cell3SizeCmd
      = fMessenger->DeclareMethodWithUnit("Cell3Size","mm",
                                  &ComptonDetectorConstruction::SetCell3Size,
                                  "Set size of the Cell3.");
    Cell3SizeCmd.SetParameterName("cellsize", true);
    Cell3SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
    Cell3SizeCmd.SetDefaultValue("2.");

    auto& Layer1ZPositionCmd
      = fMessenger->DeclareMethodWithUnit("Layer1Z","mm",
                                  &ComptonDetectorConstruction::SetLayer1ZPosition,
                                  "Set z position of the Layer1.");
    Layer1ZPositionCmd.SetParameterName("ZPosition", true);
    Layer1ZPositionCmd.SetRange("ZPosition>=-100. && ZPosition<100.");
    Layer1ZPositionCmd.SetDefaultValue("-18");

    auto& Layer2ZPositionCmd
      = fMessenger->DeclareMethodWithUnit("Layer2Z","mm",
                                  &ComptonDetectorConstruction::SetLayer2ZPosition,
                                  "Set z position of the Layer2.");
    Layer2ZPositionCmd.SetParameterName("ZPosition", true);
    Layer2ZPositionCmd.SetRange("ZPosition>=-50. && ZPosition<50.");
    Layer2ZPositionCmd.SetDefaultValue("0");

    auto& Layer3ZPositionCmd
      = fMessenger->DeclareMethodWithUnit("Layer3Z","mm",
                                  &ComptonDetectorConstruction::SetLayer3ZPosition,
                                  "Set z position of the Layer3.");
    Layer3ZPositionCmd.SetParameterName("ZPosition", true);
    Layer3ZPositionCmd.SetRange("ZPosition>=-100. && ZPosition<100.");
    Layer3ZPositionCmd.SetDefaultValue("18");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
