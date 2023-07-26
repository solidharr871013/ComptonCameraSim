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
/// \file ComptonPrimaryGeneratorAction.cc
/// \brief Implementation of the ComptonPrimaryGeneratorAction class

#include "ComptonPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "random"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonPrimaryGeneratorAction::ComptonPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fMessenger(nullptr),
  fParticleGun(nullptr),
  xpos(0.),
  ypos(40),
  zpos(-40.)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  DefineCommands();

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(662*keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonPrimaryGeneratorAction::~ComptonPrimaryGeneratorAction()
{
    delete  fMessenger;
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.






    //the following code is for generating two seperated source, the source is aparted in
    // the distance of source_divergence in x-axis.


    //G4int randNumber = rand()%2;
    G4double randNumber = G4UniformRand();

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);

    G4double cellNo = 20, cellSide = 2.*mm;
    G4double source_x = 10*cm, source_y = 0*cm, source_z = -50*cm;
    G4double detectDim = 0.5*cellNo*cellSide;

    G4double targetX = detectDim*(2*G4UniformRand()-1),
             targetY = detectDim*(2*G4UniformRand()-1),
             targetZ = 25*mm;




/*********************************/
    if(randNumber<0.7){

        source_x=0*cm;
        source_y=0*cm;
        source_z=-50*cm;

        /****cover the detector, any direction****

        G4double philimitL = std::sqrt(source_x*source_x+source_y*source_y)+detectDim;
        G4double phicosLimL = -1*source_z/(std::sqrt(philimitL*philimitL+source_z*source_z));

        std::uniform_real_distribution<double> phiUniDistr(phicosLimL,1);

        std::uniform_real_distribution<double> thetaUniDistr(0,2*pi);

        G4double cos_theta = std::cos(thetaUniDistr(generator)),
                 sin_theta = std::sin(thetaUniDistr(generator));

        G4double cos_phi = phiUniDistr(generator),
                 sin_phi = std::sqrt(1-cos_phi*cos_phi);
        G4double x_direction = 1.*sin_phi*cos_theta,
                 y_direction = 1.*sin_phi*sin_theta,
                 z_direction = 1.*cos_phi;

        ********/

        G4double x_direction = targetX-source_x,
                 y_direction = targetY-source_y,
                 z_direction = targetZ-source_z;

        fParticleGun->SetParticleEnergy(662*keV);

        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x_direction,y_direction,z_direction));

        fParticleGun->SetParticlePosition(G4ThreeVector(source_x,source_y,source_z));

        fParticleGun->GeneratePrimaryVertex(anEvent);

    }
    else if (randNumber>=0.7) {

        /****cover the detector, any direction****

        G4double philimitH = std::sqrt(source_x*source_x+source_y*source_y)+detectDim,
                 philimitL = std::sqrt(source_x*source_x+source_y*source_y)-detectDim;

        G4double phicosLimH = -1*source_z/(std::sqrt(philimitH*philimitH+source_z*source_z)),
                 phicosLimL = -1*source_z/(std::sqrt(philimitL*philimitL+source_z*source_z));

        std::uniform_real_distribution<double> phiUniDistr(phicosLimL,phicosLimH);

        G4double thetalimitL_x = source_x+0.5*cellNo*cellSide,
                 thetalimitH_x = source_x-0.5*cellNo*cellSide,
                 thetalimitL_y = source_y-0.5*cellNo*cellSide,
                 thetalimitH_y = source_y+0.5*cellNo*cellSide;

        G4double thetacosLimH =
                thetalimitH_x/(std::sqrt(thetalimitH_x*thetalimitH_x+thetalimitH_y*thetalimitH_y)),
                 thetacosLimL =
                thetalimitL_x/(std::sqrt(thetalimitL_x*thetalimitL_x+thetalimitL_y*thetalimitL_y));

        std::uniform_real_distribution<double> thetaUniDistr(thetacosLimH,thetacosLimL);


        G4double cos_theta = -1*thetaUniDistr(generator),
                 sin_theta = -1*std::sqrt(1-cos_theta*cos_theta);

        G4double cos_phi = phiUniDistr(generator),
                 sin_phi = std::sqrt(1-cos_phi*cos_phi);

        G4double x_direction = 1.*sin_phi*cos_theta,
                 y_direction = 1.*sin_phi*sin_theta,
                 z_direction = 1.*cos_phi;

        ********/

        G4double x_direction = targetX-source_x,
                 y_direction = targetY-source_y,
                 z_direction = targetZ-source_z;

        fParticleGun->SetParticleEnergy(662*keV);

        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x_direction,y_direction,z_direction));

        fParticleGun->SetParticlePosition(G4ThreeVector(source_x,source_y,source_z));

        fParticleGun->GeneratePrimaryVertex(anEvent);
    }


/*********************/
}

void ComptonPrimaryGeneratorAction::DefineCommands(){
    fMessenger = new G4GenericMessenger(this,
                                        "/ComptonCamera/gunpos/",
                                        "gun Control");


    auto& xposCmd
      = fMessenger->DeclareMethodWithUnit("xpos","cm",
                                  &ComptonPrimaryGeneratorAction::setPosX,
                                  "Set x position for source.");
    xposCmd.SetParameterName("xPosition", true);
    xposCmd.SetRange("xPosition>=-10000. && xPosition<10000.");
    xposCmd.SetDefaultValue("0");

    auto& yposCmd
      = fMessenger->DeclareMethodWithUnit("ypos","cm",
                                  &ComptonPrimaryGeneratorAction::setPosY,
                                  "Set y position for source.");
    yposCmd.SetParameterName("yPosition", true);
    yposCmd.SetRange("yPosition>=-10000. && yPosition<10000.");
    yposCmd.SetDefaultValue("40");

    auto& zposCmd
      = fMessenger->DeclareMethodWithUnit("zpos","cm",
                                  &ComptonPrimaryGeneratorAction::setPosZ,
                                  "Set z position for source.");
    zposCmd.SetParameterName("zPosition", true);
    zposCmd.SetRange("zPosition>=-10000. && zPosition<10000.");
    zposCmd.SetDefaultValue("-40.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

