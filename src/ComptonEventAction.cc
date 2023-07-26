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
/// \file ComptonEventAction.cc
/// \brief Implementation of the ComptonEventAction class

#include "ComptonEventAction.hh"
#include "ComptonRunAction.hh"
#include "ComptonSteppingAction.hh"
#include "ComptonAnalysis.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4THitsMap.hh"

#include "ComptonLayer1Hit.hh"
#include "ComptonLayer1SD.hh"

#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

#include "random"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl;
      G4Exception("ComptonEventAction::EndOfEventAction()",
                  "ComptonCode001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl;
    G4Exception("ComptonEventAction::EndOfEventAction()",
                "ComptonCode001", JustWarning, msg);
  }
  return hc;
}
}

ComptonEventAction::ComptonEventAction(ComptonRunAction* runAction)
: G4UserEventAction(),
  fMessenger(nullptr),
  fRunAction(runAction),
  fEdep(0.),
  f3LayerCoincidence(0),
  fTotalNumber(0),
  fOrdered3(0),
  fOrdered2(0),
  f2LayerCoincidence(0),
  fLayer1ID(-1),
  fLayer2ID(-1),
  fLayer3ID(-1),
  f2TotalE_low(0.642),
  f2TotalE_high(0.682),
  f2ScatterE_low(0.01),
  f2ScatterE_high(0.165)
{
    DefineCommands();
     G4RunManager::GetRunManager()->SetPrintProgress(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ComptonEventAction::~ComptonEventAction()
{ delete fMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonEventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
  f3LayerCoincidence = 0;
  fTotalNumber = 0;
  fOrdered3 = 0;
  fOrdered2 = 0;
  f2LayerCoincidence = 0;

  if(fLayer1ID == -1 ){
      auto sdManager = G4SDManager::GetSDMpointer();

      fLayer1ID = sdManager->GetCollectionID("Layer1/LayerColl");
      fLayer2ID = sdManager->GetCollectionID("Layer2/LayerColl");
      fLayer3ID = sdManager->GetCollectionID("Layer3/LayerColl");
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComptonEventAction::EndOfEventAction(const G4Event* event)
{   
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();

    G4double eThreshold = 1*keV;

    auto hcLayer1 = GetHC(event, fLayer1ID);
    if(!hcLayer1) return;
    auto hcLayer2 = GetHC(event, fLayer2ID);
    if(!hcLayer2) return;
    auto hcLayer3 = GetHC(event, fLayer3ID);
    if(!hcLayer3) return;

/*******************************************************************/
 //Three Layer coincidence

    
    // if the photon enter the first layer, total number added, used for efficiency 
    // calculation.
    if(hcLayer1->GetSize()!=0){

        ++fTotalNumber;
        fRunAction->AddTotalNumber(fTotalNumber);

        

        if(hcLayer2->GetSize()!=0 && hcLayer3->GetSize()!=0){

            G4int eventID =
                  G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

            G4double EnergyDepositeInLayer1 = 0, EnergyDepositeInLayer2 = 0,
                     EnergyDepositeInLayer3 = 0,
                     firstHitTime1 = 0, firstHitTime2 = 0, firstHitTime3 = 0;

                // for recording the energy desposited in the layer, 
                // and get the time of the first interaction in the layer
            for(unsigned long itr = 0; itr != hcLayer1->GetSize(); ++itr){

                auto hit = static_cast<ComptonLayer1Hit*>(hcLayer1->GetHit(itr));

                G4double bufferTime = 0;
                if(hit->GetEdep()!=0){

                    EnergyDepositeInLayer1 += hit->GetEdep();

                    if(firstHitTime1==0){firstHitTime1 = hit->GetTime();}
                    else {
                        bufferTime = hit->GetTime();
                        if(firstHitTime1>=bufferTime){
                           firstHitTime1 = bufferTime;
                        }
                    }
                }

            }



            for(unsigned long itr = 0; itr != hcLayer2->GetSize(); ++itr){

                auto hit = static_cast<ComptonLayer1Hit*>(hcLayer2->GetHit(itr));

                G4double bufferTime = 0;

                if(hit->GetEdep()!=0){

                    EnergyDepositeInLayer2 += hit->GetEdep();

                    if(firstHitTime2==0){firstHitTime2 = hit->GetTime();}
                    else {

                        bufferTime = hit->GetTime();
                        if(firstHitTime2>=bufferTime){
                            firstHitTime2 = bufferTime;
                        }
                    }
                }

            }



            for(unsigned long itr = 0; itr != hcLayer3->GetSize(); ++itr){

                auto hit = static_cast<ComptonLayer1Hit*>(hcLayer3->GetHit(itr));

                G4double bufferTime = 0;

                if(hit->GetEdep()!=0){

                    EnergyDepositeInLayer3 += hit->GetEdep();

                    if(firstHitTime3==0){firstHitTime3 = hit->GetTime();}
                    else {

                        bufferTime = hit->GetTime();

                        if(firstHitTime3>=bufferTime){

                            firstHitTime3 = bufferTime;
                        }
                    }
                }

            }

            //G4cout << EnergyDepositeInLayer1 << " " << EnergyDepositeInLayer2 << " " << EnergyDepositeInLayer3 << G4endl;

            G4double eDep1_fluc = AddFluction(EnergyDepositeInLayer1),
                     eDep2_fluc = AddFluction(EnergyDepositeInLayer2),
                     eDep3_fluc = AddFluction(EnergyDepositeInLayer3);

            if(eDep1_fluc>10*keV && eDep2_fluc<f2ScatterE_high
                    && eDep2_fluc>10*keV && eDep3_fluc>10*keV){
                 
                for (unsigned long iter=0; iter!=hcLayer1->GetSize(); ++iter) {
                    auto hit = static_cast<ComptonLayer1Hit*>(hcLayer1->GetHit(iter));
                        
                   // if(hit->GetParticleName()=="gamma"
                        //    && hit->GetProcessName()!="Transportation"){
    
                        analysis->FillNtupleIColumn(0,0,eventID);
    
                        analysis->FillNtupleSColumn(0,1,hit->GetParticleName());
                        analysis->FillNtupleDColumn(0,2,eDep1_fluc);
                        analysis->FillNtupleDColumn(0,3,hit->GetScatteringPos().x());
                        analysis->FillNtupleDColumn(0,4,hit->GetScatteringPos().y());
                        analysis->FillNtupleDColumn(0,5,hit->GetScatteringPos().z());
                        analysis->FillNtupleDColumn(0,6,hit->GetTime());
                        analysis->FillNtupleSColumn(0,7,hit->GetProcessName());
                        analysis->AddNtupleRow(0);
    
                    //}
                }
                    
                for (unsigned long iter=0; iter!=hcLayer2->GetSize(); ++iter) {
                    auto hit = static_cast<ComptonLayer1Hit*>(hcLayer2->GetHit(iter));
                        
                    //if(hit->GetParticleName()=="gamma"
                      //      && hit->GetProcessName()!="Transportation"){
                        analysis->FillNtupleIColumn(0,0,eventID);
    
                        analysis->FillNtupleSColumn(0,8,hit->GetParticleName());
                        analysis->FillNtupleDColumn(0,9,eDep2_fluc);
                        analysis->FillNtupleDColumn(0,10,hit->GetScatteringPos().x());
                        analysis->FillNtupleDColumn(0,11,hit->GetScatteringPos().y());
                        analysis->FillNtupleDColumn(0,12,hit->GetScatteringPos().z());
                        analysis->FillNtupleDColumn(0,13,hit->GetTime());
                        analysis->FillNtupleSColumn(0,14,hit->GetProcessName());
                        analysis->AddNtupleRow(0);
                    //}
                }
                    
                for (unsigned long iter=0; iter!=hcLayer3->GetSize(); ++iter) {
                    auto hit = static_cast<ComptonLayer1Hit*>(hcLayer3->GetHit(iter));
                        
                    //if(hit->GetParticleName()=="gamma"
                      //      && hit->GetProcessName()!="Transportation"){
                        analysis->FillNtupleIColumn(0,0,eventID);
    
                        analysis->FillNtupleSColumn(0,15,hit->GetParticleName());
                        analysis->FillNtupleDColumn(0,16,eDep3_fluc);
                        analysis->FillNtupleDColumn(0,17,hit->GetScatteringPos().x());
                        analysis->FillNtupleDColumn(0,18,hit->GetScatteringPos().y());
                        analysis->FillNtupleDColumn(0,19,hit->GetScatteringPos().z());
                        analysis->FillNtupleDColumn(0,20,hit->GetTime());
                        analysis->FillNtupleSColumn(0,21,hit->GetProcessName());
                        analysis->AddNtupleRow(0);
                  //  }
                }

                    
                ++f3LayerCoincidence;
                fRunAction->Add3LayerCoincidence(f3LayerCoincidence);
                    
                    
                if(firstHitTime1<firstHitTime2
                    && firstHitTime2<firstHitTime3){

                    ++fOrdered3;
                    fRunAction->AddOrdered3(fOrdered3);

                }

            }
            analysis->FillH1(0,eDep1_fluc);
            analysis->FillH1(1,eDep2_fluc);
            analysis->FillH1(2,eDep3_fluc);

            //analysis->FillH2(1,eDep1_fluc,eDep2_fluc);
            //analysis->FillH2(2,eDep2_fluc,eDep3_fluc);
        }
    }

/**************************************************************/

/*****************************************************************/
//two Layer coincidence
//in this coincidence, the third Layer is set as absorber, and the first and second Layer
//is scatterer.
    if(hcLayer3->GetSize()!=0){

        G4double  EdepInLayer1 = 0, EdepInLayer2 = 0, EdepInLayer3 = 0;
        G4double EdepInTwoLayer = 0,
                 firstHitTime1 = 0, firstHitTime2 = 0, firstHitTime3 = 0;

        G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

        for(unsigned long itr = 0; itr != hcLayer1->GetSize(); ++itr){

            G4double bufferTime = 0.;
            auto hit = static_cast<ComptonLayer1Hit*>(hcLayer1->GetHit(itr));

            if(hit->GetEdep()!=0){
                EdepInLayer1 += hit->GetEdep();
                bufferTime = hit->GetTime();
                if(firstHitTime1==0){firstHitTime1 = hit->GetTime();}
                else {
                    if(bufferTime<=firstHitTime1){
                        firstHitTime1 = bufferTime;
                    }
                }
            }

        }

        for(unsigned long itr = 0; itr != hcLayer2->GetSize(); ++itr){

            G4double bufferTime = 0.;
            auto hit = static_cast<ComptonLayer1Hit*>(hcLayer2->GetHit(itr));
            if(hit->GetEdep()!=0){
                EdepInLayer2 += hit->GetEdep();
                bufferTime = hit->GetTime();
                if(firstHitTime2==0){firstHitTime2 = hit->GetTime();}
                else {
                    if(bufferTime<=firstHitTime2){
                        firstHitTime2 = bufferTime;
                    }
                }
            }
        }

        for(unsigned long itr = 0; itr != hcLayer3->GetSize(); ++itr){

            G4double bufferTime = 0.;
            auto hit = static_cast<ComptonLayer1Hit*>(hcLayer3->GetHit(itr));
            if(hit->GetEdep()!=0){
                EdepInLayer3 += hit->GetEdep();
                bufferTime = hit->GetTime();
                if(firstHitTime3==0){firstHitTime3 = hit->GetTime();}
                else {
                    if(bufferTime<=firstHitTime3){
                        firstHitTime3 = bufferTime;
                    }
                }
            }

        }
        //G4cout << EdepInLayer1 << " " << EdepInLayer2 << " " << EdepInLayer3 << G4endl;



        G4double eDep1_fluc = AddFluction(EdepInLayer1),
                 eDep2_fluc = AddFluction(EdepInLayer2),
                 eDep3_fluc = AddFluction(EdepInLayer3);


        if(eDep1_fluc > eThreshold && eDep2_fluc == 0.*keV && eDep3_fluc > eThreshold){

            EdepInTwoLayer = eDep1_fluc + eDep3_fluc;

            if(EdepInTwoLayer<f2TotalE_high && EdepInTwoLayer>f2TotalE_low
                    && eDep1_fluc<f2ScatterE_high
                    && eDep1_fluc>f2ScatterE_low){

                for (unsigned long iter=0; iter!=hcLayer1->GetSize(); ++iter) {
                    auto hit = static_cast<ComptonLayer1Hit*>(hcLayer1->GetHit(iter));

                   // if(hit->GetParticleName()=="gamma"
                      //      && hit->GetProcessName()!="Transportation"){

                        analysis->FillNtupleIColumn(1,0,eventID);

                        analysis->FillNtupleSColumn(1,1,hit->GetParticleName());
                        analysis->FillNtupleDColumn(1,2,eDep1_fluc);
                        analysis->FillNtupleDColumn(1,3,hit->GetScatteringPos().x());
                        analysis->FillNtupleDColumn(1,4,hit->GetScatteringPos().y());
                        analysis->FillNtupleDColumn(1,5,hit->GetScatteringPos().z());
                        analysis->FillNtupleDColumn(1,6,hit->GetTime());
                        analysis->FillNtupleSColumn(1,7,hit->GetProcessName());
                        analysis->AddNtupleRow(1);
                    //}


                }

                for (unsigned long iter=0; iter!=hcLayer3->GetSize(); ++iter) {
                    auto hit = static_cast<ComptonLayer1Hit*>(hcLayer3->GetHit(iter));

                   // if(hit->GetParticleName()=="gamma"
                     //       && hit->GetProcessName()!="Transportation"){
                        analysis->FillNtupleIColumn(1,0,eventID);

                        analysis->FillNtupleSColumn(1,8,hit->GetParticleName());
                        analysis->FillNtupleDColumn(1,9,eDep3_fluc);
                        analysis->FillNtupleDColumn(1,10,hit->GetScatteringPos().x());
                        analysis->FillNtupleDColumn(1,11,hit->GetScatteringPos().y());
                        analysis->FillNtupleDColumn(1,12,hit->GetScatteringPos().z());
                        analysis->FillNtupleDColumn(1,13,hit->GetTime());
                        analysis->FillNtupleSColumn(1,14,hit->GetProcessName());
                        analysis->AddNtupleRow(1);
                   // }


                }
                ++f2LayerCoincidence;
                fRunAction->Add2LayerCoincidence(f2LayerCoincidence);

                if(firstHitTime1<firstHitTime3){
                    ++fOrdered2;
                    fRunAction->AddOrdered2(fOrdered2);
                }
            }

            analysis->FillH1(3,eDep1_fluc);
            analysis->FillH1(4,eDep3_fluc);

            //analysis->FillH2(0,eDep1_fluc,eDep3_fluc);

        }
        else if(eDep1_fluc == 0. && eDep2_fluc > eThreshold
                && eDep3_fluc > eThreshold){

            EdepInTwoLayer = eDep2_fluc + eDep3_fluc;

            if(EdepInTwoLayer<f2TotalE_high && EdepInTwoLayer>f2TotalE_low
                    && eDep2_fluc<f2ScatterE_high
                    && eDep2_fluc>f2ScatterE_low){

                for (unsigned long iter=0; iter!=hcLayer2->GetSize(); ++iter) {
                    auto hit = static_cast<ComptonLayer1Hit*>(hcLayer2->GetHit(iter));

                   // if(hit->GetParticleName()=="gamma"
                    //        && hit->GetProcessName()!="Transportation"){

                        analysis->FillNtupleIColumn(1,0,eventID);

                        analysis->FillNtupleSColumn(1,1,hit->GetParticleName());
                        analysis->FillNtupleDColumn(1,2,eDep2_fluc);
                        analysis->FillNtupleDColumn(1,3,hit->GetScatteringPos().x());
                        analysis->FillNtupleDColumn(1,4,hit->GetScatteringPos().y());
                        analysis->FillNtupleDColumn(1,5,hit->GetScatteringPos().z());
                        analysis->FillNtupleDColumn(1,6,hit->GetTime());
                        analysis->FillNtupleSColumn(1,7,hit->GetProcessName());
                        analysis->AddNtupleRow(1);
                   // }


                }

                for (unsigned long iter=0; iter!=hcLayer3->GetSize(); ++iter) {
                    auto hit = static_cast<ComptonLayer1Hit*>(hcLayer3->GetHit(iter));

                  //  if(hit->GetParticleName()=="gamma"
                   //         && hit->GetProcessName()!="Transportation"){
                        analysis->FillNtupleIColumn(1,0,eventID);

                        analysis->FillNtupleSColumn(1,8,hit->GetParticleName());
                        analysis->FillNtupleDColumn(1,9,eDep3_fluc);
                        analysis->FillNtupleDColumn(1,10,hit->GetScatteringPos().x());
                        analysis->FillNtupleDColumn(1,11,hit->GetScatteringPos().y());
                        analysis->FillNtupleDColumn(1,12,hit->GetScatteringPos().z());
                        analysis->FillNtupleDColumn(1,13,hit->GetTime());
                        analysis->FillNtupleSColumn(1,14,hit->GetProcessName());
                        analysis->AddNtupleRow(1);
                   // }


                }

                ++f2LayerCoincidence;
                fRunAction->Add2LayerCoincidence(f2LayerCoincidence);

                if(firstHitTime2<firstHitTime3){
                    ++fOrdered2;
                    fRunAction->AddOrdered2(fOrdered2);
                }
            }

                analysis->FillH1(3,eDep2_fluc);
                analysis->FillH1(4,eDep3_fluc);

                //analysis->FillH2(0,eDep2_fluc,eDep3_fluc);

        }
        else {
            EdepInTwoLayer = 0.;
        }

    }

/*****************************************************************/
    //G4int test = ComptonEventAction::GetTotalNumber();
    //G4cout << "now the total number is " << test << G4endl;


  // accumulate statistics in run action
   // if(fEdep>eThreshold){
     //  fRunAction->AddEdep(fEdep);
       //analysis->FillH1(0,fEdep);
   // }

    if(hcLayer2){

        G4double edepInLayer2 = 0.;

        for (unsigned long itr=0; itr != hcLayer2->GetSize(); ++itr) {

            auto hit = static_cast<ComptonLayer1Hit*>(hcLayer2->GetHit(itr));
            edepInLayer2 += hit->GetEdep();
        }
        if(edepInLayer2>eThreshold){
            analysis->FillH1(0,AddFluction(edepInLayer2));
        }
    }

}

G4double ComptonEventAction::AddFluction(G4double val){

    //calculate the coefficient of the resolution, based on 662keV
    G4double Resolution662 =0.08,
             coefficient = Resolution662*std::sqrt(0.662);

    G4double valVar = coefficient*std::sqrt(val)/(2*std::sqrt(2*std::log(2)));
    //G4double valVar = 0.08*val/(2*std::sqrt(2*std::log(2)));


    std::random_device rd;
    std::mt19937 generator(rd());

    std::normal_distribution<G4double> distribution(val,valVar);
    G4double var_new = distribution(generator);

    return var_new;
}


void ComptonEventAction::Set2TotalELow(G4double val){f2TotalE_low = val;}
void ComptonEventAction::Set2TotalEHigh(G4double val){f2TotalE_high = val;}
void ComptonEventAction::Set2ScatterELow(G4double val){f2ScatterE_low = val;}
void ComptonEventAction::Set2ScatterEHigh(G4double val){f2ScatterE_high = val;}

void ComptonEventAction::DefineCommands(){
    fMessenger = new G4GenericMessenger(this,
                                        "/ComptonCamera/event/",
                                        "Event Control");


    auto& ThresTotalEHighCmd
      = fMessenger->DeclareMethodWithUnit("ThresTotalEHigh","MeV",
                                  &ComptonEventAction::Set2TotalEHigh,
                                  "Set high threshold of total edep in a 2coincident event.");
    ThresTotalEHighCmd.SetParameterName("Energy", true);
    ThresTotalEHighCmd.SetRange("Energy>=0. && Energy<100.");
    ThresTotalEHighCmd.SetDefaultValue("0.710");

    auto& ThresTotalELowCmd
      = fMessenger->DeclareMethodWithUnit("ThresTotalELow","MeV",
                                  &ComptonEventAction::Set2TotalELow,
                                  "Set low threshold of total edep in a 2coincident event.");
    ThresTotalELowCmd.SetParameterName("Energy", true);
    ThresTotalELowCmd.SetRange("Energy>=0. && Energy<100.");
    ThresTotalELowCmd.SetDefaultValue("0.615");

    auto& ThresScatterEHighCmd
      = fMessenger->DeclareMethodWithUnit("ThresScatterEHigh","MeV",
                                  &ComptonEventAction::Set2ScatterEHigh,
                                  "Set high threshold of Scatter layer edep in a 2coincident event.");
    ThresScatterEHighCmd.SetParameterName("Energy", true);
    ThresScatterEHighCmd.SetRange("Energy>=0. && Energy<100.");
    ThresScatterEHighCmd.SetDefaultValue("0.165");

    auto& ThresScatterELowCmd
      = fMessenger->DeclareMethodWithUnit("ThresScatterELow","MeV",
                                  &ComptonEventAction::Set2ScatterELow,
                                  "Set Low threshold of Scatter layer edep in a 2coincident event.");
    ThresScatterELowCmd.SetParameterName("Energy", true);
    ThresScatterELowCmd.SetRange("Energy>=0. && Energy<100.");
    ThresScatterELowCmd.SetDefaultValue("0.01");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
