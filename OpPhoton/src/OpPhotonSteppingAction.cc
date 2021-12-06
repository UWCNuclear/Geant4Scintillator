#include "OpPhotonSteppingAction.hh"
#include "OpPhotonEventAction.hh"
#include "OpPhotonAnalysis.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "OpPhotonRunAction.hh"

#include "G4SystemOfUnits.hh"

OpPhotonSteppingAction::OpPhotonSteppingAction()
: fComptonE(0.)
, fComptonDepth(0.)
, fComptonTime(0.)
, fNbFrontPhot(0)
, fV_timeFront()
, fNbBackPhot(0)
, fV_timeBack()
{}

void OpPhotonSteppingAction::BeginOfEventAction()
{
  fComptonE = -1.;
  fComptonDepth = -1.;
  fComptonTime = -1.;

  fNbFrontPhot = 0;
  fV_timeFront.clear();
  fNbBackPhot = 0;
  fV_timeBack.clear();
}

void OpPhotonSteppingAction::EndOfEventAction()
{  
  if (fComptonE > 0.)
  {
    // The time of the front/back end detector is assumed to be the time of 
    // arrival of the 20th/100th optical photon (need to be experimentally 
    // verified).
    // To select the 20th/100th optical gamma arriving, the following
    // sorts the element of the vector in increasing order:
    sort(fV_timeFront.begin(), fV_timeFront.end());
    sort(fV_timeBack.begin(), fV_timeBack.end());

    G4double timeFront20 = -1.;
    G4double timeFront100 = -1.;
    G4double timeBack20 = -1.;
    G4double timeBack100 = -1.;

    // Checks if there was at least 20 optical photons reaching either end 
    // during the event
    if (fV_timeFront.size() >= 20 && fV_timeBack.size() >= 20)
    {
      timeFront20 = fV_timeFront.at(19);
      timeBack20 = fV_timeBack.at(19);
    }

    // Checks if there was at least 100 optical photons reaching either end 
    // during the event
    if (fV_timeFront.size() >= 100 && fV_timeBack.size() >= 100)
    {
      timeFront100 = fV_timeFront.at(99);
      timeBack100 = fV_timeBack.at(99);
    }

    G4AnalysisManager* mgr = G4AnalysisManager::Instance();
    mgr->FillNtupleDColumn(0, fComptonE);
    mgr->FillNtupleDColumn(1, fComptonDepth);
    mgr->FillNtupleDColumn(2, fNbFrontPhot);
    mgr->FillNtupleDColumn(3, timeFront20);
    mgr->FillNtupleDColumn(4, timeFront100);
    mgr->FillNtupleDColumn(5, fNbBackPhot);
    mgr->FillNtupleDColumn(6, timeBack20);
    mgr->FillNtupleDColumn(7, timeBack100);
    mgr->AddNtupleRow();
  }
}

void OpPhotonSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step->GetTrack();

  G4String particleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();
  
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4StepPoint* postStepPoint = step->GetPostStepPoint();

  if (postStepPoint->GetPhysicalVolume() == nullptr) return;
  const G4String postVolume = postStepPoint->GetPhysicalVolume()->GetName();

  if (postStepPoint->GetProcessDefinedStep() == nullptr) return;
  const G4String postProcess =
   postStepPoint->GetProcessDefinedStep()->GetProcessName();

  // Records the annihilation gamma Compton scattering.
  if (step->GetTrack()->GetTrackID() == 1 && postVolume == "Stick"
      && postProcess == "compt")
  {
    fComptonE = (preStepPoint->GetKineticEnergy()
                 - postStepPoint->GetKineticEnergy()) / keV;

    fComptonDepth = postStepPoint->GetPosition().z() / cm;
    fComptonTime = step->GetTrack()->GetGlobalTime();
    
    // WARNING: This stops the tracking of the annihilation gamma
    // We are focusing on its first interaction.
    track->SetTrackStatus(fStopAndKill);
  }

  // Checks if the optical gamma is reaching either stick end
  if (particleName == "opticalphoton" && fComptonE > 0.)
  {
    G4double TOA = (step->GetTrack()->GetGlobalTime() - fComptonTime) / ns;

    if (postVolume == "StickFrontEnd")
    {
      // Count the numbre of optical photons reaching the front end
      fNbFrontPhot++;
      fV_timeFront.push_back(TOA);
    } 
      
    if (postVolume == "StickBackEnd")
    {
      // Count the numbre of optical photons reaching the front end
      fNbBackPhot++;
      fV_timeBack.push_back(TOA);
    }

    return;
  }
}