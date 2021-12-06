#include "G4Timer.hh"

#include "OpPhotonRunAction.hh"
#include "OpPhotonAnalysis.hh"

#include "G4Run.hh"

OpPhotonRunAction::OpPhotonRunAction()
 : G4UserRunAction(),
   fTimer(0)
{
  fTimer = new G4Timer;

  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in OpPhotonAnalysis.hh
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  G4cout << "Using " << mgr->GetType() << G4endl;

  mgr->SetVerboseLevel(0);
  mgr->SetNtupleMerging(true);

  // Data stored in the EventAction
  mgr->CreateNtuple("OpPhoton", "Compton and optical photon data");
  mgr->CreateNtupleDColumn("ComptonE");
  mgr->CreateNtupleDColumn("ComptonDepth");
  mgr->CreateNtupleDColumn("nbFrontPhot");
  mgr->CreateNtupleDColumn("timeFront20");
  mgr->CreateNtupleDColumn("timeFront100");
  mgr->CreateNtupleDColumn("nbBackPhot");
  mgr->CreateNtupleDColumn("timeBack20");
  mgr->CreateNtupleDColumn("timeBack100");
  mgr->FinishNtuple();
}

OpPhotonRunAction::~OpPhotonRunAction()
{
  delete fTimer;
  delete G4AnalysisManager::Instance();
}

void OpPhotonRunAction::BeginOfRunAction(const G4Run*)
{
    // Get analysis manager
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();

  // Open an output file
  G4String fileName = "8sides_10cm_0975";
  mgr->OpenFile(fileName);

  fTimer->Start();
}

void OpPhotonRunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent()
         << " " << *fTimer << G4endl;
  
  // Print histogram statistics
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();

  // Save histograms & ntuple
  mgr->Write();
  mgr->CloseFile();
}
