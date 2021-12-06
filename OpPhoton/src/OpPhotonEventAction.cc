// #include "G4EventManager.hh"
// #include "Randomize.hh"

#include "OpPhotonEventAction.hh"

#include "OpPhotonRunAction.hh"
#include "OpPhotonVSteppingAction.hh"

OpPhotonEventAction::OpPhotonEventAction(
  OpPhotonVSteppingAction* onePhotonSteppingAction)
:fpOpPhotonVSteppingAction(onePhotonSteppingAction)
{}

OpPhotonEventAction::~OpPhotonEventAction()
{}

void OpPhotonEventAction::BeginOfEventAction(const G4Event*)
{
  fpOpPhotonVSteppingAction->BeginOfEventAction();
}

void OpPhotonEventAction::EndOfEventAction(const G4Event*)
{
  fpOpPhotonVSteppingAction->EndOfEventAction();
}
