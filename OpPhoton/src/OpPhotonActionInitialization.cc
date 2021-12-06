#include "OpPhotonActionInitialization.hh"
#include "OpPhotonPrimaryGeneratorAction.hh"
#include "OpPhotonRunAction.hh"
#include "OpPhotonEventAction.hh"
#include "OpPhotonSteppingAction.hh"
#include "OpPhotonStackingAction.hh"

OpPhotonActionInitialization::OpPhotonActionInitialization()
 : G4VUserActionInitialization()
{}

OpPhotonActionInitialization::~OpPhotonActionInitialization()
{}

void OpPhotonActionInitialization::BuildForMaster() const
{
  SetUserAction(new OpPhotonRunAction);
}

void OpPhotonActionInitialization::Build() const
{
  OpPhotonSteppingAction* steppingAction = new OpPhotonSteppingAction();

  OpPhotonEventAction* eventAction = new OpPhotonEventAction(steppingAction);

  SetUserAction(new OpPhotonPrimaryGeneratorAction);
  SetUserAction(new OpPhotonRunAction);
  SetUserAction(eventAction);
  SetUserAction(steppingAction);
  SetUserAction(new OpPhotonStackingAction);  
}
