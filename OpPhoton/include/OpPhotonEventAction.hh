#ifndef OpPhotonEventAction_h
#define OpPhotonEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class OpPhotonRunAction;
class OpPhotonVSteppingAction;

class OpPhotonEventAction : public G4UserEventAction
{
  public:
  
    OpPhotonEventAction(OpPhotonVSteppingAction*);
    virtual ~OpPhotonEventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    
  private:
    
    OpPhotonVSteppingAction* fpOpPhotonVSteppingAction;
};

#endif

    
