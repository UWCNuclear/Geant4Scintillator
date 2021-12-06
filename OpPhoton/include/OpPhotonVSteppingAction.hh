// An intermediate inheritable class for stepping action.  It introduces a
// BeginOfEventAction and an EndOfEventAction - very useful.  So any
// stepping action in this project should inherit.
//
// BeginOfEventAction and EndOfEventAction are called from Tangle2EventAction.
#ifndef OpPhotonVSteppingAction_hh
#define OpPhotonVSteppingAction_hh

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class OpPhotonVSteppingAction : public G4UserSteppingAction
{
public:
  virtual void BeginOfEventAction() {};
  virtual void EndOfEventAction() {};
};

#endif
