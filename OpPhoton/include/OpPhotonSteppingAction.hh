#ifndef AnnihilationPhotonsSteppingAction_hh
#define AnnihilationPhotonsSteppingAction_hh

#include "OpPhotonVSteppingAction.hh"

#include "G4ThreeVector.hh"

#include <vector> 

class OpPhotonSteppingAction : public OpPhotonVSteppingAction
{
  public:
    OpPhotonSteppingAction();
    virtual void BeginOfEventAction();
    virtual void UserSteppingAction(const G4Step*);
    virtual void EndOfEventAction();

  private:

    // These are used to remember quantities from call to call of 
    // UserSteppingAction
    G4double fComptonE;
    G4double fComptonDepth;
    G4double fComptonTime;

    G4int fNbFrontPhot;
    std::vector<G4double> fV_timeFront;
    G4int fNbBackPhot;
    std::vector<G4double> fV_timeBack;
    std::vector<G4int> fV_trackID;
};

#endif
