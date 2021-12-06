#ifndef OpPhotonStackingAction_H
#define OpPhotonStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpPhotonStackingAction : public G4UserStackingAction
{
  public:
    OpPhotonStackingAction();
    virtual ~OpPhotonStackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();

  private:
    G4int fScintillationCounter;
    G4int fCerenkovCounter;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
