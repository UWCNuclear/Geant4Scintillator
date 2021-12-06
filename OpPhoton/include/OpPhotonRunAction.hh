#ifndef OpPhotonRunAction_h
#define OpPhotonRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Timer;

class OpPhotonRunAction : public G4UserRunAction
{
public:
  OpPhotonRunAction();
  virtual ~OpPhotonRunAction();

public:
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  private:
    G4Timer* fTimer;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*OpPhotonRunAction_h*/
