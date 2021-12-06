#ifndef OpPhotonDetectorConstruction_h
#define OpPhotonDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpPhotonDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    OpPhotonDetectorConstruction();
    virtual ~OpPhotonDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();

    G4double CalculateReflectivity(G4double, G4double);

  private:
    G4double fExpHall_halfX;
    G4double fExpHall_halfY;
    G4double fExpHall_halfZ;

    G4double fStick_halfX;
    G4double fStick_halfY;
    G4double fStick_halfZ;

    G4double fWrapper_halfThk;
    G4double fWrapperElem_halfZ;
    G4int fNbWrapperDiv;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*OpPhotonDetectorConstruction_h*/
