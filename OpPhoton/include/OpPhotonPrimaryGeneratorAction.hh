#ifndef OpPhotonPrimaryGeneratorAction_h
#define OpPhotonPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class OpPhotonPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    OpPhotonPrimaryGeneratorAction();
    virtual ~OpPhotonPrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun* fParticleGun;
};

#endif /*OpPhotonPrimaryGeneratorAction_h*/
