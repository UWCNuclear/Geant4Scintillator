#include "OpPhotonPrimaryGeneratorAction.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"


OpPhotonPrimaryGeneratorAction::OpPhotonPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), 
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  //default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.*ns);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm ,0.*cm, -10.*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fParticleGun->SetParticleEnergy(511.*keV);
}


OpPhotonPrimaryGeneratorAction::~OpPhotonPrimaryGeneratorAction()
{
  delete fParticleGun;
}


void OpPhotonPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent);
}