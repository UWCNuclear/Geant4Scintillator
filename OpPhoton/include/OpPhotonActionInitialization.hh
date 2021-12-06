#ifndef OpPhotonActionInitialization_h
#define OpPhotonActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class OpPhotonActionInitialization : public G4VUserActionInitialization
{
  public:
    OpPhotonActionInitialization();
    virtual ~OpPhotonActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif
