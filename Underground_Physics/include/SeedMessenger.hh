#ifndef SEED_MESSENGER_H
#define SEED_MESSENGER_H

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4Tokenizer.hh"

class SeedMessenger: public G4UImessenger {
  public:
  SeedMessenger() {
    directory = new G4UIdirectory("/seed/");
    directory->SetGuidance("Set random number seeds.");

    seedCmd = new G4UIcmdWithAString("/seed/setSeeds", this);
    seedCmd->SetGuidance("Initialize the random number generator with integer seed stream.");
    seedCmd->SetGuidance("Number of integers should be more than 1.");
    seedCmd->SetGuidance("Actual number of integers to be used depends on the individual random number engine.");
    #ifdef G4MULTITHREADED
    seedCmd->SetGuidance("This command sets the seeds for the master thread.");
    #endif
    seedCmd->SetParameterName("IntArray", false);
    seedCmd->AvailableForStates(G4State_PreInit, G4State_Idle, G4State_GeomClosed);
    seedCmd->SetToBeBroadcasted(false);
  }
  ~SeedMessenger() {
    delete directory;
    delete seedCmd;
  }
  virtual void SetNewValue(G4UIcommand* command, G4String newValue) {
    if( command==seedCmd ) {
      G4Tokenizer next(newValue);
      G4int idx=0;
      G4long seeds[100];
      G4String vl;
      while(!(vl=next()).empty()) {
        seeds[idx] = StoL(vl);
        idx++;
      }
      if(idx<2) {
        G4cerr << "/seed/setSeeds should have at least two integers. Command ignored." << G4endl;
      }
      else {
        seeds[idx] = 0;
        G4Random::setTheSeeds(seeds);
      }
    }
  }
  protected:
  inline G4long StoL(G4String str) const {
    long vl;
    const char* t = str;
    std::istringstream is(t);
    is >> vl;
    return vl;
  }
  private:
  G4UIdirectory* directory;
  G4UIcmdWithAString* seedCmd;
};
#endif
