//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// DetectorConstruction header
// --------------------------------------------------------------

#ifndef DMXDetectorConstruction_h
#define DMXDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class G4UserLimits;

class DMXScintSD;
class DMXPmtSD;

class DMXDetectorMessenger;

class DMXDetectorConstruction : public G4VUserDetectorConstruction 
{
public:

  DMXDetectorConstruction();
  ~DMXDetectorConstruction();

public:

  G4VPhysicalVolume* Construct();
  void ConstructSDandField();

  void SetRoomEnergyCut(G4double);
  void SetEnergyCut(G4double);
  void SetTimeCut(G4double);
  void SetRoomTimeCut(G4double);
 
private:

  void DefineMaterials();

  G4UserLimits*    theUserLimitsForRoom; 
  G4UserLimits*    theUserLimitsForDetector; 
  //  G4UserLimits*    theUserLimitsForXenon; 

  G4double         theMaxTimeCuts;
  G4double         theMaxStepSize;
  G4double         theDetectorStepSize;
  G4double         theMinEkine;
  G4double         theRoomMinEkine;
  
  G4double         theRoomTimeCut;


#include "DMXDetectorMaterial.ihh"  // materials used

  G4double sourceZ;

  G4LogicalVolume*   world_log;        // pointers
  G4VPhysicalVolume* world_phys;  



  G4LogicalVolume*   lab_log;
  G4VPhysicalVolume* lab_phys;  


 
 




  /*
  G4VPhysicalVolume* vesseltop_phys1;
  G4VPhysicalVolume* vesselbottom_phys2;
  G4LogicalVolume*   GXe_log;
  G4VPhysicalVolume* GXe_phys;  
  G4LogicalVolume*   gaslag_log;
  G4VPhysicalVolume* gaslag_phys;  
  */
  G4LogicalVolume*   LXe_log; 
  G4VPhysicalVolume* LXe_phys; 
  /*
  G4LogicalVolume*   liqLag_log; 
  G4VPhysicalVolume* liqLag_phys;  
  */
  //G4LogicalVolume*   pmt_log;   
  //G4VPhysicalVolume* pmt_phys; 
  G4LogicalVolume*   phcath_log;


  G4Cache<DMXScintSD*> LXeSD; //pointer to sensitive detectors
  G4Cache<DMXPmtSD*> pmtSD;

  // pointer to the Detector Messenger:
  DMXDetectorMessenger*  detectorMessenger;

};

#endif

