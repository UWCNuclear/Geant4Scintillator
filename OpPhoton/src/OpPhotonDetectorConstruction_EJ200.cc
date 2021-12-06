#include "OpPhotonDetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"

#include <vector> 

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Polyhedra.hh"

OpPhotonDetectorConstruction::OpPhotonDetectorConstruction()
 : G4VUserDetectorConstruction()
{
  // Half the dimensions
  fExpHall_halfX = fExpHall_halfY = fExpHall_halfZ = 1.*m;
  
  fStick_halfX = fStick_halfY = 5*mm; 
// fStick_halfZ = 2.5*cm;
  fStick_halfZ = 5.*cm;
//  fStick_halfZ = 7.5*cm;
  
  fWrapper_halfThk = 3.*mm;  
//  fWrapper_halfThk = 1.5*mm;

  // The wrapper is divided in 1 mm elements that will each be assigned a 
  // different reflectivity 
  fWrapperElem_halfZ = .5*mm; 
  fNbWrapperDiv = fStick_halfZ / fWrapperElem_halfZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpPhotonDetectorConstruction::~OpPhotonDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* OpPhotonDetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

  // Elements for the EJ-200 and the PTFE
  G4Element* carbon = new G4Element("Carbon", "C", 6., 12.011*g/mole);
  G4Element* hydrogen = new G4Element("Hydrogen", "H", 1., 1.00794*g/mole);
  G4Element* fluorine = new G4Element("Fluorine", "F", 9., 19.00*g/mole);

  // Plastic scintillators: EJ-200 of Eljen Technology
  G4Material* EJ200 = new G4Material("EJ200", 1.023*g/cm3, 2);
  EJ200->AddElement(carbon, .915);
  EJ200->AddElement(hydrogen, .085);

  // Polytetrafluoroethylene (PTFE) wrapping
  G4Material* PTFE = new G4Material("PTFE", 2.2*g/cm3, 2);
  PTFE->AddElement(carbon, 2);
  PTFE->AddElement(fluorine, 4);

  // ------------ Generate & Add Material Properties Table ------------
  //
  // EJ-200 properties come from Eljen Technology document 
  // "General purpose plastic scintillator EJ-200, EJ-204, EJ-208, EJ-212"

  G4double photonEnergy[] = {2.48*eV, 2.49*eV, 2.50*eV, 2.52*eV,
                             2.54*eV, 2.56*eV, 2.57*eV, 2.59*eV,
                             2.61*eV, 2.63*eV, 2.65*eV, 2.67*eV,
                             2.68*eV, 2.69*eV, 2.71*eV, 2.73*eV,
                             2.75*eV, 2.76*eV, 2.78*eV, 2.79*eV,
                             2.81*eV, 2.82*eV, 2.84*eV, 2.86*eV,
                             2.88*eV, 2.90*eV, 2.91*eV, 2.92*eV,
                             2.94*eV, 2.95*eV, 2.96*eV, 2.97*eV,
                             2.98*eV, 2.99*eV, 3.00*eV, 3.02*eV,
                             3.03*eV, 3.04*eV, 3.06*eV, 3.07*eV,
                             3.09*eV, 3.10*eV};

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  G4double refractiveIndex1[] = {1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30, 1.30, 1.30,
                                 1.30, 1.30};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] = {3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m, 3.8*m, 3.8*m,
                           3.8*m, 3.8*m};

  assert(sizeof(absorption) == sizeof(photonEnergy));

  G4double scintilFast[] = {0.06, 0.06, 0.08, 0.09,
                            0.11, 0.14, 0.16, 0.19,
                            0.23, 0.28, 0.33, 0.39,
                            0.41, 0.43, 0.46, 0.49,
                            0.53, 0.56, 0.60, 0.65,
                            0.71, 0.76, 0.84, 0.89,
                            0.94, 0.99, 1.00, 1.00,
                            0.95, 0.90, 0.81, 0.73,
                            0.61, 0.52, 0.40, 0.32,
                            0.25, 0.18, 0.10, 0.05,
                            0.02, 0.01};

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1, nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("FASTCOMPONENT", photonEnergy, scintilFast, nEntries)
        ->SetSpline(true);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE", 0.);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
  myMPT1->AddConstProperty("YIELDRATIO", 1.);

  G4cout << "Water G4MaterialPropertiesTable" << G4endl;
  myMPT1->DumpTable();

  EJ200->SetMaterialPropertiesTable(myMPT1);

  // Birks constant for the EJ-200
  // Data from Tajudin et al. Appl. Radiat. Isot. 159 (2020) 109086
  EJ200->GetIonisation()->SetBirksConstant(0.156*mm/MeV);

  // Air
  G4double refractiveIndex2[] = {1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00,
                                 1.00, 1.00, 1.00, 1.00};

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntries);

  G4cout << "Air G4MaterialPropertiesTable" << G4endl;
  myMPT2->DumpTable();

  air->SetMaterialPropertiesTable(myMPT2);

  // Setup geometry

  // The experimental Hall
  G4Box* expHall_solid = new G4Box("ExpHall", fExpHall_halfX, fExpHall_halfY,
                                   fExpHall_halfZ);

  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_solid, air,
                                                     "ExpHall", 0, 0, 0);

  G4VPhysicalVolume* expHall_phys = new G4PVPlacement(0, G4ThreeVector(), 
                                                      expHall_log, "ExpHall", 
                                                      0, false, 0);
  // Plastic scintillator stick

// area of 4 faces = 2x*2y
// area of 1 face = pi*r^2
// area of circle segment = (theta*pi/360 - sin(theta)/2)*r^2
// area of 3 faces equilateral = 3^0.5/4*x^2 = pi*r^2 - 3*(theta*pi/360 - sin(theta)/2)*r^2 = r^2(pi - 3*(theta*pi/360 - sin(theta)/2))
// area of 6 faces = 3*3^0.5/2*x^2 = 3/2*x*h = pi*r^2 - 6*(theta*pi/360 - sin(theta)/2)*r^2 = r^2(pi - 6*(theta*pi/360 - sin(theta)/2))
// area of 3 faces right = 0.5*a*b = 0.5*a^2/tan(theta)   using tan(theta) = a/b


G4double area4 = 2 * fStick_halfX * 2 * fStick_halfY;

// SQUARE STICK
//  G4Box* stick_solid = new G4Box("Stick", fStick_halfX, fStick_halfY, fStick_halfZ);
//G4double xx = 0.;
//G4double yy = 0.;


// CYLINDER ROD
/*G4double cylrad = sqrt(area4 / pi);
  G4VSolid* stick_solid = new G4Tubs("Stick", 0., cylrad, fStick_halfZ, 0., 2*pi);
G4double xx = 0.;
G4double yy = 0.;
*/

// PRISMS
// From https://geant4-forum.web.cern.ch/t/need-advice-how-to-create-hexagonal-prism/5302/5
/*const G4int nsect = 3; // 3, 4, 6
std::vector<G4TwoVector> polygon(nsect);
G4double ang = twopi / nsect;
G4double dR = sqrt(area4 / (pi - nsect*(ang - sin(ang))/2));
//G4double dR = fStick_halfX;

for (G4int i = 0; i < nsect; ++i)
{
	G4double phi = i * ang;
	G4double cosphi = std::cos(phi);
	G4double sinphi = std::sin(phi);
	polygon[i].set(dR * cosphi, dR * sinphi);
}

G4TwoVector offsetA(0., 0.), offsetB(0., 0.);
G4double scaleA = 1., scaleB = 1.;
G4VSolid* stick_solid = new G4ExtrudedSolid("Extruded", polygon, fStick_halfZ, offsetA, scaleA, offsetB, scaleB);
G4double xx = 0.;
G4double yy = 0.;
*/

// NOT EQUILATERAL TRIANGLES
G4double theta = 45.*deg; // 45, 30 degrees
/*G4RotationMatrix* rm = new G4RotationMatrix();
rm->rotateZ(theta);
G4double side3y = sqrt(area4 * 2 * tan(theta))/2.;
G4double side3x = side3y / tan(theta);
G4VSolid* box1 = new G4Box("Box1", side3x, side3y, fStick_halfZ);
G4VSolid* box2 = new G4Box("Box2", side3x*2., side3y, fStick_halfZ);
G4VSolid* stick_solid = new G4SubtractionSolid("Box1-Box2", box1, box2, rm, G4ThreeVector(side3y*cos(theta),side3y*sin(theta),0.));
G4double xx = side3y*cos(theta)/sqrt(2.);
G4double yy = side3y*sin(theta)/sqrt(2.);
*/

// Polyhedra  // All prisms
G4double side3x = sqrt(area4 * 2 * tan(theta))/2.;
  const G4int nrz2 = 4;
  G4double zz2[nrz2] = {-fStick_halfZ, -fStick_halfZ, fStick_halfZ, fStick_halfZ}, rr2[nrz2] = { 0., side3x, side3x, 0. };
  G4VSolid* stick_solid = new G4Polyhedra("Stick", 0., 2.*pi, 8, nrz2, rr2, zz2);
G4double xx = 0.;
G4double yy = 0.;


// Polyhedra  // Right angle triangle
/*  const G4int nrz2 = 4;
  G4double zz2[nrz2] = {-fStick_halfZ, -fStick_halfZ, fStick_halfZ, fStick_halfZ}, rr2[nrz2] = { 0., side3x, side3x, 0. }; 
  G4VSolid* stick_solid = new G4Polyhedra("Stick", 0., pi, 2, nrz2, rr2, zz2);
G4double xx = side3x*cos(theta)/sqrt(2.);
G4double yy = side3x*cos(theta)/sqrt(2.);
*/

  G4LogicalVolume* stick_log = new G4LogicalVolume(stick_solid, EJ200, "Stick",
                                                   0, 0, 0);

// Translation of NOT EQUILATERAL TRIANGLES
  G4VPhysicalVolume* stick_phys =
   new G4PVPlacement(0, G4ThreeVector(xx, yy, fStick_halfZ), stick_log, "Stick", expHall_log, false, 0, true);  

  
  // Stick ends to collect the optical photons. So far these are filled with air.
  // For more realisitc simulations, they can be filled with SiPM material.
  G4double stickEnd_halfZ = 1.*cm;

  G4Box* stickEnd_solid = new G4Box("StickEnd", fStick_halfX, fStick_halfY, 
                                    stickEnd_halfZ);

  G4LogicalVolume* stickEnd_log = new G4LogicalVolume(stickEnd_solid, air,
                                                      "StickEnd", 0, 0, 0);

  G4VPhysicalVolume* stickFront_phys =
   new G4PVPlacement(0, G4ThreeVector(0., 0., -stickEnd_halfZ),
                     stickEnd_log, "StickFrontEnd", expHall_log, false, 0,
                     true);

  G4VPhysicalVolume* stickBack_phys =
   new G4PVPlacement(0, 
                     G4ThreeVector(0., 0., 2. * fStick_halfZ + stickEnd_halfZ),
                     stickEnd_log, "StickBackEnd", expHall_log, false, 0, true);

  // PTFE wrapper

  // The wrapper is divided in 1 mm long elements that will each be assigned
  // a specific reflectivity.
  // Wrapper element = small box with an extrusion of the shape of the stick
  // cross section.

  // Unextruded wrapper
  G4Box* unextrWrapper_solid = new G4Box("UnextrWrapper",
                                         fWrapper_halfThk + fStick_halfX,
                                         fWrapper_halfThk + fStick_halfY,
                                         fWrapperElem_halfZ);

  // Stick shaped extrusion in the wrapper
  G4VSolid* extrWrapper_solid =
   new G4SubtractionSolid("ExtrWrapper", unextrWrapper_solid, stick_solid, 0,
                          G4ThreeVector(xx, yy, 0.));

  G4LogicalVolume* wrapper_log = new G4LogicalVolume(extrWrapper_solid, PTFE,
                                                     "Wrapper", 0, 0, 0);

  std::vector<G4VPhysicalVolume*> wrapper_phys;
  G4double posZ;
  std::vector<G4double> wrapperElemPosZ;
 
  // To place the wrapper elements
  for (G4int i = 0; i < fNbWrapperDiv; i++)
  {
    posZ = (2 * i + 1) * fWrapperElem_halfZ;
    
    // Store the positions of the wrapper elements. These is used later to
    // calulate their reflectivity.
    wrapperElemPosZ.push_back(posZ);

    wrapper_phys.push_back(
      new G4PVPlacement(0, G4ThreeVector(0., 0., wrapperElemPosZ.at(i)),
                        wrapper_log, "Wrapper", expHall_log, false, 0, true));
  }

  // ------------- Surfaces --------------

  G4double ePhoton[] = {2.48*eV, 3.10*eV};
  const G4int num = sizeof(ePhoton)/sizeof(G4double);

  // EJ-200 -> PTFE
  std::vector<G4OpticalSurface*> scintWrap;
  std::vector<G4MaterialPropertiesTable*> scintWrapProperty;
  
  for (G4int i = 0; i < fNbWrapperDiv; i++)
  {
    char opSurfacName[64];
    sprintf(opSurfacName, "ScintWrap%i", i);

    scintWrap.push_back(new G4OpticalSurface(opSurfacName));

    new G4LogicalBorderSurface(opSurfacName, stick_phys, wrapper_phys.at(i),
                               scintWrap.at(i));

    scintWrap.at(i)->SetType(dielectric_metal);
    scintWrap.at(i)->SetFinish(polished);
    scintWrap.at(i)->SetModel(glisur);

    G4double reflectivityAtZ = CalculateReflectivity(wrapperElemPosZ.at(i),
                                                     2 * fStick_halfZ);
    
    G4cout << "Stick depth (cm) = " << wrapperElemPosZ.at(i) / cm
           << ", reflectivity = " << reflectivityAtZ << G4endl;

    G4double reflectivity[] = {reflectivityAtZ, reflectivityAtZ};
    assert(sizeof(reflectivity) == sizeof(ePhoton));
    G4double efficiency[] = {1., 1.};
    assert(sizeof(efficiency) == sizeof(ePhoton));

    scintWrapProperty.push_back(new G4MaterialPropertiesTable());

    scintWrapProperty.at(i)->AddProperty("REFLECTIVITY", ePhoton, reflectivity, num);
    scintWrapProperty.at(i)->AddProperty("EFFICIENCY", ePhoton, efficiency, num);
    scintWrap.at(i)->SetMaterialPropertiesTable(scintWrapProperty.at(i));
  }

  // EJ-200 -> Air
  G4OpticalSurface* stickToAir = new G4OpticalSurface("StickAir");

  new G4LogicalBorderSurface("StickToFront", stick_phys, stickFront_phys,
                             stickToAir);

  new G4LogicalBorderSurface("StickToBack", stick_phys, stickBack_phys,
                             stickToAir);

  stickToAir->SetType(dielectric_dielectric);
  stickToAir->SetFinish(polished);
  stickToAir->SetModel(glisur);

  G4double reflectivityStickToAir[] = {0., 0.};
  assert(sizeof(reflectivityStickToAir) == sizeof(ePhoton));
  G4double efficiencyStickToAir[] = {1., 1.};
  assert(sizeof(efficiencyStickToAir) == sizeof(ePhoton));

  G4MaterialPropertiesTable* stickToAirProperty =
   new G4MaterialPropertiesTable();

  stickToAirProperty->AddProperty("REFLECTIVITY", ePhoton,
                                  reflectivityStickToAir, num);

  stickToAirProperty->AddProperty("EFFICIENCY", ePhoton,
                                  efficiencyStickToAir, num);
  
  stickToAir->SetMaterialPropertiesTable(stickToAirProperty);

  // For the setup visualisation
  expHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  stickEnd_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  //stick_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes* stick_vis = new G4VisAttributes(G4Colour(1., 1., 0., .5));
  stick_vis->SetDaughtersInvisible(false);
  stick_vis->SetForceSolid(true);
  stick_log->SetVisAttributes(stick_vis);

  //wrapper_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes* wrapper_vis = new G4VisAttributes(G4Colour(.7, .7, .7, .3));
  wrapper_vis->SetDaughtersInvisible(false);
  wrapper_vis->SetForceSolid(true);
  wrapper_log->SetVisAttributes(wrapper_vis);

  return expHall_phys;
}

G4double OpPhotonDetectorConstruction::CalculateReflectivity(
                                                      G4double wrapperElemPosZ,
                                                      G4double stickLength)
{
  G4double reflectivity_min = .975;
  G4double reflectivity_max = .975;
//  G4double reflectivity_max = .975;

  G4double reflectivity = 0.;

  // From the stick entrance to the middle, the reflectivity inceases linearly
  // from 0.9 to 0.975
  if (wrapperElemPosZ <= stickLength / 2.)
   reflectivity =
    2. * (reflectivity_max - reflectivity_min) * wrapperElemPosZ / stickLength
    + reflectivity_min;

  // From the stick middle to the exit, the reflectivity deceases linearly from
  // 0.975 to 0.9
  if (wrapperElemPosZ > stickLength / 2.)
   reflectivity =
    2. * (- reflectivity_max + reflectivity_min)
    * (wrapperElemPosZ - stickLength / 2.) / stickLength
    + reflectivity_max;

  return reflectivity;
}
