//-----------------------
//Simple Central Detector: simple resolution parameterization
//Perfect efficiency and PID.
//Acceptance from [-4,4] in eta
//
//Simple forward Detectors for Hadrons
//Perfect efficiency+PID+Resolution (for now)
//Acceptance from [4.4,6.5] in eta (~3-25mrad)
//-----------------------

//To Build with ROOT6:
//>>root
//>>detector_full.C+
//-----------------------

#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/Acceptance.h"
#include "eicsmear/smear/Device.h"
#include "eicsmear/smear/Detector.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/PerfectID.h"

//Convert pseudorapidity (eta) to polar angle (theta) in radians.
//Make use of TLorentzVector to do eta-to-theta conversion.
double etaToTheta(const double eta) {
  TLorentzVector v;
  v.SetPtEtaPhiM(1., eta, 0., 0.);
  return v.Theta();
}

//Convert and angle in degrees to one in radians.
double degreesToRadians(double degrees) {
  return degrees / 180. * TMath::Pi();
}

//Build Detector
Smear::Detector BuildDetector(){

  //Create Devices
  Smear::Device energy(Smear::kE,"0.2 * TMath::Sqrt(E)");
  Smear::Device momentum(Smear::kP,"0.01 * P");
  Smear::Device theta(Smear::kTheta,"0.01 / ( P * TMath::Sqrt(TMath::Sin(theta)) )");
  Smear::Device phi(Smear::kPhi,"0.01");
  
  //Central Detector Acceptance
  Smear::Acceptance::Zone central(etaToTheta(4.),etaToTheta(-4.));

  energy.Accept.AddZone(central);
  momentum.Accept.AddZone(central);
  theta.Accept.AddZone(central);
  phi.Accept.AddZone(central);

  //Forward Detector
  Smear::Acceptance::Zone forward(etaToTheta(6.5),etaToTheta(4.4));

  Smear::Device forward_p(Smear::kP,"0");
  Smear::Device forward_E(Smear::kE,"0");
  Smear::Device forward_theta(Smear::kTheta,"0");
  Smear::Device forward_phi(Smear::kPhi,"0");

  forward_p.Accept.SetGenre(Smear::kHadronic);
  forward_E.Accept.SetGenre(Smear::kHadronic);
  forward_phi.Accept.SetGenre(Smear::kHadronic);
  forward_theta.Accept.SetGenre(Smear::kHadronic);

  forward_p.Accept.AddZone(forward);
  forward_E.Accept.AddZone(forward);
  forward_theta.Accept.AddZone(forward);
  forward_phi.Accept.AddZone(forward); 

  //PID performance is unparameterised as of now
  Smear::PerfectID pid;
  pid.Accept.AddZone(central);
  pid.Accept.AddZone(forward);
  
  //Create the detector and add devices
  Smear::Detector det;
  det.AddDevice(energy);
  det.AddDevice(momentum);
  det.AddDevice(theta);
  det.AddDevice(phi);
  det.AddDevice(pid);
  det.AddDevice(forward_p);
  det.AddDevice(forward_E);
  det.AddDevice(forward_phi);
  det.AddDevice(forward_theta);

  det.SetEventKinematicsCalculator("NM JB DA");

  return det;

}
