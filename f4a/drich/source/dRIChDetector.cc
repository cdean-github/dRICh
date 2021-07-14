#include "dRIChDetector.h"
#include "ci_DRICH_Config.hh"
#include "g4dRIChOptics.hh"

#include <fun4all/Fun4AllBase.h>

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>

#include <G4Box.hh>
#include <G4Color.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>
#include <G4VisAttributes.hh>
#include <G4tgbVolumeMgr.hh>

#include <cmath>
#include <iostream>

class G4VSolid;
class PHCompositeNode;

/*
using std::cerr;
using std::cout;
using std::endl;
using std::string;
*/
// ---------------------------------------------------
dRIChDetector::dRIChDetector(PHG4Subsystem *subsys, PHCompositeNode *Node,
                             PHParameters *parameters, const std::string &dnam)
    : PHG4Detector(subsys, Node, dnam), m_Params(parameters)
    {}

// ---------------------------------------------------
int dRIChDetector::IsInDetector(G4VPhysicalVolume *volume) const 
{
  std::set<G4VPhysicalVolume *>::const_iterator iter =
      m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end()) 
  {
    return 1;
  }
  return 0;
}

// ---------------------------------------------------
void dRIChDetector::ConstructMe(G4LogicalVolume *logicWorld) 
{

  // get env vars
  char *drichHome = std::getenv("DRICH_HOME");
  if (drichHome == NULL) 
  {
    std::cerr << "[+] ERROR: DRICH_HOME env var not set" << std::endl;
    return;
  }

  // load model text file and configuration
  ci_DRICH_Config cfg;
  cfg.model_file = "/sphenix/user/cdean/ECCE/test/dRICh/f4a/drich/source/drich-g4model.txt";
  std::cout << "[+] MODEL TEXT FILE: " << cfg.model_file << std::endl;
  // - check existence
  std::ifstream mf(cfg.model_file.data());
  if (!mf.is_open()) 
  {
    std::cerr << "[+] ERROR in dRIChDetector: cannot find MODEL TEXT FILE" << std::endl;
    return;
  }
  mf.close();

  // graphics
  G4VisAttributes *vis = new G4VisAttributes(G4Color(0., 0., 0.9, 0.5));
  vis->SetForceSolid(true);
  vis->SetVisibility(true);
  vis->SetLineWidth(1);
  vis->SetForceSolid(true);
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  // logicWorld->SetVisAttributes(vis);

  // build detector by text file
  std::cout << "[+] read model text file" << std::endl;
  G4tgbVolumeMgr *volmgr = G4tgbVolumeMgr::GetInstance();
  volmgr->AddTextFile(cfg.model_file);
  std::cout << "[+] construct detector from text file" << std::endl;
  G4VPhysicalVolume *vesselPhysVol = volmgr->ReadAndConstructDetector();
  std::cout << "[+] detector summary" << std::endl;
  volmgr->DumpSummary();
  std::cout << "[+] detector G4Solid list" << std::endl;
  volmgr->DumpG4SolidList();

  // material optical properties (see shared header g4dRIChOptics.hh)
  // - aerogel
  auto aeroPO = new g4dRIChAerogel("ciDRICHaerogelMat");
  aeroPO->setOpticalParams(cfg.aerOptModel); // mode=3: use experimental data
  // - acrylic filter
  printf("[+] Acrylic Wavelength Threshold : %f nm\n", cfg.filter_thr / (nm));
  auto acryPO = new g4dRIChFilter("ciDRICHfilterMat");
  acryPO->setOpticalParams(cfg.filter_thr);
  // - gas radiator
  auto gasPO = new g4dRIChGas("ciDRICHgasMat");
  gasPO->setOpticalParams();
  // - photo sensors
  auto photoSensor = new g4dRIChPhotosensor("ciDRICHpsst");
  photoSensor->setOpticalParams("ciDRICH");
  // - mirror (simular to photosensor, but different params)
  auto mirror = new g4dRIChMirror("ciDRICHmirror");
  mirror->setOpticalParams("ciDRICH");

  // add to logical world
  logicWorld->AddDaughter(vesselPhysVol);

  // activate volumes, for hit readout
  this->ActivateVolumeTree(vesselPhysVol);

  return;
}

// ---------------------------------------------------
// recursively add detectors to active volume list, descending the tree
// from `volu`
// - use the "activation filter" to decide for which volumes to save hits
// - the petal number is added to `m_PetalMap`, which together
//   with the copy number, provides a unique ID for each photo sensor
void dRIChDetector::ActivateVolumeTree(G4VPhysicalVolume *volu, G4int petal) 
{

  // get objects
  G4String voluName = volu->GetName();
  G4int voluCopyNo = volu->GetCopyNo();
  G4LogicalVolume *logi = volu->GetLogicalVolume();

  // obtain petal number
  if (voluName.contains("petal"))
  {
    petal = voluCopyNo;
  }

  // activation filter: use this to decide which volumes to save
  // hits for, i.e., which volumes are "active"
  // TODO: need to decide what volume we want to be active
  G4bool activate = voluName.contains("psst") || voluName.contains("vessel");
  activate = true; // override
  if (activate) 
  {
    if (Verbosity() >= Fun4AllBase::VERBOSITY_SOME)  
    {
      std::cout << "[+] activate " << voluName << " petal " << petal << " copy "
           << voluCopyNo << std::endl;
    }
    m_PhysicalVolumesSet.insert(volu);
    m_PetalMap.insert(std::pair<G4VPhysicalVolume *, G4int>(volu, petal));
  }

  // loop over daughters
  G4int nd = logi->GetNoDaughters();
  for (int d = 0; d < nd; d++) 
  {
    this->ActivateVolumeTree(logi->GetDaughter(d), petal);
  }
}

// ---------------------------------------------------
// get petal number
int dRIChDetector::GetPetal(G4VPhysicalVolume *volu) 
{
  int petalNum;
  try 
  {
    petalNum = m_PetalMap.at(volu);
  } 
  catch (const std::out_of_range &ex) 
  {
    std::cerr << "ERROR in dRIChDetector: cannot find petal associated with volume"
         << std::endl;
    return -1;
  }
  return petalNum;
}

// ---------------------------------------------------
// get PSST number
int dRIChDetector::GetPSST(G4VPhysicalVolume *volu) 
{
  return volu->GetName().contains("psst") ? volu->GetCopyNo() : 0;
}

// ---------------------------------------------------
void dRIChDetector::Print(const std::string &what) const 
{
  std::cout << "dRICh Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME") 
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}
