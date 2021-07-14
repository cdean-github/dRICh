//-----------------------------------------------
// dRICh event tree module
//-----------------------------------------------

#ifndef DRICHTREE_H
#define DRICHTREE_H

#include <Geant4/G4String.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Rtypes.h>
#include <TString.h>
#include <fun4all/SubsysReco.h>

// class Fun4AllHistoManager; //---
class PHCompositeNode;
class TFile;
class TTree;
// class TH1; //---
// class TH2; //---

class dRIChTree : public SubsysReco 
{
  public:
    // constructor
    dRIChTree(const std::string &name = "dRIChTree",
              const std::string &fname = "dRIChTree.root");

    // destructor
    virtual ~dRIChTree();

    // SubsysReco processing methods:
    int Init(PHCompositeNode *);
    int process_event(PHCompositeNode *);
    int End(PHCompositeNode *);

  private:
    std::string m_outfileN; // output file name
    TFile *m_outfile;       // output file
    TTree *m_tree;          // the dRICh tree

    // TH1 *m_histo; //---
    // Fun4AllHistoManager *m_hm; // histogram manager //---

    // data accessors
    void getHits(PHCompositeNode *topNode);
    // void getHEPMCTruth(PHCompositeNode *topNode); //---

    void vectorToArray(G4ThreeVector vec, Double_t *arr);
    void initTrees();
    // void resetVars(); //---

    //------------------------------
    // tree variables
    //------------------------------

    Int_t evnum;
    Int_t trackID;
    char hitType[128];
    char hitSubtype[128];
    Int_t petal, psst, pdg;
    char particleName[128];
    char process[256];
    Int_t parentID;
    Double_t hitPos[3]; // [xyz]
    Double_t hitP[3];
    Double_t hitPdir[3];
    Double_t hitVtxPos[3];
    Double_t hitVtxPdir[3];
    Double_t deltaT;
    Double_t edep;
};

#endif // DRICHTREE_H
