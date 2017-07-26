//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jul  1 14:41:11 2017 by ROOT version 5.34/34
// from TTree U102/Envent
// found on file: UBTF_1000.root
//////////////////////////////////////////////////////////

#ifndef UBTF_SelectorAnalysis_h
#define UBTF_SelectorAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <fstream>
#include <TGraph.h>
#include <iostream>
#include <TProfile.h>
#include "TF1.h"


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class UBTF_SelectorAnalysis : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   std::vector<TH1F*> fhVec;
//Lecce's Histograms
   TH1F           *hM2noCut;
   TH1F           *hM2leadingECALcl;
   TH1F           *hM2subleadECALcl;
   TH1F           *hM2subsublECALcl;
   TH1F           *hM2noCutAt0;
   TH1F           *hMnoCut;
   TH1F           *hMleadingECALcl;

   TH1F           *hM2clusterCut;
   TH1F           *hMclusterCut;

   TH1F           *energyCut;
   TH1F           *radiusCut;
   TH1F           *TimeCut;
   TH1F           *energyPositronPhotonCut;
   TH1F           *DeltaTime;
   TH1F           *TimeError;
   TH1F           *DeltaPVetoTrEne;
   TH1F           *PVetoTrEne1Error;
   TH1F           *controllo;
   TH1F           *ETotBremEv;
   TH1F           *ETotBremEv1;
   TH2F           *pVetoEClvsFinger;
   TH2F           *pVetoEClvsEphoton;
   TH2F           *Ephotonvsfinger;
   TH2F           *NfingervsPVetoTrEne;
   TH2F           *NfingervsPVetoTrEneBremm;
   TH2F           *TPhotonVsTVeto;
   TH2F           *PVetoFingVsTveto;
   TH2F           *xClustervseCluster;
   TH2F           *yClustervseCluster;
   TH2F           *xClustervsm2Cluster;
   TH2F           *yClustervsm2Cluster;
   TH2F           *eClustervsm2Cluster;
   TH2F           *xClustervsyCluster;
   TProfile       *TPhotonVsTVetoCut;
   TProfile       *PVetoTrEnevsNfingerProfile;

   // event counters
   int            fNProcessedEvents;
   int            fNEvAtLeastOneCluster;
   int            fNEv2Clusters;
   int            fNEv3Clusters;
   int            fNEv4Clusters;
   int            fNEvOneCluster;
   int            fNEvAtLeastOneClusterAboveThr;
   int            fNEv2ClustersAboveThr;
   int            fNEv3ClustersAboveThr;
   int            fNEv4ClustersAboveThr;
   int            fNEvOneClusterAboveThr;
   int            fNEvEphotonPass;
   int            fNEvRadiusPass;
   int            fNEvFingerPass;
   int            fNEvTimePass;
   int            NBremmEvT;
   int            NBremmEv;
   int            NPassBremm;

   int            fAprimeMass;


   Double_t       photonE;
   Double_t       radius;
   Double_t       Emin;
   Double_t       Emax;
   Double_t       thethamax;
   Double_t       thethamin;
   double         Rmin;
   double         Rmax;
   ofstream       out1;
   double         NEvFingerPassVector[20];
   double         Intercetta;
   double         IntercettaError;
   double         CoeffAngolare;
   double         CoeffAngolareError;
   int            sogliaCut;
   int            preSoglia;


   double         EvsAngle(double theta, double mA );
   void           SaveTrainingData(double m2Aprime, double ECluster, double XCluster, double YCluster, int go=0);

   bool            save;
   bool            m_getCorrelations;
   Double_t        fECALthresholdE;
   bool            m_useECALthreshold;



   // Declaration of leaf types
   Int_t           Nevent;
   Double_t        ETot;
   Double_t        IDProc;
   Double_t        PBeam;
   Double_t        PPrim;
   Double_t        XBeam;
   Double_t        YBeam;
   Int_t           NClusters;
   Int_t           NTracks;
   Int_t           NHEPVetoTracks;
   Int_t           NPVetoTracks;
   Int_t           NEVetoTracks;
   Int_t           NSAC;
   Int_t           NCal;
   Int_t           NLAV;
   Int_t           NTarget;
   Double_t        ESAC[100];
   Double_t        TSAC[100];
   Double_t        PTypeSAC[100];
   Double_t        XSAC[100];
   Double_t        YSAC[100];
   Int_t           SACCh[100];
   Double_t        EPartCal[20];
   Double_t        TPartCal[20];
   Int_t           PTypePartCal[20];
   Double_t        XPartCal[20];
   Double_t        YPartCal[20];
   Double_t        ECluster[20];
   Double_t        QCluster[20];
   Double_t        XCluster[20];
   Double_t        YCluster[20];
   Double_t        ThCluster[20];
   Double_t        M2Cluster[20];
   Double_t        TCluster[20];
   Double_t        NClusCells[20];
   Double_t        ETarget;
   Double_t        TTarget;
   Double_t        XTarget;
   Double_t        YTarget;
   Double_t        HEPVetoTrEne[100];
   Int_t           HEPVetoNFing[100];
   Double_t        HEPVetoTrTime[100];
   Double_t        HEPVetoFingE[100];
   Double_t        HEPVetoX[100];
   Double_t        HEPVetoY[100];
   Int_t           HEPVetoClIndex[100];
   Double_t        HEPVetoECl[100][10];
   Double_t        HEPVetoTimeCl[100][10];
   Double_t        PVetoTrEne[100];
   Int_t           PVetoNFing[100];
   Double_t        PVetoTrTime[100];
   Double_t        PVetoFingE[100];
   Double_t        PVetoX[100];
   Double_t        PVetoY[100];
   Double_t        PVetoBarE[100];
   Double_t        PVetoBarT[100];
   Int_t           PVetoClIndex[100];
   Double_t        PVetoECl[100][10];
   Double_t        PVetoTimeCl[100][10];
   Double_t        EVetoTrEne[100];
   Int_t           EVetoNFing[100];
   Double_t        EVetoTrTime[100];
   Double_t        EVetoFingE[100];
   Double_t        EVetoX[100];
   Double_t        EVetoY[100];
   Int_t           EVetoClIndex[100];
   Double_t        EVetoECl[100][10];
   Double_t        EVetoTimeCl[100][10];

   // List of branches
   TBranch        *b_Nevent;   //!
   TBranch        *b_ETot;   //!
   TBranch        *b_IDProc;   //!
   TBranch        *b_PBeam;   //!
   TBranch        *b_PPrim;   //!
   TBranch        *b_XBeam;   //!
   TBranch        *b_YBeam;   //!
   TBranch        *b_NClusters;   //!
   TBranch        *b_NTracks;   //!
   TBranch        *b_NHEPVetoTracks;   //!
   TBranch        *b_NPVetoTracks;   //!
   TBranch        *b_NEVetoTracks;   //!
   TBranch        *b_NSAC;   //!
   TBranch        *b_NCal;   //!
   TBranch        *b_NLAV;   //!
   TBranch        *b_NTarget;   //!
   TBranch        *b_ESAC;   //!
   TBranch        *b_TSAC;   //!
   TBranch        *b_PTypeSAC;   //!
   TBranch        *b_XSAC;   //!
   TBranch        *b_YSAC;   //!
   TBranch        *b_SACCh;   //!
   TBranch        *b_CalE;   //!
   TBranch        *b_CalT;   //!
   TBranch        *b_CalPType;   //!
   TBranch        *b_CalX;   //!
   TBranch        *b_CalY;   //!
   TBranch        *b_ECluster;   //!
   TBranch        *b_QCluster;   //!
   TBranch        *b_XCluster;   //!
   TBranch        *b_YCluster;   //!
   TBranch        *b_ThCluster;   //!
   TBranch        *b_M2Cluster;   //!
   TBranch        *b_TCluster;   //!
   TBranch        *b_NClusCells;   //!
   TBranch        *b_ETarget;   //!
   TBranch        *b_TTarget;   //!
   TBranch        *b_XTarget;   //!
   TBranch        *b_YTarget;   //!
   TBranch        *b_NTHEPVetoTrkEne;   //!
   TBranch        *b_NTHEPVetoTrkFinger;   //!
   TBranch        *b_NTHEPVetoTrkTime;   //!
   TBranch        *b_NTHEPVetoFingE;   //!
   TBranch        *b_NTHEPVetoX;   //!
   TBranch        *b_NTHEPVetoY;   //!
   TBranch        *b_NTHEPVetoClIndex;   //!
   TBranch        *b_NTHEPVetoECl;   //!
   TBranch        *b_NTHEPVetoTimeCl;   //!
   TBranch        *b_NTPVetoTrkEne;   //!
   TBranch        *b_NTPVetoTrkFinger;   //!
   TBranch        *b_NTPVetoTrkTime;   //!
   TBranch        *b_NTPVetoFingE;   //!
   TBranch        *b_NTPVetoX;   //!
   TBranch        *b_NTPVetoY;   //!
   TBranch        *b_NTPVetoBarEnergy;   //!
   TBranch        *b_NTPVetoBarTime;   //!
   TBranch        *b_NTPVetoClIndex;   //!
   TBranch        *b_NTPVetoECl;   //!
   TBranch        *b_NTPVetoTimeCl;   //!
   TBranch        *b_NTEVetoTrkEne;   //!
   TBranch        *b_NTEVetoTrkFinger;   //!
   TBranch        *b_NTEVetoTrkTime;   //!
   TBranch        *b_NTEVetoFingE;   //!
   TBranch        *b_NTEVetoX;   //!
   TBranch        *b_NTEVetoY;   //!
   TBranch        *b_NTEVetoClIndex;   //!
   TBranch        *b_NTEVetoECl;   //!
   TBranch        *b_NTEVetoTimeCl;   //!

   UBTF_SelectorAnalysis(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~UBTF_SelectorAnalysis() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   void            PrintCounters();
   void 	   setUbosonMass(double mass){fAprimeMass = mass;}	


   ClassDef(UBTF_SelectorAnalysis,0);
};

#endif

#ifdef UBTF_SelectorAnalysis_cxx
void UBTF_SelectorAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Nevent", &Nevent, &b_Nevent);
   fChain->SetBranchAddress("ETot", &ETot, &b_ETot);
   fChain->SetBranchAddress("IDProc", &IDProc, &b_IDProc);
   fChain->SetBranchAddress("PBeam", &PBeam, &b_PBeam);
   fChain->SetBranchAddress("PPrim", &PPrim, &b_PPrim);
   fChain->SetBranchAddress("XBeam", &XBeam, &b_XBeam);
   fChain->SetBranchAddress("YBeam", &YBeam, &b_YBeam);
   fChain->SetBranchAddress("NClusters", &NClusters, &b_NClusters);
   fChain->SetBranchAddress("NTracks", &NTracks, &b_NTracks);
   fChain->SetBranchAddress("NHEPVetoTracks", &NHEPVetoTracks, &b_NHEPVetoTracks);
   fChain->SetBranchAddress("NPVetoTracks", &NPVetoTracks, &b_NPVetoTracks);
   fChain->SetBranchAddress("NEVetoTracks", &NEVetoTracks, &b_NEVetoTracks);
   fChain->SetBranchAddress("NSAC", &NSAC, &b_NSAC);
   fChain->SetBranchAddress("NCal", &NCal, &b_NCal);
   fChain->SetBranchAddress("NLAV", &NLAV, &b_NLAV);
   fChain->SetBranchAddress("NTarget", &NTarget, &b_NTarget);
   fChain->SetBranchAddress("ESAC", ESAC, &b_ESAC);
   fChain->SetBranchAddress("TSAC", TSAC, &b_TSAC);
   fChain->SetBranchAddress("PTypeSAC", PTypeSAC, &b_PTypeSAC);
   fChain->SetBranchAddress("XSAC", XSAC, &b_XSAC);
   fChain->SetBranchAddress("YSAC", YSAC, &b_YSAC);
   fChain->SetBranchAddress("SACCh", SACCh, &b_SACCh);
   fChain->SetBranchAddress("EPartCal", EPartCal, &b_CalE);
   fChain->SetBranchAddress("TPartCal", TPartCal, &b_CalT);
   fChain->SetBranchAddress("PTypePartCal", PTypePartCal, &b_CalPType);
   fChain->SetBranchAddress("XPartCal", XPartCal, &b_CalX);
   fChain->SetBranchAddress("YPartCal", YPartCal, &b_CalY);
   fChain->SetBranchAddress("ECluster", ECluster, &b_ECluster);
   fChain->SetBranchAddress("QCluster", QCluster, &b_QCluster);
   fChain->SetBranchAddress("XCluster", XCluster, &b_XCluster);
   fChain->SetBranchAddress("YCluster", YCluster, &b_YCluster);
   fChain->SetBranchAddress("ThCluster", ThCluster, &b_ThCluster);
   fChain->SetBranchAddress("M2Cluster", M2Cluster, &b_M2Cluster);
   fChain->SetBranchAddress("TCluster", TCluster, &b_TCluster);
   fChain->SetBranchAddress("NClusCells", NClusCells, &b_NClusCells);
   fChain->SetBranchAddress("ETarget", &ETarget, &b_ETarget);
   fChain->SetBranchAddress("TTarget", &TTarget, &b_TTarget);
   fChain->SetBranchAddress("XTarget", &XTarget, &b_XTarget);
   fChain->SetBranchAddress("YTarget", &YTarget, &b_YTarget);
   fChain->SetBranchAddress("HEPVetoTrEne", HEPVetoTrEne, &b_NTHEPVetoTrkEne);
   fChain->SetBranchAddress("HEPVetoNFing", HEPVetoNFing, &b_NTHEPVetoTrkFinger);
   fChain->SetBranchAddress("HEPVetoTrTime", HEPVetoTrTime, &b_NTHEPVetoTrkTime);
   fChain->SetBranchAddress("HEPVetoFingE", HEPVetoFingE, &b_NTHEPVetoFingE);
   fChain->SetBranchAddress("HEPVetoX", HEPVetoX, &b_NTHEPVetoX);
   fChain->SetBranchAddress("HEPVetoY", HEPVetoY, &b_NTHEPVetoY);
   fChain->SetBranchAddress("HEPVetoClIndex", HEPVetoClIndex, &b_NTHEPVetoClIndex);
   fChain->SetBranchAddress("HEPVetoECl", HEPVetoECl, &b_NTHEPVetoECl);
   fChain->SetBranchAddress("HEPVetoTimeCl", HEPVetoTimeCl, &b_NTHEPVetoTimeCl);
   fChain->SetBranchAddress("PVetoTrEne", PVetoTrEne, &b_NTPVetoTrkEne);
   fChain->SetBranchAddress("PVetoNFing", PVetoNFing, &b_NTPVetoTrkFinger);
   fChain->SetBranchAddress("PVetoTrTime", PVetoTrTime, &b_NTPVetoTrkTime);
   fChain->SetBranchAddress("PVetoFingE", PVetoFingE, &b_NTPVetoFingE);
   fChain->SetBranchAddress("PVetoX", PVetoX, &b_NTPVetoX);
   fChain->SetBranchAddress("PVetoY", PVetoY, &b_NTPVetoY);
   fChain->SetBranchAddress("PVetoBarE", PVetoBarE, &b_NTPVetoBarEnergy);
   fChain->SetBranchAddress("PVetoBarT", PVetoBarT, &b_NTPVetoBarTime);
   fChain->SetBranchAddress("PVetoClIndex", PVetoClIndex, &b_NTPVetoClIndex);
   fChain->SetBranchAddress("PVetoECl", PVetoECl, &b_NTPVetoECl);
   fChain->SetBranchAddress("PVetoTimeCl", PVetoTimeCl, &b_NTPVetoTimeCl);
   fChain->SetBranchAddress("EVetoTrEne", EVetoTrEne, &b_NTEVetoTrkEne);
   fChain->SetBranchAddress("EVetoNFing", EVetoNFing, &b_NTEVetoTrkFinger);
   fChain->SetBranchAddress("EVetoTrTime", EVetoTrTime, &b_NTEVetoTrkTime);
   fChain->SetBranchAddress("EVetoFingE", EVetoFingE, &b_NTEVetoFingE);
   fChain->SetBranchAddress("EVetoX", EVetoX, &b_NTEVetoX);
   fChain->SetBranchAddress("EVetoY", EVetoY, &b_NTEVetoY);
   fChain->SetBranchAddress("EVetoClIndex", EVetoClIndex, &b_NTEVetoClIndex);
   fChain->SetBranchAddress("EVetoECl", EVetoECl, &b_NTEVetoECl);
   fChain->SetBranchAddress("EVetoTimeCl", EVetoTimeCl, &b_NTEVetoTimeCl);
}

Bool_t UBTF_SelectorAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef UBTF_SelectorAnalysis_cxx
