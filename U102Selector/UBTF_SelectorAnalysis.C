#define UBTF_SelectorAnalysis_cxx
#include "UBTF_SelectorAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>

template<class T> float sign(T t)
{
	if (t<0.) return -1.;
	else       return 1.;
}

void UBTF_SelectorAnalysis::Begin(TTree * /*tree*/)
{
	fNProcessedEvents=0					   ;
	fNEvAtLeastOneCluster=0			   ;
	fNEv2Clusters=0						     ;
	fNEv3Clusters=0							   ;
	fNEv4Clusters=0                ;
	fNEvOneCluster=0               ;
	fNEvAtLeastOneClusterAboveThr=0;
	fNEv2ClustersAboveThr=0        ;
	fNEv3ClustersAboveThr=0        ;
	fNEv4ClustersAboveThr=0        ;
	fNEvOneClusterAboveThr=0       ;


	fNEvEphotonPass=0;
	fNEvRadiusPass=0;
	fNEvFingerPass=0;
	fNEvTimePass=0;

	sogliaCut=0;
	preSoglia=0;

	NBremmEvT=0;
	NBremmEv=0;
	NPassBremm=0;

	/*thethamax=70*pow(10,-3);//rad
	thethamin=20*pow(10,-3);//rad*/

	Rmax=(2.5+10)*2.1;
	Rmin=(2.5+2)*2.1;
	thethamax=atan(Rmax/300);
	thethamin=atan(Rmin/300);
	fAprimeMass=10.;//MeV
	fECALthresholdE = 50.;//MeV

	save=false;
	m_getCorrelations=false;
	m_useECALthreshold=true;


	// container of all histograms

	Double_t minMissM2 = -200.;
	Double_t maxMissM2 =  600.;
	Int_t    nBinMissM2 = 200;
	Double_t minMissM  = -30.;
	Double_t maxMissM  =  30.;
	Int_t    nBinMissM = 200;

	hM2noCut          = new TH1F("MissM2noCut","MissM2noCut", nBinMissM2, minMissM2, maxMissM2);
	hMnoCut           = new TH1F("MissMnoCut", "MissMnoCut" , nBinMissM , minMissM , maxMissM );
	fhVec.push_back(hM2noCut);
	fhVec.push_back(hMnoCut);
	hM2leadingECALcl  = new TH1F("MissM2leadingECALcl","MissM2leadingECALcl", nBinMissM2, minMissM2, maxMissM2);
	hMleadingECALcl   = new TH1F("MissMleadingECALcl" ,"MissMleadingECALcl" , nBinMissM , minMissM , maxMissM );
	fhVec.push_back(hM2leadingECALcl);
	fhVec.push_back(hMleadingECALcl);
	hM2subleadECALcl  = new TH1F("MissM2subleadECALcl","MissM2subleadECALcl", nBinMissM2, minMissM2, maxMissM2);
	hM2subsublECALcl  = new TH1F("MissM2subsublECALcl","MissM2subsublECALcl", nBinMissM2, minMissM2, maxMissM2);
	hM2noCutAt0       = new TH1F("MissM2noCutAt0","MissM2noCutAt0",           nBinMissM2, minMissM2, maxMissM2);
	fhVec.push_back(hM2subleadECALcl);
	fhVec.push_back(hM2subsublECALcl);
	fhVec.push_back(hM2noCutAt0);
	hM2clusterCut     = new TH1F("MissM2clusterCut", "MissM2clusterCut", nBinMissM2, minMissM2, maxMissM2);
	hMclusterCut      = new TH1F("MissMclusterCut",  "MissMclusterCut",  nBinMissM , minMissM , maxMissM );
	fhVec.push_back(hM2clusterCut);
	fhVec.push_back(hMclusterCut);


	energyCut = new TH1F("energyCut", "energyCut", 1200, -600, 600);
	radiusCut = new TH1F("radiusCut", "radiusCut", 1200, -600, 600);
	TimeCut = new TH1F("TimeCut", "TimeCut", 1200, -600, 600);
	energyPositronPhotonCut = new TH1F("energyPositronPhotonCut", "energyPositronPhotonCut", 1200, -600, 600);
	TimeError = new TH1F("TimeError", "TimeError", 10, 0, 1);
	PVetoTrEne1Error = new TH1F("PVetoTrEne1Error", "PVetoTrEne1Error", 10, 0, 1);
	controllo= new TH1F("controllo","controllo",160,0,16);
	ETotBremEv= new TH1F("ETotBremEv","ETotBremEv",100,500,800);
	ETotBremEv1= new TH1F("ETotBremEv1","ETotBremEv1",100,500,800);
	DeltaTime= new TH1F("DeltaTime","DeltaTime",60,-30,30);
	DeltaPVetoTrEne= new TH1F("DeltaPVetoTrEne","DeltaPVetoTrEne",120,-60,60);
	pVetoEClvsFinger = new TH2F("pVetoEClvsFinger", "pVetoEClvsFinger", 100,0,10, 120,0,120);
	pVetoEClvsEphoton = new TH2F("pVetoEClvsEphoton", "pVetoEClvsEphoton", 100,0,10, 400,0,220);
	Ephotonvsfinger = new TH2F("Ephotonvsfinger", "Ephotonvsfinger", 400,0,220, 150,0,150);
	TPhotonVsTVeto = new TH2F("TPhotonVsTVeto", "TPhotonVsTVeto", 100,0,60, 100,0,60);
	PVetoFingVsTveto = new TH2F("PVetoFingVsTveto","PVetoFingVsTveto",100,0.,60.,110,0.,109.);
	NfingervsPVetoTrEne = new TH2F("NfingervsPVetoTrEne", "NfingervsPVetoTrEne" ,100,0,100, 600,0,600);
	NfingervsPVetoTrEneBremm = new TH2F("NfingervsPVetoTrEneBremm", "NfingervsPVetoTrEneBremm" ,100,0,100, 600,0,600);
	xClustervseCluster = new TH2F("xClustervseCluster","xClustervseCluster",600,-300,300,500,0,500);
	yClustervseCluster=new TH2F("yClustervseCluster","yClustervseCluster",600,-300,300,500,0,500);
	xClustervsm2Cluster=new TH2F("xClustervsm2Cluster","xClustervsm2Cluster",600,-300,300,1200,-600,600);
	yClustervsm2Cluster=new TH2F("yClustervsm2Cluster","yClustervsm2Cluster",600,-300,300,1200,-600,600);
	eClustervsm2Cluster=new TH2F("eclustervsm2Cluster","eclustervsm2Cluster",500,0,500,1200,-600,600);
	xClustervsyCluster=new TH2F("xClustervsyCluster","xClustervsyCluster",600,-300,300,600,-300,300);

	if (save==true)
	{
		out1.open("Data.txt");
		out1<<"topology: 4 5 1"<< std::endl;
	}

	TString option = GetOption();

}


double UBTF_SelectorAnalysis::EvsAngle(double theta, double mA )
{
	double     me=0.5; //MeV
	double     EBeam=550;//MeV
	double     Ephoton = (2*me*EBeam-mA*mA)/(2*me+theta*theta*EBeam);

	return Ephoton;

}

void UBTF_SelectorAnalysis::SaveTrainingData(double m2Aprime, double ECluster, double XCluster, double YCluster, int go)
{
	if (save==false) return;
	else
	{
		if (abs(ECluster)>1e-6 && abs(XCluster)>0.01 && abs(YCluster)>0.01)
		{
			out1<<"in: "<< ECluster<< "  " << XCluster << "  " <<YCluster << "  "<< m2Aprime << endl;
			if (go==0) out1<<"out: 0"<< endl;
			else out1<<"out: 1"<< endl;
		}
	}
}


void UBTF_SelectorAnalysis::SlaveBegin(TTree * /*tree*/)
{
	// The SlaveBegin() function is called after the Begin() function.
	// When running with PROOF SlaveBegin() is called on each slave server.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();

}

Bool_t UBTF_SelectorAnalysis::Process(Long64_t entry)
{
	Int_t iEntry = GetEntry(entry);

	bool debug0=false;
	bool debug1=false;

	fNProcessedEvents++;
	if ( (fNProcessedEvents%1000)==0 ) debug0=true;
	if (debug0) cout<<"-------------------------------- ev n. "<<fNProcessedEvents<<" ientry = "<<iEntry<<endl;
	if (NClusters<1)
	{
		SaveTrainingData(M2Cluster[0], ECluster[0], XCluster[0],  YCluster[0], 0);
		if (debug1 || debug0) cout<<"---- this event has < 1 ECAL cluster; skipping"<<endl;
		return kTRUE;
	}
	++fNEvAtLeastOneCluster;

	int nclAboveEthr     =  0;
	int iClLeading       = -1;
	int iClsubLeading    = -1;
	int iClsubsubLeading = -1;

	double EClLeading       = 0.;
	double EClsubLeading    = 0.;
	double EClsubsubLeading = 0.;


	for (unsigned int jc=0; jc<(unsigned int)NClusters; ++jc)
	{
		hM2noCut->Fill(M2Cluster[jc]);
		hMnoCut ->Fill( sign(M2Cluster[jc]) * sqrt( fabs( M2Cluster[jc] ) ));
		if (jc==0) hM2noCutAt0 ->Fill(  M2Cluster[jc] );

		if (ECluster[jc] > fECALthresholdE)
		{
			++nclAboveEthr;
			if (nclAboveEthr==1)
			{
				EClLeading = ECluster[jc];
				iClLeading = jc;
			}
			else if (nclAboveEthr==2)
			{
				if (ECluster[jc]>EClLeading)
				{
					EClsubLeading = EClLeading;
					iClsubLeading = iClLeading;
					EClLeading = ECluster[jc];
					iClLeading = jc;
				}
				else
				{
					EClsubLeading = ECluster[jc];
					iClsubLeading = jc;
				}
			}
			else if (nclAboveEthr==3)
			{
				if (ECluster[jc]>EClLeading)
				{
					EClsubsubLeading = EClsubLeading;
					iClsubsubLeading = iClsubLeading;

					EClsubLeading = EClLeading;
					iClsubLeading = iClLeading;

					EClLeading = ECluster[jc];
					iClLeading = jc;
				}
				else if (ECluster[jc]>EClsubLeading)
				{
					EClsubsubLeading = EClsubLeading;
					iClsubsubLeading = iClsubLeading;
					EClsubLeading = ECluster[jc];
					iClsubLeading = jc;
				}
				else if (ECluster[jc]>EClsubsubLeading)
				{
					EClsubsubLeading = ECluster[jc];
					iClsubsubLeading = jc;
				}
			}
			else
			{
				cout<<"More than 3 clusters above threshold !!!!!!!!!!!"<<endl;
			}
		}
	}
	if (nclAboveEthr>0)
	{
		fNEvAtLeastOneClusterAboveThr++;
		if (nclAboveEthr>1)
		{
			fNEv2ClustersAboveThr++;
			if (nclAboveEthr>2)
			{
				fNEv3ClustersAboveThr++;
				if (nclAboveEthr>3)
				{
					fNEv4ClustersAboveThr++;
				}
			}
		}
		else
		{
			fNEvOneClusterAboveThr++;
		}
  }
	if (nclAboveEthr>0)
	{
		hM2leadingECALcl->Fill(M2Cluster[iClLeading])		;
		hMleadingECALcl ->Fill(sign(M2Cluster[iClLeading]) * sqrt( fabs( M2Cluster[iClLeading] ) ))		;

		if (nclAboveEthr>1)
		{
			hM2subleadECALcl->Fill(M2Cluster[iClsubLeading])		;
			if (nclAboveEthr>2)
			{
				hM2subsublECALcl->Fill(M2Cluster[iClsubsubLeading])		;
			}
		}

	}

	if (NClusters>1)
	//if (nclAboveEthr>1)
	{
		SaveTrainingData(M2Cluster[0], ECluster[0], XCluster[0],  YCluster[0], 0);
		if (debug1 || debug0) cout<<"---- this event has > 1 ECAL cluster; skipping"<<endl;
		switch(NClusters)
		{
			case(2):
			fNEv2Clusters++;
			break;
			case(3):
			fNEv3Clusters++;
			break;
			case(4):
			fNEv4Clusters++;
			break;
		}
		//return kTRUE;
	}
	else
	{
		++fNEvOneCluster;
	}

	Int_t iCL=-1;
	if (m_useECALthreshold)
	{
		if (nclAboveEthr!=1)
		{
			return kTRUE;
		}
		else
		{
			iCL = iClLeading;
		}
  }
	else
	{
  	if (NClusters>1) {
			return kTRUE;
		}
		else
		{
			iCL = 0;
		}
	}
	if (iCL<0)
	{
		std::cout<<"ERROR iCL <0 "<<std::endl;
		return kTRUE;
	}
	photonE = ECluster[iCL];
	Double_t missM2 = M2Cluster[iCL];
	Double_t missM  = sign(missM2)*sqrt(fabs(M2Cluster[iCL]));
	hM2clusterCut->Fill(missM2);
	hMclusterCut ->Fill(missM );

	/*
	*
	* riempimento istogrammi di controllo
	*
	*
	*/

	if(ECluster[iCL]>0.01)
	{
		xClustervseCluster->Fill(XCluster[iCL],ECluster[iCL]);
		yClustervseCluster->Fill(YCluster[iCL],ECluster[iCL]);
		eClustervsm2Cluster->Fill(ECluster[iCL],M2Cluster[iCL]);
	}
	xClustervsm2Cluster->Fill(XCluster[iCL],M2Cluster[iCL]);
	yClustervsm2Cluster->Fill(YCluster[iCL],M2Cluster[iCL]);
	xClustervsyCluster->Fill(XCluster[iCL],YCluster[iCL]);




	for (int i=0; i<100; i++)
	{
		if(abs(PVetoFingE[i])>0.01)
		{
			//cout<<entry<< "    #entry   con PVetoECl  " << i << " i-simo(finger) con j(iCluster):  "
			//    <<j<< " pari a  " << PVetoECl[i][j] <<endl;
			pVetoEClvsFinger->Fill(PVetoFingE[i],PVetoNFing[i]);
			pVetoEClvsEphoton->Fill(PVetoFingE[i], ECluster[iCL]);
			Ephotonvsfinger->Fill(ECluster[iCL],PVetoNFing[i]);
			//cout<<"N.Finger  "<< i <<"  -sima è  "<< PVetoNFing[i]<< endl;

		}
	}



	Emin = EvsAngle(thethamax, fAprimeMass);
	Emax = EvsAngle(thethamin, fAprimeMass);
	//cout << "energia:   minima    "<< Emin<< "    massima:      "<<Emax<< endl;
	if (photonE<Emin || photonE>Emax)
	{
		SaveTrainingData(M2Cluster[iCL], ECluster[iCL], XCluster[iCL],  YCluster[iCL], 0);
		if (debug1 || debug0) cout<<"---- this event has ECAL cluster energy = "<<photonE<<" out of range ("<<Emin<<","<<Emax<<") for M_Uboson="<<fAprimeMass<<"; skipping"<<endl;
		return kTRUE;
	}
	energyCut->Fill(M2Cluster[iCL]);
	++fNEvEphotonPass;


	radius=sqrt(XCluster[iCL]*XCluster[iCL]+YCluster[iCL]*YCluster[iCL]);
	if (radius < Rmin*10 || radius>Rmax*10)
	{
		SaveTrainingData(M2Cluster[iCL], ECluster[iCL], XCluster[iCL],  YCluster[iCL], 0);
		if (debug1 || debug0) cout<<"---- this event has ECAL cluster out of FV; skipping"<<endl;
		return kTRUE;
	}
	radiusCut->Fill(M2Cluster[iCL]);
	++fNEvRadiusPass;

	for (int i=0; i<NPVetoTracks; i++)
	{
		if (PVetoTrEne[i]>10.)   NfingervsPVetoTrEne->Fill(PVetoNFing[i],PVetoTrEne[i]); // for e+ E>10 MeV
		TPhotonVsTVeto->Fill(PVetoTrTime[i],TCluster[iCL]); //
		PVetoFingVsTveto->Fill(PVetoTrTime[i],PVetoNFing[i]); //
	}


	/*
	intercetta TPhotonVsTVetoCut 26.2071+-0.460039
	coeff angolare TPhotonVsTVetoCut 0.186076+-0.0175328*/
	double LinFitP0=5.96;
	double LinFitP1=1.01;
	double LinFitP0Err=0.460;
	double LinFitP1Err=0.0175;

	//parte relativa al fit di NfingervsPVetoTrEneProfile
	double a= 14.699;
	double aErr= 0.133;
	double b= 0.849;
	double bErr= 0.012;
	double c= 0.039663;
	double cErr= 0.00017;

	int nbremmBACO=0;
	int nbremm=0;

	int nbremmT=0;
	for (int iFing=0; iFing<100; iFing++)
	{
		for (int iCluster =0; iCluster<10; iCluster++)
		{
			preSoglia++;
			controllo->Fill(PVetoECl[iFing][iCluster]);
			if(PVetoECl[iFing][iCluster]<1.) continue;  //soglia di 1MeV
			sogliaCut++;
			double TVeto=PVetoTimeCl[iFing][iCluster];     //Controllo sul tempo Delta t < 2ns è bremmst.
			double TExpBremmFoton=LinFitP0+TVeto*LinFitP1;
			double TExpBremmFotonError=pow(LinFitP0Err*LinFitP0Err+TVeto*TVeto*LinFitP1Err*LinFitP1Err,0.5);
			DeltaTime->Fill(TExpBremmFoton-TCluster[iCL]);
			TimeError->Fill(TExpBremmFotonError);
			if(abs(TExpBremmFoton-TCluster[iCL])<=2.)
			{
				nbremmT++;                           // numero di hit nel pVeto in tempo con il cluster in ECAL
				//ora devo fare il controllo sull'energia
				double ERaccolta=PVetoECl[iFing][iCluster];
				double FingEneFoton= a + b*PVetoNFing[iFing]+c*PVetoNFing[iFing]*PVetoNFing[iFing];
				double Fing2PKinE  = a + b*double(iFing)+c*double(iFing)*double(iFing);

				double FingEneFotonErr= pow(aErr*aErr+bErr*bErr*PVetoNFing[iFing]+cErr*cErr*PVetoNFing[iFing]*PVetoNFing[iFing],0.5);
				double DeltaTrEne  = 550. - FingEneFoton - ECluster[iCL];
				double DeltaTrEneSS= 550. - Fing2PKinE - ECluster[iCL];

				if (debug1 || debug0)std::cout<<" in time positron hit and cluster - nfing, pEneFromFing, E(ECAL), pEne+E(ECAL) "<<PVetoNFing[iFing]<<" "<<FingEneFoton<<" "<<ECluster[iCL]<<" "<<DeltaTrEne<<std::endl;
				if (debug1 || debug0)std::cout<<" in time positron hit and clusterSS:nfing, pEneFromFing, E(ECAL), pEne+E(ECAL) "<<iFing<<" "<<Fing2PKinE<<" "<<ECluster[iCL]<<" "<<DeltaTrEneSS<<std::endl;
				DeltaPVetoTrEne->Fill(DeltaTrEne);
				PVetoTrEne1Error->Fill(FingEneFotonErr);
				if (debug1 || debug0)cout << "somma delle energie " << ECluster[iCL]+ FingEneFoton << endl;
				if(ECluster[iCL]+FingEneFoton > 500 && ECluster[iCL]+FingEneFoton <650)
				{
					if (debug1 || debug0)cout <<"Efoton+FingEneProton " <<ECluster[iCL]+FingEneFoton<<endl;
					nbremmBACO++;
					// break;
				}
				if(ECluster[iCL]+Fing2PKinE > 500 && ECluster[iCL]+Fing2PKinE <650)
				{
					if (debug1 || debug0)cout <<"Efoton+Fing2PKinE " <<ECluster[iCL]+Fing2PKinE<<endl;
					nbremm++;
					// break;
				}
			}



		}

	}

	bool PassBremm=true;
	if (debug1 || debug0)cout << "nbremmT " << nbremmT << " nbremmBACO "  << nbremmBACO << " nbremm "<<nbremm<<endl;
	if(nbremmT>0)
	{
		NBremmEvT++;
		if (nbremm>0)
		{
			NBremmEv++;
			PassBremm=false;
		}
	}
	if(!PassBremm) return kTRUE;
	NPassBremm ++;
	TimeCut->Fill(M2Cluster[iCL]);


	//  energyPositronPhotonCut->Fill(M2Cluster[iCL]);
if (debug1 || debug0)	cout <<"--------------------------------------------------------------"<<endl;



	return kTRUE;

}

void UBTF_SelectorAnalysis::SlaveTerminate()
{
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.

}

void UBTF_SelectorAnalysis::Terminate()
{
	cout<<" "<<endl;
	cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%% In UBTF_SelectorAnalysis::Terminate %%%%%%%%%%"<<endl;
	cout<<" "<<endl;


	TFile *outFile = new TFile("Histogram.root", "RECREATE");
	if(!outFile->IsOpen())
	{
		cout<<"Histogram.root non aperto correttamente"<<endl;
		return;
	}
	else
	{
		for (std::vector<TH1F*>::const_iterator it=fhVec.begin(); it!=fhVec.end(); ++it)
		{
			std::cout<<"Writing to output root file histogram named "<<(*it)->GetName()<<std::endl;
			(*it)->Write();
		}
	}
	/***controllo***/
	/*
	xClustervseCluster->Write();
	yClustervseCluster->Write();
	xClustervsm2Cluster->Write();
	yClustervsm2Cluster->Write();
	eClustervsm2Cluster->Write();
	xClustervsyCluster->Write();

	hMnoCut->Write();
	clusterCut->GetXaxis()->SetTitle("m2A'(MeV2)");
	clusterCut->Write();
	energyCut->GetXaxis()->SetTitle("m2A'(MeV2)");
	energyCut->Write();
	radiusCut->GetXaxis()->SetTitle("m2A'(MeV2)");
	radiusCut->Write();
	TimeCut->GetXaxis()->SetTitle("m2A'(MeV2)");
	TimeCut->Write();
	energyPositronPhotonCut->GetXaxis()->SetTitle("m2A'(MeV2)");
	energyPositronPhotonCut->Write();
	//controllo->GetXaxis()->SetTitle("m2A'(MeV)");
	//controllo->Write();

	pVetoEClvsFinger->GetXaxis()->SetTitle("PVetoECluster");
	pVetoEClvsFinger->GetYaxis()->SetTitle("finger");
	pVetoEClvsFinger->Write();
	pVetoEClvsEphoton->GetXaxis()->SetTitle("PVetoECluster");
	pVetoEClvsEphoton->GetYaxis()->SetTitle("Ephoton");
	pVetoEClvsEphoton->Write();
	Ephotonvsfinger->GetXaxis()->SetTitle("Ephoton");
	Ephotonvsfinger->GetYaxis()->SetTitle("finger");
	Ephotonvsfinger->Write();
	TPhotonVsTVeto->GetXaxis()->SetTitle("timePVetoECl");
	TPhotonVsTVeto->GetYaxis()->SetTitle("timePhoton");
	TPhotonVsTVeto->Write();
*/
	if (m_getCorrelations) {
 	TPhotonVsTVetoCut=TPhotonVsTVeto->ProfileX();
	TPhotonVsTVetoCut->Fit("pol1");
	TCanvas *cTIME = new TCanvas("cTimeCorrelation","cTimeCorrelation");
	TPhotonVsTVetoCut->Draw();
	cTIME->SaveAs("cTimeCorrelation.png");
	TF1 *fit = (TF1*)TPhotonVsTVetoCut->GetFunction("pol1");
	Intercetta=fit->GetParameter(0);
	IntercettaError=fit->GetParError(0);
	cout<< "intercetta "<< Intercetta << "+-" << IntercettaError  << endl;
	CoeffAngolare=fit->GetParameter(1);
	CoeffAngolareError=fit->GetParError(1);
	cout<< "coeff angolare "<< CoeffAngolare << "+-" << CoeffAngolareError << endl;
	TPhotonVsTVetoCut->GetXaxis()->SetTitle("timePVetoECl");
	TPhotonVsTVetoCut->GetYaxis()->SetTitle("timePhoton");
	TPhotonVsTVetoCut->Write();


	PVetoTrEnevsNfingerProfile= NfingervsPVetoTrEne->ProfileX();
	if (PVetoTrEnevsNfingerProfile->GetEntries()==0) return ;
	PVetoTrEnevsNfingerProfile->Fit("pol2");
	TF1 *fit0 = (TF1*)PVetoTrEnevsNfingerProfile->GetFunction("pol2");
	cout<<"a= "<<fit0->GetParameter(0)<<endl;
	cout<<"aErr= "<<fit0->GetParError(0)<<endl;
	cout<<"b= "<<fit0->GetParameter(1)<<endl;
	cout<<"bErr= "<<fit0->GetParError(1)<<endl;
	cout<<"c= "<<fit0->GetParameter(2)<<endl;
	cout<<"cErr= "<<fit0->GetParError(2)<<endl;
	TCanvas *cEFING = new TCanvas("cEnergyFingerCorrelation","cEnergyFingerCorrelation");
	PVetoTrEnevsNfingerProfile->Draw();
	cEFING->SaveAs("cEnergyFingerCorrelation.png");

	PVetoTrEnevsNfingerProfile->Write();
	}

/*
	NfingervsPVetoTrEne->Write();

	DeltaPVetoTrEne->Write();
	PVetoFingVsTveto->Write();
	PVetoTrEne1Error->Write();
	DeltaTime->Write();
	TimeError->Write();
	controllo->Write();
*/

	outFile->ls();
	outFile->Close();
	delete outFile;
	if (save==true) out1.close();


	PrintCounters();
	cout<<" "<<endl;
	cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of UBTF_SelectorAnalysis::Terminate %%%%%%%%%%"<<endl;
	cout<<" "<<endl;


}
void UBTF_SelectorAnalysis::PrintCounters()
{
	cout<<" "<<endl;
	cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%% In UBTF_SelectorAnalysis::PrintCounters %%%%%%%%%%"<<endl;
	cout<<" "<<endl;

	cout<<" fNProcessedEvents =========== "<<fNProcessedEvents<<std::endl;
	cout<<" fNEvAtLeastOneCluster ======== "<<fNEvAtLeastOneCluster<<endl;
	cout<<" fNEvOneCluster        -------- "<<fNEvOneCluster<<endl;
	cout<<" fNEv2Clusters         -------- "<<fNEv2Clusters<<endl;
	cout<<" fNEv3Clusters         -------- "<<fNEv3Clusters<<endl;
	cout<<" fNEv >3 Clusters      -------- "<<fNEv4Clusters<<endl;
	cout<<" fNEvAtLeastOneClusterAboveThr= "<<fNEvAtLeastOneClusterAboveThr<<" ---- Threshold = "<<fECALthresholdE<<" MeV"<<endl;
	cout<<" fNEvOneClusterAboveThr   ----- "<<fNEvOneClusterAboveThr<<endl;
	cout<<" fNEv2ClustersAboveThr    ----- "<<fNEv2ClustersAboveThr<<endl;
	cout<<" fNEv3ClustersAboveThr    ----- "<<fNEv3ClustersAboveThr<<endl;
	cout<<" fNEv >3 ClustersAboveThr ----- "<<fNEv4ClustersAboveThr<<endl;


	cout<<" fNEvEphotonPass    =========== "<<fNEvEphotonPass<<std::endl;
	cout<<" fNEvRadiusPass     =========== "<<fNEvRadiusPass<<std::endl;
	cout<<" fNEvFingerPass     =========== "<<fNEvFingerPass<<std::endl;
	cout<<" fNEvTimePass       =========== "<<fNEvTimePass<< endl;
	cout<<" NBremmEvT  ----- " << NBremmEvT
	<< "/n NBremmEv   ----- " << NBremmEv
	<< "/n NPassBremm ----- " << NPassBremm << endl;
	//   cout<<"prima "<<preSoglia << "soglia Cut  "<< sogliaCut << " time Cut " << NEvTimePass
	//        << "energy positron e photon " << fNEvFingerPass << endl;
	cout<<" "<<endl;
	cout<<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of UBTF_SelectorAnalysis::PrintCounters %%%%%%%%%%"<<endl;
	cout<<" "<<endl;
	return;
}
