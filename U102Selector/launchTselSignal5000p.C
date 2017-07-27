#include "TChain.h"
//#include "UBTF_SelectorAnalysis.h"
// now new cahnegs in the local develo prach 

void launchTselSignal5000p(){

//gROOT->ProcessLine(".L UBTF_selector.C++");
TChain *f = new TChain("U102");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170720_12_12_52/UBTF.root");

//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_01_33/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_05_02/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_04_39/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_06_41/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_05_35/UBTF.root");

f->Add("/einstein2/stefania/padme/prod/prodUboson5000p_20170725_17_10_38/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/prodUboson5000p_20170725_17_10_52/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/prodUboson5000p_20170725_17_10_44/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/prodUboson5000p_20170725_17_10_56/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/prodUboson5000p_20170725_17_10_49/UBTF.root");

UBTF_SelectorAnalysis a;

f->Process(&a);

}
