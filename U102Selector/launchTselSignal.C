#include "TChain.h"
//#include "UBTF_SelectorAnalysis.h"

void launchTselSignal(){

//gROOT->ProcessLine(".L UBTF_selector.C++");
TChain *f = new TChain("U102");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170720_12_12_52/UBTF.root");

//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_01_33/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_05_02/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_04_39/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_06_41/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/prodUboson_20170724_17_05_35/UBTF.root");

f->Add("/einstein2/stefania/padme/prod/prodUboson_20170725_16_31_30/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/prodUboson_20170725_16_33_32/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/prodUboson_20170725_16_33_42/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/prodUboson_20170725_16_33_37/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/prodUboson_20170725_16_33_46/UBTF.root");

UBTF_SelectorAnalysis a;

f->Process(&a);

}
