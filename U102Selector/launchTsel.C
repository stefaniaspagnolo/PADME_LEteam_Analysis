#include "TChain.h"
//#include "UBTF_SelectorAnalysis.h"
/// introducing a cahneg into the master branch !!!!!!

void launchTsel(){

//gROOT->ProcessLine(".L UBTF_selector.C++");
TChain *f = new TChain("U102");
//f->Add("/einstein2/stefania/padme/prod/test_20170707_17_19/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/test_20170707_17_20/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/test_20170707_17_21/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/test_20170707_17_22/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/test_20170707_17_43/UBTF.root");
//f->Add("/einstein2/stefania/padme/prod/test_20170707_17_44/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170707_17_22/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170707_17_19/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170707_17_20/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170707_17_21/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170707_17_43/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170707_17_44/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170708_07_19_33384evs/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_13_43/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_50/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_29_36/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_36/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_53/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_32/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_03/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_29_32/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_29_21/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_29/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_25/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_30_55/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_44/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_39/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_16/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_20/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_35_09/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_29/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_11/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_47/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_24/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_30_59/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_31_07/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_36/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_45/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_09/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_15/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_21/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_01/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_35_03/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_49/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_40/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_42/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_33/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_35_06/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_55/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_12/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_53/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_59/UBTF.root");
f->Add("/einstein2/stefania/padme/prod/test_20170711_16_34_18/UBTF.root");

UBTF_SelectorAnalysis a;

f->Process(&a);

}
