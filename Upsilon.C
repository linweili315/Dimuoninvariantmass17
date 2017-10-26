#include <TCanvas.h>
#include <TH1D.h>
#include <TH3.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TCutG.h>
#include <Math/Functor.h>
#include <TPad.h>
#include <TString.h>

#include <TMinuit.h>
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooAbsPdf.h>
#include <RooHistPdf.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooPolynomial.h>
#include <RooMCStudy.h>
#include <RooMinuit.h>
#include <RooWorkspace.h>
#include <RooConstVar.h>
#include <RooRandom.h>
#include <RooDataHist.h>
#include <RooFunctorBinding.h>
#include <RooStats/RooStatsUtils.h>

using namespace RooFit ;


void Upsilon()
{
    TFile *f = new TFile(“***.root");
    TH1F *h1 = (TH1F*)f->Get("Upsilon_Soft_Inner_AllHtlPaths/hInvMass_fEtaMin0.  _fEtaMax2.4 _Upsilon_Soft_Inner_AllHtlPaths");
 
    RooRealVar mass("mass","mass",8.7,11.0);
    RooDataHist data("data","data",mass,Import(*h1));

   //////////////upsilon (1S)//////////////////////////////////////////////
    RooRealVar cbmean("mean", "mean",9.49,9.41,9.52);
    RooRealVar cbwidth("width", "width",0.0,0.08);
    RooRealVar cbn("n","n",0,100);
    RooRealVar cbalpha("alpha","", 0,10);
    RooCBShape cball("cball","cball",mass, cbmean, cbwidth, cbalpha, cbn);

    RooRealVar cbwidth0("width0", "width0",0.0,0.15);
    RooCBShape cball0("cball0","cball0",mass, cbmean, cbwidth0, cbalpha, cbn); 

    RooRealVar Snx("Snx", "Snx",1000,0,10000);
    RooRealVar Snx1("Snx1", "Snx1",700,0,10000);
    RooAddPdf  signal("signal", "signal", RooArgList(cball,cball0), RooArgList(Snx,Snx1));

   //////////////upsilon (2S)//////////////////////////////////////////////
    RooRealVar cbmean1("mean1", "mean1",10.01,10.07);
    RooRealVar cbwidth1("width1", "width1",0.0,0.1);
    RooRealVar cbn1("n1","n1",0.0,10);
    RooRealVar cbalpha1("alpha1","", 0,10);   
    RooCBShape cball1("cball1","cball1",mass, cbmean1, cbwidth1, cbalpha1, cbn1);


   //////////////upsilon (3S)//////////////////////////////////////////////
    RooRealVar cbmean2("mean2", "mean2",10.34,10.32,10.35);
    RooRealVar cbwidth2("width2", "width2",0.0,0.1);
    RooRealVar cbn2("n2","n2",0.,10);
    RooRealVar cbalpha2("alpha2","",0,10);   
    RooCBShape cball2("cball2","cball2",mass, cbmean2, cbwidth2, cbalpha2, cbn2);

   ///////////////background ////////////////////////////////////////////////////////
    RooRealVar poly0("poly0","poly0",-1.0e100,1.0e100);
    RooRealVar poly1("poly1","poly1",-1.0e100,1.0e100);
    RooRealVar poly2("poly2","poly2",-1.0e100,1.0e100);
    RooRealVar poly3("poly3","poly3",-1.0e100,1.0e100);
    RooChebychev  chebychev("chebychev", "chebychev background ",mass, RooArgList(poly0, poly1));

   //////////////total pdf /////////////////////////////////////////////////////////////
    RooRealVar nx("nx", "nx",120000,0.,2000000);
    RooRealVar nx1("nx1", "nx1",10000,0.,300000);
    RooRealVar nx2("nx2", "nx2",5000,0.,200000);
    RooRealVar nx3("nx3", "nx3",100000,10000.,1500000);

    RooAddPdf  total("total", "total", RooArgList(signal,cball1,cball2,chebychev),RooArgList(nx,nx1,nx2,nx3));
    RooFitResult* FitResult1 = total.fitTo(data,Strategy(0),Save(true),Extended(1));
    FitResult1->Print("v");

   /////////// count sigma ////////////////////////////////////////   
    double f1 = Snx.getVal()/(Snx.getVal() + Snx1.getVal());
    double f2 = 1 - f1;
    double sigma = TMath::Sqrt( f1*cbwidth.getVal()*cbwidth.getVal() + f2*cbwidth0.getVal()*cbwidth0.getVal());
    cout << sigma <<endl;
    setTDRStyle();
    TCanvas *c2=new TCanvas("c2","c2",800,600);
    c2->cd();

    TPad* pad1 = new TPad("pad1","pad1",0.,0.,1.,1.);
    setTDRStyle();
    pad1->Draw();
    pad1->cd();    
    RooPlot* xf = mass.frame(Range(8.7,11.0));
  ///////////////////set axis////////////////////////////////////
    xf->GetXaxis()->SetTitle("#mu^{+}#mu^{-} invariant mass [GeV]");
    xf->GetYaxis()->SetTitle("Events / 20 MeV");
     TGaxis::SetMaxDigits(3);
    data.plotOn(xf);
    total.plotOn(xf ,LineColor(1));
    total.plotOn(xf ,Components("chebychev"),LineStyle(kDashed),LineColor(kBlue-7));
    total.plotOn(xf ,LineColor(1));
    total.plotOn(xf ,Components("total"),LineStyle(1),LineColor(kRed));

    TPaveText* paveTextP = new TPaveText(0.65,0.6,0.9,0.85,"NDC");
    paveTextP->AddText(Form("%s%.f%s","#sigma = ",1000*sigma," MeV"));
    paveTextP->AddText(Form("%s%.f%s","p_{T}(#mu^{+}#mu^{-}) > ",12.," GeV"));
    paveTextP->AddText(Form("%s%.1f","|#eta(#mu^{+/-})| < ", 1.5));
    paveTextP->SetTextAlign(10);
    paveTextP->SetBorderSize(0.0);
    paveTextP->SetFillStyle(0);
    paveTextP->SetTextSize(0.05);
    paveTextP->Paint();
    xf->addObject(paveTextP);
    xf->Draw();
    CMS_lumi(pad1, 4, 10);
    c2->Update();

}



