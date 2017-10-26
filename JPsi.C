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


void JPsi()
{
    TFile *f = new TFile(“****.root");
    TH1F *h1 = (TH1F*)f->Get("JPsi_Soft_Inner_AllHtlPaths/hInvMass_fEtaMin0.  _fEtaMax2.4 _JPsi_Soft_Inner_AllHtlPaths");
 
    RooRealVar mass("mass","mass",2.9,3.3);
    RooDataHist data("data","data",mass,Import(*h1));
  
   //////////////crystal ball fuction//////////////////////////////////////////////
    RooRealVar cbmean("mean", "mean",3.08,3.12);
    RooRealVar cbwidth("cbwidth", "cbwidth",0.03,0.0,0.05);
    RooRealVar cbn("n","n",0,10);
    RooRealVar cbalpha("alpha","",0.,10.);
    RooCBShape cball("cball","cball",mass, cbmean, cbwidth, cbalpha, cbn);

    RooRealVar cbwidth1("cbwidth1", "cbwidth1",0.04,0.0,0.08);
    RooCBShape cball0("cball0","cball0",mass, cbmean, cbwidth1, cbalpha, cbn);

    RooRealVar Snx("Snx", "Snx",1000,0,10000);
    RooRealVar Snx1("Snx1", "Snx1",1000,0,10000);

    RooAddPdf  signal("signal", "signal", RooArgList(cball,cball0), RooArgList(Snx,Snx1));

   ///////////////background ////////////////////////////////////////////////////////
    RooRealVar poly0("poly0","poly0",-1.0e100,1.0e100);
    RooRealVar poly1("poly1","poly1",-1.0e100,1.0e100);
    RooRealVar poly2("poly2","poly2",-1.0e100,1.0e100);
    RooRealVar poly3("poly3","poly3",-1.0e100,1.0e100);
    RooChebychev  chebychev("chebychev", "chebychev background ",mass, RooArgList(poly0, poly1));
    RooPolynomial poly("poly", "Polynomial background ",mass, RooArgList(poly0, poly1));

   //////////////total pdf /////////////////////////////////////////////////////////////
    RooRealVar nx("nx", "nx",2500000, 0, 3000000);
    RooRealVar nx1("nx1", "nx1",500000, 0,1000000);
    
    RooAddPdf  total("total", "total", RooArgList(signal,chebychev), RooArgList(nx,nx1));
    
    RooFitResult* FitResult1 = total.fitTo(data,Strategy(1),Save(true),Range(2.9,3.3),Extended(1));
    FitResult1->Print("v");

  /////////// count sigma ////////////////////////////////////////   
    double f1 = Snx.getVal()/(Snx.getVal() + Snx1.getVal());
    double f2 = 1 - f1;
    double sigma = TMath::Sqrt( f1*cbwidth.getVal()*cbwidth.getVal() + f2*cbwidth1.getVal()*cbwidth1.getVal());

  ////////////////////////////////////////////////////////////
    
    TCanvas *c2=new TCanvas("c2","c2",800,600);
    c2->cd();
    setTDRStyle();
  //////////////////////TPad ////////////////////////////////////
    TPad* pad1 = new TPad("pad1","pad1",0,0,1,1);
    pad1->Draw();
    pad1->cd(); 
    RooPlot* xf = mass.frame(Range(2.9,3.3));
    fixOverlay();
  ///////////////////set axis////////////////////////////////////
    xf->GetXaxis()->SetTitle("#mu^{+}#mu^{-} invariant mass [GeV]");
    xf->GetYaxis()->SetTitle("Events / 10 MeV");

    data.plotOn(xf);
    total.plotOn(xf ,LineColor(1));
    total.plotOn(xf ,Components("chebychev"),LineStyle(kDashed),LineColor(kBlue-7));
    total.plotOn(xf ,LineColor(1));
    total.plotOn(xf ,Components("total"),LineStyle(1),LineColor(kRed));

    TPaveText* paveTextP = new TPaveText(0.65,0.6,0.9,0.85,"NDC");
    paveTextP->AddText(Form("%s%.f%s","#sigma = ",1000*sigma," MeV"));
    paveTextP->AddText(Form("%s%.f%s","p_{T}(#mu^{+}#mu^{-}) > ",20.," GeV"));
    paveTextP->AddText(Form("%s%.2f","#left|y(#mu^{+}#mu^{-})#right| < ",1.25));
    paveTextP->SetTextAlign(11);
    paveTextP->SetBorderSize(0.0);
    paveTextP->SetFillStyle(0);
    paveTextP->SetTextSize(0.05);
    paveTextP->Paint();
    xf->addObject(paveTextP);   
    xf->Draw();

    CMS_lumi(pad1, 4, 10);

}


