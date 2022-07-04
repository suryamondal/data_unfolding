/*
  
  Compile with:
  g++ -g -pthread -m64 -Wno-deprecated -std=c++11 -I/home/surya/products/root-6.14.00/include -o anal_muon_unfold_basic_half anal_muon_unfold_basic_half.cc `root-config --cflags` `root-config --libs` -I/home/surya/products/RooUnfold/src/ -L/home/surya/products/RooUnfold/ -lRooUnfold

  Execute with:
  ./anal_muon_unfold_basic_half coorelation_half.root
  
  ******************************
  
  RooUnfoldResponse(const TH1* measured, const TH1* truth, const TH2* response, const char* name = 0, const char* title = 0)
  
  RooUnfoldResponse constructor - create from already-filled histograms.
  "response" gives the response matrix, measured X truth.
  "measured" and "truth" give the projections of "response"
  onto the X-axis and Y-axis respectively,
  but with additional entries in "measured" for measurements with
  no corresponding truth (fakes/background) and
  in "truth" for unmeasured events (inefficiency).
  "measured" and/or "truth" can be specified as 0 (1D case only) or
  an empty histograms (no entries) as a shortcut
  to indicate, respectively, no fakes and/or no inefficiency.
  
  ********************************
  
  RooUnfoldBayes(const RooUnfoldResponse* res, const TH1* meas, Int_t niter = 4, Bool_t smoothit = false, const char* name = 0, const char* title = 0)
  
  Constructor with response matrix object and measured unfolding input histogram.
  The regularisation parameter is niter (number of iterations).
  
  ********************************


*/

#include "TDecompSVD.h"
// #include "CLHEP/Vector/TwoVector.h"
// #include "CLHEP/Vector/ThreeVector.h"
// #include "CLHEP/Vector/LorentzVector.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <TStyle.h>
#include "TPostScript.h"
#include "TLegend.h"
#include "TRandom.h"

#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldResponse.h"
//#include "RooUnfHistoSvd.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldInvert.h"

/*
  +ve or -ve muon 
*/
// #define isP


bool pedsub=true;
using namespace std;
// using namespace CLHEP;

const char* namex;
char namey[300];
const char* titlex;
char titley[300];

const int nbinmx=10000;
const int nitermx=10;
int icol[12] = {2,3,3,4,4,5,5,6,6,7,7,1};


// int subtract_background(TH2D h2d_correl, TH1D mcd, double* fakerate, double* effi) {

//   int nbinx = h2d_correl->GetNbinsX();
//   int nbiny = h2d_correl->GetNbinsY();
//   if (nbinx>nbinmx || nbiny > nbinmx) {
//     cout <<"Increase nbinmx, which is "<< nbinmx<<" in subtract_background"<<endl;
//     cout <<"Reconstructed objects are not pedestal subtracted"<<endl;
//     return 0;
//   }
//   double totalgen[nbinmx]={0.};
//   double totalreco[nbinmx]={0.};
//   for (int ix=0; ix<nbinx+1; ix++) {
//     for (int iy=0; iy<nbiny+1; iy++) {
//       if (ix==nbinx && iy==nbiny) continue;
//       totalreco[ix] +=h2d_correl->GetBinContent(ix+1, iy+1);
//       if (iy==nbiny) fakerate[ix] =h2d_correl->GetBinContent(ix+1, iy+1);
      
//       totalgen[iy] +=h2d_correl->GetBinContent(ix+1, iy+1);
//       if (ix==nbinx) effi[iy] =h2d_correl->GetBinContent(ix+1, iy+1);
      
//       // cout <<"ix "<< ix<<" iy "<<iy
//       // 	   <<" binc "<<h2d_correl->GetBinContent(ix+1, iy+1)
//       // 	   <<" gen "<<totalgen[iy]<<" fake "<<fakerate[ix]
//       // 	   <<" reco "<<totalreco[ix]<<" effi "<<effi[iy]<<endl;

//       //GMA 27th Jan 2013
//       if (ix==nbinx || iy==nbiny) {
//       	h2d_correl->SetBinContent(ix+1, iy+1, 0.0);
//       	h2d_correl->SetBinError(ix+1, iy+1, 0.0);
//       }
//     } // for (int iy=0; iy<nbiny+1; iy++) {
//   } // for (int ix=0; ix<nbinx+1; ix++) {

//   for (int iy=0; iy<nbiny; iy++) {
//     effi[iy] = (totalgen[iy] - effi[iy])/max(1.0, totalgen[iy]);
//     if (effi[iy]<1.e-3) effi[iy]=1.e-3;
//   }
//   for (int ix=0; ix<nbinx; ix++){
//     fakerate[ix] /=max(1.0, totalreco[ix]);
//     // cout << " " << ix << " " << fakerate[ix] << endl;
//   }

//   if(pedsub) { 
//     for(int ix=0; ix<nbinx; ix++) {
//       mcd->SetBinContent(ix+1, (1.0 - fakerate[ix])*mcd->GetBinContent(ix+1));
//       mcd->SetBinError(ix+1, (1.0 - fakerate[ix])*mcd->GetBinError(ix+1));
//     }
//   }
//   return 1;
// }


void add_nameindex(TNamed* hist, TNamed* org, char ttl[]) {
  // void add_nameindex(TH1* hist, TH1* org, char ttl[]) {
  
  namex = org->GetName();
  sprintf(namey, "%s_%s", ttl, namex);
  hist->SetName(namey);
  
  titlex = org->GetTitle();
  sprintf(titley, "%s_%s", ttl, titlex);
  hist->SetTitle(titley);
}



int main(int argc, char**argv) {
  
  gStyle->SetPaintTextFormat("5.2e"); //4.1f    
  //   minorc_JetJPT_Gen_c1_e1->Draw("text0")

  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetStatColor(10);
  
  gStyle->SetCanvasColor(10);
  gStyle->SetOptStat(0); //1110);
  gStyle->SetOptTitle(1);

  gStyle->SetTitleW(.96);
  gStyle->SetTitleH(.06);
  gStyle->SetTitleY(0.99);
  gStyle->SetTitleX(.02);
  gStyle->SetTitleAlign(13);
  gStyle->SetTitleColor(10);
  //  gStyle->SetTitleOffset(-0.05);
  gStyle->SetTitleBorderSize(0); //1);
  gStyle->SetTitleFontSize(0.10);

  gStyle->SetPalette(1,0);
  gStyle->SetPadColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatBorderSize(1);

  gStyle->SetStatStyle(1001);
  gStyle->SetOptFit(101);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetHistLineWidth(1);

  gStyle->SetStatX(.99);
  gStyle->SetStatY(.99);
  gStyle->SetStatW(.45);
  gStyle->SetStatH(.16);
  gStyle->SetLabelSize(0.18,"XY");  
  gStyle->SetLabelOffset(0.21,"XYZ");
  gStyle->SetTitleSize(0.095,"XY");  
  gStyle->SetTitleOffset(0.06,"XYZ");
  gStyle->SetPadTopMargin(0.075); //2); //0.11 //0.09
  gStyle->SetPadBottomMargin(0.075);
  gStyle->SetPadLeftMargin(0.09);
  gStyle->SetPadRightMargin(0.10);
  //  gStyle->SetPadGridX(1);
  //  gStyle->SetPadGridY(1);
  //  gStyle->SetGridStyle(3);
  //  gStyle->SetNdivisions(606,"XY");

  gStyle->SetMarkerSize(0.44);
  gStyle->SetMarkerColor(2);
  gStyle->SetMarkerStyle(20);


  double chisq[nbinmx], chisq1[nbinmx];
  double chisq0;

  int len1 = strlen(argv[1]);
  strncpy(titley,argv[1],len1-5);
  titley[len1-5] = '\0';
  
  // sprintf(namey, "o_%s.txt", titley);
  // ofstream file_out(namey); // "outfilex.txt");
  sprintf(namey, "%s_o.ps", titley);
  TPostScript ps(namey, 111); // "outfilex.ps",111);  
  ps.Range(20,30); //ps.Range(10,20);
  ps.NewPage();  
  
  sprintf(namey, "%s_o.root", titley);
  TFile* fileOut = new TFile(namey, "recreate");

  TFile* fileIn = new TFile(argv[1]);

#ifdef isP
  TH2D* sim_response = (TH2D*)((TH2D*)fileIn->Get("sim_response_p"))->Clone();
#else
  TH2D* sim_response = (TH2D*)((TH2D*)fileIn->Get("sim_response_n"))->Clone();
#endif
  
  int nbinx = sim_response->GetNbinsX();
  int nbiny = sim_response->GetNbinsY();
  cout << " nbinx " << nbinx << " nbiny " << nbiny << endl;

#ifdef isP
  TH1D* momin_true = (TH1D*)((TH1D*)fileIn->Get("momin_sim_p"))->Clone();
#else
  TH1D* momin_true = (TH1D*)((TH1D*)fileIn->Get("momin_sim_n"))->Clone();
#endif
  momin_true->SetLineColor(1);
  momin_true->SetBinContent(0,0);
  momin_true->SetBinContent(nbiny+1,0);

#ifdef isP
  TH1D* sim_reco = (TH1D*)((TH1D*)fileIn->Get("sim_reco_p"))->Clone();
#else
  TH1D* sim_reco = (TH1D*)((TH1D*)fileIn->Get("sim_reco_n"))->Clone();
#endif
  sim_reco->SetBinContent(0,0);
  sim_reco->SetBinContent(nbinx+1,0);
  
#ifdef isP
  TH1D* data_reco = (TH1D*)((TH1D*)fileIn->Get("data_reco_p"))->Clone();
#else
  TH1D* data_reco = (TH1D*)((TH1D*)fileIn->Get("data_reco_n"))->Clone();
#endif
  // TH1D* data_reco = (TH1D*)((TH1D*)fileIn->Get("sim_reco_p"))->Clone();
  // data_reco->SetNameTitle("data_reco_p","data_reco_p");
  
  // data_reco->Scale(sim_reco->GetMaximum()/data_reco->GetMaximum());
  
  // data_reco->SetBinContent(0,0);
  // data_reco->SetBinContent(nbinx+1,0);
  // data_reco->SetBinError(0,0);
  // data_reco->SetBinError(nbinx+1,0);
  
  // for(int nx1=0; nx1<=nbinx+1; nx1++) {
  //   cout << " " << nx1 << " " << sim_reco->GetBinContent(nx1) << " " << data_reco->GetBinContent(nx1) << endl;
  // }
  
  // data_reco->SetBinContent(0,sim_reco->GetBinContent(0));
  // data_reco->SetBinContent(nbinx+1,sim_reco->GetBinContent(nbinx+1));
  // data_reco->SetBinError(0,sim_reco->GetBinError(0));
  // data_reco->SetBinError(nbinx+1,sim_reco->GetBinError(nbinx+1));
  
  TH1D *fakerate = (TH1D*)sim_reco->Clone();
  fakerate->SetNameTitle("fakerate","fakerate");
  fakerate->Add((TH1D*)sim_response->ProjectionX(),-1.);
  fakerate->Divide(sim_reco);
  // for(int ix=0; ix<nbinx; ix++) {
  //   data_reco->SetBinContent(ix+1, (1.0 - fakerate->GetBinContent(ix+1))*data_reco->GetBinContent(ix+1));
  //   data_reco->SetBinError(ix+1, (1.0 - fakerate->GetBinContent(ix+1))*data_reco->GetBinError(ix+1));
  // }
  
  // data_reco->SetBinContent(0, (1.0 - fakerate->GetBinContent(1))*data_reco->GetBinContent(0));
  // data_reco->SetBinError(0, (1.0 - fakerate->GetBinContent(1))*data_reco->GetBinError(0));
  // data_reco->SetBinContent(nbinx+1, (1.0 - fakerate->GetBinContent(nbinx))*data_reco->GetBinContent(nbinx+1));
  // data_reco->SetBinError(nbinx+1, (1.0 - fakerate->GetBinContent(nbinx))*data_reco->GetBinError(nbinx+1));
  
  // for(int nx1=0; nx1<=nbinx+1; nx1++) {
  //   cout << " " << nx1 << " " << sim_reco->GetBinContent(nx1) << " " << data_reco->GetBinContent(nx1) << endl;
  // }

  // for(int nx1=0; nx1<nbinx/2; nx1++) {
  //   // cout << " " << nx1 << " " << sim_reco->GetBinContent(nx1+1) << " " << data_reco->GetBinContent(nx1+1) << endl;
  //   data_reco->SetBinContent(nx1+1, sim_reco->GetBinContent(nx1+1));
  //   data_reco->SetBinError(nx1+1, sim_reco->GetBinError(nx1+1));
  // }
  
  
  chisq0 = 0;
  for(int nx1=0; nx1<nbinx; nx1++) {
    if(data_reco->GetBinContent(nx1+1)!=0
       && sim_reco->GetBinContent(nx1+1)!=0) {
      chisq0 += pow(data_reco->GetBinContent(nx1+1)-sim_reco->GetBinContent(nx1+1),2)/(pow(data_reco->GetBinError(nx1+1),2)+pow(sim_reco->GetBinError(nx1+1),2));
    }
  }
  cout<<" chisq0 "<<chisq0<<endl;
  
  TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
  c1->Divide(2,2,1.e-5, 1.e-5);
  c1->cd(1);
  c1->cd(2);
  data_reco->Draw("hist");
  c1->cd(3);
  momin_true->Draw();
  c1->cd(4);
  sim_response->Draw("colz");  
  c1->Update();
  ps.NewPage();  
  
  float nmc = momin_true->Integral();
  
  cout<<"Title 1 "<<endl;
  sprintf(namey, "unfold_%s", data_reco->GetName());
  sprintf(titley, "unfold_%s", data_reco->GetTitle());
  RooUnfoldResponse response = RooUnfoldResponse(sim_reco, momin_true, sim_response, namey, titley);
  // RooUnfoldResponse response = RooUnfoldResponse(0, 0, sim_response, namey, titley); // no fake, no effi
  // RooUnfoldResponse response = RooUnfoldResponse(0, momin_true, sim_response, namey, titley); // no fake
  // RooUnfoldResponse response = RooUnfoldResponse(sim_reco, 0, sim_response, namey, titley); // no effi
  cout<<"Title 2 "<<endl;
  
  TH1D* h_unfoldedbayes[100]={0};
  
  momin_true->SetLineWidth(3);
  int canSel=0;
  for(int ijk=0; ijk<nitermx; ijk++)  {
    if(canSel<4) {
      canSel++;
    } else {
      canSel = 1;
      c1->Update();
      ps.NewPage();
    }
    c1->cd(canSel);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    
    // RooUnfoldBinByBin unfoldbayes(&response,data_reco);
    // RooUnfoldSvd unfoldbayes(&response,data_reco,ijk+1);
    // RooUnfoldTUnfold unfoldbayes(&response,data_reco);
    // unfoldbayes.SetRegParm(ijk+1);
    // RooUnfoldInvert unfoldbayes(&response, data_reco);
    
    RooUnfoldBayes unfoldbayes(&response, data_reco, ijk+1, false);
    // RooUnfoldBayes unfoldbayes(&response, sim_reco, ijk+1, false);
    unfoldbayes.SetVerbose(0);
    
    h_unfoldedbayes[ijk] = (TH1D*)unfoldbayes.Hreco(RooUnfold::kCovariance)->Clone();
    // h_unfoldedbayes[ijk] = (TH1D*)unfoldbayes.Hreco(RooUnfold::kErrors)->Clone();
    // h_unfoldedbayes[ijk] = (TH1D*)unfoldbayes.Hreco(RooUnfold::kNoError)->Clone();
    // h_unfoldedbayes[ijk] = (TH1D*)unfoldbayes.Hreco(RooUnfold::kCovToy)->Clone();
    
    cout<<"Title 6 "<<ijk<<endl;
    char namez[100];
    
    sprintf(namez, "bayes_h_i%i", ijk);
    h_unfoldedbayes[ijk]->SetNameTitle(namez,namez);
    // h_unfoldedbayes[ijk]->Scale(nmc/max(1.,h_unfoldedbayes[ijk]->Integral()));

    h_unfoldedbayes[ijk]->SetLineColor(2);
    h_unfoldedbayes[ijk]->DrawNormalized("hist");
    
    momin_true->DrawNormalized("E1:same");
    
    chisq[ijk] = 0;
    chisq1[ijk] = 0;
    for(int nx1=0; nx1<nbinx; nx1++) {
      if(h_unfoldedbayes[ijk]->GetBinContent(nx1+1)!=0
	 && momin_true->GetBinContent(nx1+1)!=0) {
	chisq[ijk] += pow(h_unfoldedbayes[ijk]->GetBinContent(nx1+1)-momin_true->GetBinContent(nx1+1),2)/(pow(h_unfoldedbayes[ijk]->GetBinError(nx1+1),2)+pow(momin_true->GetBinError(nx1+1),2));
	chisq1[ijk] += pow(h_unfoldedbayes[ijk]->GetBinContent(nx1+1)-momin_true->GetBinContent(nx1+1),2);
      }
    } // for(int nx1=0; nx1<nbinx; nx1++) {
    
    fileOut->cd();
    h_unfoldedbayes[ijk]->Write();
  } // for(int ijk=0; ijk<nitermx; ijk++)  {
  c1->Update();
  
  ps.Close();

  double xval[nitermx];
  double yval1[nitermx];
  TGraph* gr1;
  int npts = 0;
  for(int jk=0; jk<nitermx+1; jk++) {
    xval[jk] = jk;
    yval1[jk] = chisq[jk];
    npts++;
  }
  sprintf(namey,"gr1_chi2");
  gr1 = new TGraph(npts,xval,yval1);
  gr1->SetName(namey);
  gr1->SetTitle(namey);
  gr1->SetMarkerSize(1);
  gr1->SetMarkerStyle(22);
  gr1->SetDrawOption("P");
  
  data_reco->Write();
  sim_reco->Write();
  momin_true->Write();
  sim_response->Write();
  fakerate->Write();
  gr1->Write();
  
  fileOut->Close();
  cout<<"Title 7 "<<endl;
  
  cout<<"0 ";
  cout<<"& "<<chisq0<<" ";
  cout<<"\\\\"<<endl;
  for(int ij=0; ij<nitermx; ij++) {
    cout<<ij+1<<" ";
    cout<<"& "<<chisq[ij]<<" ";
    cout<<"& "<<chisq1[ij]<<" ";
    cout<<"\\\\"<<endl;
  }
  
}

