/*

source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.00/x86_64-slc5-gcc43-dbg/root/bin/thisroot.csh

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/physics/RooUnfold-1.1.1

RooUnfoldResponse response(h_recoinput[indexx], h_mcgeninput[indexx], h2d_mccorrel, namey, titley);

rm anal_evt_unfold

make anal_evt_unfold

./anal_evt_unfold | tee test_aap8cor_alpgen.log

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


bool pedsub=true;
using namespace std;
// using namespace CLHEP;

const char* namex;
char namey[300];
const char* titlex;
char titley[300];

const int nbinmx=1000;
const int nitermx=10;
int icol[12] = {2,3,3,4,4,5,5,6,6,7,7,1};


// int subtract_background(TH2D* h2d_correl, TH1D* data, TH1D* mc, TH1D* mc2, TH1D* mc3, TH1D* mc4, TH1D* mc5, double* fakerate, double* effi) {
int subtract_background(TH2D* h2d_correl, TH1D* mcd, double* fakerate, double* effi) {
  int nbinx = h2d_correl->GetNbinsX();
  int nbiny = h2d_correl->GetNbinsY();
  if (nbinx>nbinmx || nbiny > nbinmx) {
    cout <<"Increase nbinmx, which is "<< nbinmx<<" in subtract_background"<<endl;
    cout <<"Reconstructed objects are not pedestal subtracted"<<endl;
    return 0;
  }
  double totalgen[nbinmx]={0.};
  double totalreco[nbinmx]={0.};
  for (int ix=0; ix<nbinx+1; ix++) {
    for (int iy=0; iy<nbiny+1; iy++) {
      if (ix==nbinx && iy==nbiny) continue;
      totalreco[ix] +=h2d_correl->GetBinContent(ix+1, iy+1);
      if (iy==nbiny) fakerate[ix] =h2d_correl->GetBinContent(ix+1, iy+1);
      
      totalgen[iy] +=h2d_correl->GetBinContent(ix+1, iy+1);
      if (ix==nbinx) effi[iy] =h2d_correl->GetBinContent(ix+1, iy+1);
      
      // cout <<"ix "<< ix<<" iy "<<iy
      // 	   <<" binc "<<h2d_correl->GetBinContent(ix+1, iy+1)
      // 	   <<" gen "<<totalgen[iy]<<" fake "<<fakerate[ix]
      // 	   <<" reco "<<totalreco[ix]<<" effi "<<effi[iy]<<endl;

      //GMA 27th Jan 2013
      if (ix==nbinx || iy==nbiny) {
      	h2d_correl->SetBinContent(ix+1, iy+1, 0.0);
      	h2d_correl->SetBinError(ix+1, iy+1, 0.0);
      }
    }

    // SNM
    // for (int iy=0; iy<nbiny; iy++) {
    //   effi[iy] = (totalgen[iy] - effi[iy])/max(1.0, totalgen[iy]);
    //   if (effi[iy]<1.e-3) effi[iy]=1.e-3;
    // }
    // for (int ix=0; ix<nbinx; ix++){
    //   // cout <<"fake "<< ix<<" "<< fakerate[ix]<<" "<<totalreco[ix]<<" "<<mc->GetBinContent(ix+1)<<" "<<mc2->GetBinContent(ix+1)<<" "<<mc->GetBinError(ix+1)<<" "<<mc2->GetBinError(ix+1)<<endl;
      
    //   fakerate[ix] /=max(1.0, totalreco[ix]);
    //   // cout << " " << ix << " " << fakerate[ix] << endl;
    // }
    
  }

  // SNM
  for (int iy=0; iy<nbiny; iy++) {
    effi[iy] = (totalgen[iy] - effi[iy])/max(1.0, totalgen[iy]);
    if (effi[iy]<1.e-3) effi[iy]=1.e-3;
  }
  for (int ix=0; ix<nbinx; ix++){
    fakerate[ix] /=max(1.0, totalreco[ix]);
    // cout << " " << ix << " " << fakerate[ix] << endl;
  }

  if(pedsub) { 
    for(int ix=0; ix<nbinx; ix++) {
      // mcd->SetBinContent(ix+1, (1.0 - fakerate[ix])*mcd->GetBinContent(ix+1));
      // mcd->SetBinError(ix+1, (1.0 - fakerate[ix])*mcd->GetBinError(ix+1));
    }
  }
  return 1;
}

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


  char mcrootfile[100];
  char datarootfile[100];
  char hist_mc[100];
  char hist_reco[100];
  char hist_unfol[100];


  //24th Jan 2013
  //Number of iterations/regularisation paramters in bayes and SVD 
  double dataerr[100]; //be sure that 
  double datval[100]; 
  double mxdaterr=0.0;
  double relsumerr=0.0;
  double relsum2err=0.0;

  double chisq[1000];
  double chisq0;
  double delchisq[1000];
  //  const char* ttlx[4]={"mc","mc2", "mc3", "data"};
  int    mxdatbin=-1;
  //  TCanvas *c1;
  //  char namey[100];
  //  char titley[100];
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
  
  TH2D* trkmm_mc_vs_momin = (TH2D*)((TH2D*)fileIn->Get("sim_response"))->Clone();
  // TH2D* trkmm_mc_vs_momin = (TH2D*)((TH2D*)fileIn->Get("sim_response_p"))->Clone();
  
  int nbinx = trkmm_mc_vs_momin->GetNbinsX();
  int nbiny = trkmm_mc_vs_momin->GetNbinsY();
  cout << " nbinx " << nbinx << " nbiny " << nbiny << endl;

  TH1D* momin_true = (TH1D*)((TH1D*)fileIn->Get("momin_true"))->Clone();
  // TH1D* momin_true = (TH1D*)((TH1D*)fileIn->Get("momin_sim_p"))->Clone();
  momin_true->SetLineColor(2);
  
  TH1D* sim_reco = (TH1D*)((TH1D*)fileIn->Get("sim_reco"))->Clone();
  // TH1D* sim_reco = (TH1D*)((TH1D*)fileIn->Get("sim_reco_p"))->Clone();
  
  // TH1D* data_gen = (TH1D*)((TH1D*)fileIn->Get("data_gen"))->Clone();
  
  TH1D* data_reco = (TH1D*)((TH1D*)fileIn->Get("data_reco"))->Clone();
  // TH1D* data_reco = (TH1D*)((TH1D*)fileIn->Get("data_reco_p"))->Clone();
  data_reco->Scale(sim_reco->GetMaximum()/data_reco->GetMaximum());
  
  
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
  data_reco->Draw();
  c1->cd(3);
  momin_true->Draw();
  c1->cd(4);
  trkmm_mc_vs_momin->Draw("colz");  
  c1->Update();
  ps.NewPage();  
  double fakerate[1000]={0.};     //This array size should be larger than histograme bin size
  double efficiency[1000]={0.};
  
  // int isub = subtract_background(trkmm_mc_vs_momin,data_reco,fakerate, efficiency);
  
  float nmc = momin_true->Integral();
  
  cout<<"Title 1 "<<endl;
  sprintf(namey, "unfold_%s", data_reco->GetName());
  sprintf(titley, "unfold_%s", data_reco->GetTitle());
  RooUnfoldResponse response = RooUnfoldResponse(sim_reco, momin_true, trkmm_mc_vs_momin, namey, titley);
  // RooUnfoldResponse response = RooUnfoldResponse(0, 0, trkmm_mc_vs_momin, namey, titley);
  // RooUnfoldResponse response = RooUnfoldResponse(0, momin_true, trkmm_mc_vs_momin, namey, titley);
  // RooUnfoldResponse response = RooUnfoldResponse(sim_reco, 0, trkmm_mc_vs_momin, namey, titley);
  cout<<"Title 2 "<<endl;
  
  TH1D* h_unfoldedbayes[100]={0};
  TH1D* cp_unfoldedbayes[100]={0};
  
  int nsvd = 15;
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
    momin_true->DrawNormalized("E1");

    cout << "\t before" << endl;

    // RooUnfoldBinByBin unfoldbayes(&response,data_reco);
    // RooUnfoldSvd unfoldbayes(&response,data_reco);
    // RooUnfoldTUnfold unfoldbayes(&response,data_reco);
    
    RooUnfoldBayes unfoldbayes(&response, data_reco, ijk+1, false);
    // RooUnfoldBayes unfoldbayes(&response, sim_reco, ijk+1, false);
    unfoldbayes.SetVerbose(0);

    cout << "\t middle" << endl;

    h_unfoldedbayes[ijk] = (TH1D*)unfoldbayes.Hreco(RooUnfold::kCovariance)->Clone();
    // h_unfoldedbayes[ijk] = (TH1D*)unfoldbayes.Hreco(RooUnfold::kErrors)->Clone();
    // h_unfoldedbayes[ijk] = (TH1D*)unfoldbayes.Hreco(RooUnfold::kNoError)->Clone();
    // h_unfoldedbayes[ijk] = (TH1D*)unfoldbayes.Hreco(RooUnfold::kCovToy)->Clone();

    cout << "\t end" << endl;
    
    cout<<"Title 6 "<<ijk<<endl;
    char namez[100];
    
    sprintf(namez, "bayes_h_i%i", ijk);
    h_unfoldedbayes[ijk]->SetNameTitle(namez,namez);
    h_unfoldedbayes[ijk]->Scale(nmc/max(1.,h_unfoldedbayes[ijk]->Integral())); 
    h_unfoldedbayes[ijk]->DrawNormalized("sames:hist");
    
    chisq[ijk] = 0;
    for(int nx1=0; nx1<nbinx; nx1++) {
      if(h_unfoldedbayes[ijk]->GetBinContent(nx1+1)!=0
      	 && momin_true->GetBinContent(nx1+1)!=0) {
      	chisq[ijk] += pow(h_unfoldedbayes[ijk]->GetBinContent(nx1+1)-momin_true->GetBinContent(nx1+1),2)/(pow(h_unfoldedbayes[ijk]->GetBinError(nx1+1),2)+pow(momin_true->GetBinError(nx1+1),2));
      }
      // if(h_unfoldedbayes[ijk]->GetBinContent(nx1+1)!=0
      // 	 && data_gen->GetBinContent(nx1+1)!=0) {
      // 	chisq[ijk] += pow(h_unfoldedbayes[ijk]->GetBinContent(nx1+1)-data_gen->GetBinContent(nx1+1),2)/(pow(h_unfoldedbayes[ijk]->GetBinError(nx1+1),2)+pow(data_gen->GetBinError(nx1+1),2));
      // }
      // if(h_unfoldedbayes[ijk]->GetBinContent(nx1+1)!=0
      // 	 && sim_reco->GetBinContent(nx1+1)!=0) {
      // 	chisq[ijk] += pow(h_unfoldedbayes[ijk]->GetBinContent(nx1+1)-sim_reco->GetBinContent(nx1+1),2)/(pow(h_unfoldedbayes[ijk]->GetBinError(nx1+1),2)+pow(sim_reco->GetBinError(nx1+1),2));
      // }
    } // for(int nx1=0; nx1<nbinx; nx1++) {
    fileOut->cd();
    h_unfoldedbayes[ijk]->Write();
  } // for(int ijk=0; ijk<nitermx; ijk++)  {
  c1->Update();

  ps.Close();
  
  double xval[nitermx];
  double yval1[nitermx];
  double yval2[nitermx];
  TGraph* gr1[2];

  TMultiGraph* mgr1 = new TMultiGraph("mgr1_chi2","#chisq^(2)");
  TMultiGraph* mgr2 = new TMultiGraph("mgr2_dchi2","#Delta #chisq^(2)");
  int npts = 0;
  for(int jk=0; jk<nitermx+1; jk++) {
    xval[jk] = jk;
    yval1[jk] = chisq[jk];
    yval2[jk] = chisq[jk] - chisq0;
    delchisq[jk] = yval2[jk];
    npts++;
  }
  sprintf(namey,"gr1_chi2");
  gr1[0] = new TGraph(npts,xval,yval1);
  gr1[0]->SetName(namey);
  gr1[0]->SetTitle(namey);
  gr1[0]->SetMarkerSize(1);
  gr1[0]->SetMarkerStyle(22);
  gr1[0]->SetDrawOption("P");

  sprintf(namey,"gr2_dchi2");
  gr1[1] = new TGraph(npts,xval,yval2);
  gr1[1]->SetName(namey);
  gr1[1]->SetTitle(namey);
  gr1[1]->SetMarkerSize(1);
  gr1[1]->SetMarkerStyle(22);
  gr1[1]->SetDrawOption("P");
  
  // ps.NewPage();
  // TCanvas* c2 = new TCanvas("c2", "c2", 600, 800);
  // c2->Divide(1,2);//,1.e-5, 1.e-5);
  // c2->cd(1);
  // mgr1->Draw("AP");
  // c2->cd(2);
  // mgr2->Draw("AP");
  
  // c2->Update();
  // ps.Close();

  // effi_clo->Write();
  sim_reco->Write();
  data_reco->Write();
  // data_gen->Write();
  momin_true->Write();
  trkmm_mc_vs_momin->Write();
  
  gr1[0]->Write();
  gr1[1]->Write();
  
  mgr1->Write();
  mgr2->Write();
  
  //    fileOut->Write();
  fileOut->Close();
  cout<<"Title 7 "<<endl;
  
  cout<<"0 ";
  cout<<"& "<<chisq0<<" ";
  cout<<"\\\\"<<endl;
  for(int ij=0; ij<nitermx; ij++) {
    cout<<ij+1<<" ";
    cout<<"& "<<chisq[ij]<<" ";
    cout<<"\\\\"<<endl;
  }
  
}

