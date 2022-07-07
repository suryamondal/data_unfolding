

void unfoldingCorsika_Half() {
  

  TFile *f_simh = new TFile("outputrootfiles/corsika76300_FLUKA_SIBYLL_3dFlux_20220105av_trg5of8_20220105aa_750Lakh_20211205aq_o.root","READ");
  TFile *f_sim  = new TFile("outputrootfiles/corsika76300_FLUKA_SIBYLL_3dFlux_20220105av_trg5of8_20220105aa_750Lakh_20211205aq_o.root","READ");
  TFile *f_data = new TFile("outputrootfiles/BRPCv4t_evtraw_201812_20211205aq_o.root","READ");

  const int nlayer = 10;
  
  const double maxPosDev = 1.;
  const double chi2ndf = 7.5;
  const double momerrdiv = 0.1;
  const int minHit = 6;

  double GasChamberSizeX = 2.*87;    // in cm
  double GasChamberSizeY = 2.*91.75; // in cm
  
  const double detBinWidth = 0.1;
  const double genBinWidth = 0.1;
  const double xminDet =  0.75;
  const double xmaxDet =  6.0;
  const double xminGen =  0.85;
  const double xmaxGen =  3.6;
  const int nDet = (xmaxDet-xminDet)/detBinWidth;
  const int nGen = (xmaxGen-xminGen)/genBinWidth;
  
  double detEdge[2][1000], genEdge[2][1000];
  for(int ij=0;ij<nDet+1;ij++) {
    // detEdge[0][ij] = pow(10.,log10(xminDet)+ij*(log10(xmaxDet)-log10(xminDet))/nDet);
    detEdge[0][ij] = xminDet+ij*(xmaxDet-xminDet)/nDet;
    cout << " det " << detEdge[0][ij] << endl;
  }
  for(int ij=0;ij<nGen+1;ij++) {
    // genEdge[0][ij] = pow(10.,log10(xminGen)+ij*(log10(xmaxDet)-log10(xminGen))/nGen);
    genEdge[0][ij] = xminGen+ij*(xmaxGen-xminGen)/nGen;
    cout << " gen " << genEdge[0][ij] << endl;
  }
  for(int ij=0;ij<nDet+1;ij++) {
    detEdge[1][ij] = -detEdge[0][nDet-ij];
    cout << " det " << detEdge[1][ij] << endl;
  }
  for(int ij=0;ij<nGen+1;ij++) {
    genEdge[1][ij] = -genEdge[0][nGen-ij];
    cout << " gen " << genEdge[1][ij] << endl;
  }
  
  TH1D *hres[2], *hsim[2], *hin[2], *hout[2], *hdata[2], *hinEff[2];
  TH2D *hinout[2];
    
  for(int ijh=0;ijh<2;ijh++) {
    hres[ijh]   = new TH1D(TString::Format("reco_res_%s",(ijh==0?"p":"n")),
			   TString::Format("reco_res_%s",(ijh==0?"p":"n")),
			   nGen,-8.,8.);
    hsim[ijh]   = new TH1D(TString::Format("momin_sim_%s",(ijh==0?"p":"n")),
			   TString::Format("momin_sim_%s with inefficiency",(ijh==0?"p":"n")),
			   nGen,genEdge[0]);
    hin[ijh]    = new TH1D(TString::Format("sim_gen_%s",(ijh==0?"p":"n")),
			   TString::Format("sim_gen_%s",(ijh==0?"p":"n")),
			   nGen,genEdge[0]);
    hout[ijh]   = new TH1D(TString::Format("sim_reco_%s",(ijh==0?"p":"n")),
			   TString::Format("sim_reco_%s with background",(ijh==0?"p":"n")),
			   nDet,detEdge[0]);
    hinout[ijh] = new TH2D(TString::Format("sim_response_%s",(ijh==0?"p":"n")),
			   TString::Format("sim_response_%s",(ijh==0?"p":"n")),
			   nDet,detEdge[0],
			   nGen,genEdge[0]);
    hdata[ijh]  = new TH1D(TString::Format("data_reco_%s",(ijh==0?"p":"n")),
			   TString::Format("data_reco_%s",(ijh==0?"p":"n")),
			   nDet,detEdge[0]);
  } // for(int ijh=0;ijh<2;ijh++) {
  
  TH2D *f_thephi_gen = new TH2D("thephi_gen","thephi_gen",90,-TMath::Pi(),TMath::Pi(),90,0.5*TMath::Pi(),TMath::Pi());
  TH2D *f_thephi_trg = new TH2D("thephi_trg","thephi_trg",90,-TMath::Pi(),TMath::Pi(),90,0.5*TMath::Pi(),TMath::Pi());
  TH2D *f_thephi_sr  = new TH2D("thephi_sr" ,"thephi_sr" ,90,-TMath::Pi(),TMath::Pi(),90,0.5*TMath::Pi(),TMath::Pi());
  
  TH1D *f_hchi2ndf_sim = new TH1D("chi2ndf_sim","chi2ndf_sim",200,0,10);
  TH1D *f_hchi2ndf_data = new TH1D("chi2ndf_data","chi2ndf_data",200,0,10);
  TH1D *f_hndf_sim = new TH1D("ndf_sim","ndf_sim",6,4.5,10.5);
  TH1D *f_hndf_data = new TH1D("ndf_data","ndf_data",6,4.5,10.5);
  TH1D *f_phiout_sim = new TH1D("phiout_sim","phiout_sim",180,-TMath::Pi(),TMath::Pi());
  TH1D *f_theout_sim = new TH1D("theout_sim","theout_sim",180,0,TMath::Pi());
  TH1D *f_phiout_data = new TH1D("phiout_data","phiout_data",180,-TMath::Pi(),TMath::Pi());
  TH1D *f_theout_data = new TH1D("theout_data","theout_data",180,0,TMath::Pi());

  TH1D *hsep = new TH1D("recosep","recosep",1000,0.,0.1);
  
  const double detBinWidth1 = 0.05;
  const double genBinWidth1 = 0.05;
  const double xminDet1 = -6.;
  const double xmaxDet1 =  6.;
  const double xminGen1 = -4;
  const double xmaxGen1 =  4.;
  const int nDet1 = (xmaxDet1-xminDet1)/detBinWidth1;
  const int nGen1 = (xmaxGen1-xminGen1)/genBinWidth1;
  
  TH1D *f_hsim = new TH1D("momin_sim","momin_sim",nGen1,xminGen1,xmaxGen1);
  TH1D *f_hin = new TH1D("sim_gen","sim_gen",nGen1,xminGen1,xmaxGen1);
  TH1D *f_hout = new TH1D("sim_reco","sim_reco",nDet1,xminDet1,xmaxDet1);
  TH2D *f_hinout = new TH2D("sim_response","sim_response",nDet1,xminDet1,xmaxDet1,nGen1,xminGen1,xmaxGen1);
  TH1D *f_hdata = new TH1D("data_reco","data_reco",nDet1,xminDet1,xmaxDet1);
  
  TH1D *f_hsim1[nlayer-minHit+1], *f_hin1[nlayer-minHit+1], *f_hout1[nlayer-minHit+1];
  TH2D *f_hinout1[nlayer-minHit+1];
  TH1D *f_hdata1[nlayer-minHit+1];
  for(int nl=minHit;nl<=nlayer;nl++) {
    f_hin1[nl-minHit] = new TH1D(TString::Format("sim_gen_nl%i",nl),
				 TString::Format("sim_gen_nl%i",nl),
				 nGen1,xminGen1,xmaxGen1);
    f_hout1[nl-minHit] = new TH1D(TString::Format("sim_reco_nl%i",nl),
				  TString::Format("sim_reco_nl%i",nl),
				  nDet1,xminDet1,xmaxDet1);
    f_hinout1[nl-minHit] = new TH2D(TString::Format("sim_response_nl%i",nl),
				    TString::Format("sim_response_nl%i",nl),
				    nDet1,xminDet1,xmaxDet1,nGen1,xminGen1,xmaxGen1);
    f_hdata1[nl-minHit] = new TH1D(TString::Format("data_reco_nl%i",nl),
				   TString::Format("data_reco_nl%i",nl),
				   nDet1,xminDet1,xmaxDet1);
  }
  
  const double lenBinW = 0.1;	// in m
  const int lenBinC = 5./lenBinW;
  TH1D *f_hsim2[lenBinC], *f_hin2[lenBinC], *f_hout2[lenBinC];
  TH2D *f_hinout2[lenBinC];
  TH1D *f_hdata2[lenBinC];
  for(int nl=0;nl<lenBinC;nl++) {
    f_hin2[nl] = new TH1D(TString::Format("sim_gen_len%i",nl),
			  TString::Format("sim_gen_len%i",nl),
			  nGen1,xminGen1,xmaxGen1);
    f_hout2[nl] = new TH1D(TString::Format("sim_reco_len%i",nl),
			   TString::Format("sim_reco_len%i",nl),
			   nDet1,xminDet1,xmaxDet1);
    f_hinout2[nl] = new TH2D(TString::Format("sim_response_len%i",nl),
			     TString::Format("sim_response_len%i",nl),
			     nDet1,xminDet1,xmaxDet1,nGen1,xminGen1,xmaxGen1);
    f_hdata2[nl] = new TH1D(TString::Format("data_reco_len%i",nl),
			    TString::Format("data_reco_len%i",nl),
			    nDet1,xminDet1,xmaxDet1);
  }
  
  Int_t nhits_finder;
  Double_t momin, momout, momerr, chi2, theout, phiout, thein, phiin, recosep;
  Double_t trkLen;

  f_simh->cd();
  
  TTree *MomTreeh = (TTree*)f_simh->Get("SimTree");
  MomTreeh->SetBranchAddress("momin",&momin);
  MomTreeh->SetBranchAddress("thein",&thein);
  MomTreeh->SetBranchAddress("phiin",&phiin);
  Long64_t nentriesh = MomTreeh->GetEntries();
  for(Long64_t ij=0;ij<nentriesh;ij++) {
    if(ij%1000000==0) {
      cout << " reading sim1 " << nentriesh-ij << endl;}
    MomTreeh->GetEntry(ij);
    f_hsim->Fill(momin);
    if(fabs(momin)>2.) {
      f_thephi_gen->Fill(phiin,thein);}
    for(int ijh=0;ijh<2;ijh++) {
      if(ijh==1) {momin*=-1;}
      hsim[ijh]->Fill(momin);
    } // for(int ijh=0;ijh<2;ijh++) {
  }   // for(Long64_t ij=0;ij<nentriesh;ij++) {

  f_simh->cd();
  
  TTree *MomTreet = (TTree*)f_simh->Get("TriggerTree");
  MomTreet->SetBranchAddress("momin",&momin);
  MomTreet->SetBranchAddress("thein",&thein);
  MomTreet->SetBranchAddress("phiin",&phiin);
  Long64_t nentriest = MomTreet->GetEntries();
  for(Long64_t ij=0;ij<nentriest;ij++) {
    if(ij%1000000==0) {
      cout << " reading trg " << nentriest-ij << endl;}
    MomTreet->GetEntry(ij);
    if(fabs(momin)>2.) {
      f_thephi_trg->Fill(phiin,thein);
      f_thephi_sr->Fill(phiin,thein);}
  }   // for(Long64_t ij=0;ij<nentriest;ij++) {

  for(int ijx=0;ijx<f_thephi_sr->GetNbinsX();ijx++) {
    for(int ijy=0;ijy<f_thephi_sr->GetNbinsY();ijy++) {
      double srval  = f_thephi_sr->GetBinContent(ijx+1,ijy+1);
      double genval = f_thephi_gen->GetBinContent(ijx+1,ijy+1);
      double tmpval = genval==0?0:srval/genval;
      double tmperr = genval==0?0:tmpval*(1./sqrt(srval)+1./sqrt(genval));
      f_thephi_sr->SetBinContent(ijx+1,ijy+1,tmpval);
      f_thephi_sr->SetBinError(ijx+1,ijy+1,tmperr);
    }}

  double sum_sr = 0, sum_sr_tst = 0;
  double sum_sr_err = 0;
  for(int ijx=0;ijx<f_thephi_sr->GetNbinsX();ijx++) {
    for(int ijy=0;ijy<f_thephi_sr->GetNbinsY();ijy++) {
      double ymin = f_thephi_sr->GetYaxis()->GetBinLowEdge(ijy+1);
      double ymax = f_thephi_sr->GetYaxis()->GetBinUpEdge(ijy+1);
      double xbnw = f_thephi_sr->GetXaxis()->GetBinWidth(ijx+1);
      double tmpintval = xbnw*(TMath::Cos(ymin)-TMath::Cos(ymax));
      sum_sr_tst += tmpintval;
      double tmpval = f_thephi_sr->GetBinContent(ijx+1,ijy+1)*tmpintval;
      sum_sr += tmpval;
      if(tmpval!=0) {sum_sr_err += tmpintval*f_thephi_sr->GetBinError(ijx+1,ijy+1);}
    }}
  cout<<" sum_sr "<<sum_sr<<" sum_sr_tst "<<sum_sr_tst<<endl;
  cout<<" sum_sr_err "<<sum_sr_err<<endl;
  cout<<" sum_sr_err wted "<<sum_sr_err/f_thephi_sr->GetSumOfWeights()<<endl;
  
  f_sim->cd();
  
  TTree *MomTree = (TTree*)f_sim->Get("MomTree");
  MomTree->SetBranchAddress("momin",&momin);
  MomTree->SetBranchAddress("nhits_finder",&nhits_finder);
  MomTree->SetBranchAddress("chi2",&chi2);
  MomTree->SetBranchAddress("trkLen",&trkLen);
  MomTree->SetBranchAddress("momout",&momout);
  MomTree->SetBranchAddress("momerr",&momerr);
  MomTree->SetBranchAddress("theout",&theout);
  MomTree->SetBranchAddress("phiout",&phiout);
    
  int sentries = MomTree->GetEntries();
  for(int ij=0;ij<sentries;ij++) {
    MomTree->GetEntry(ij);
    if(ij%1000000==0) {
      cout << " reading sim " << sentries-ij << endl;
    }
    
    if(chi2>0) {
      f_hchi2ndf_sim->Fill(chi2/(nhits_finder-5.));
      f_hndf_sim->Fill(nhits_finder);
      f_theout_sim->Fill(theout);
      f_phiout_sim->Fill(phiout);
    }
    
    if(1
       && nhits_finder>=minHit
       && chi2/(nhits_finder-5.)<chi2ndf
       ) {
      
      f_hin->Fill(momin);
      f_hout->Fill(momout);
      f_hin1[nhits_finder-minHit]->Fill(momin);
      f_hout1[nhits_finder-minHit]->Fill(momout);
      if(trkLen>0 && trkLen<5.) {
	f_hin2[int(trkLen/lenBinW)]->Fill(momin);
	f_hout2[int(trkLen/lenBinW)]->Fill(momout);}
      if(momin>xminGen1 && momin<xmaxGen1
	 && momout>xminDet1 && momout<xmaxDet1) {
	f_hinout->Fill(momout,momin);
	f_hinout1[nhits_finder-minHit]->Fill(momout,momin);
	if(trkLen>0 && trkLen<5.) {
	  f_hinout2[int(trkLen/lenBinW)]->Fill(momout,momin);}
      }
      
      for(int ijh=0;ijh<2;ijh++) {
	if(ijh==1) {momin*=-1; momout*=-1;}
	if(momin>xminGen && momin<xmaxGen) {
	  hres[ijh]->Fill((momin-momout)/momin);
	}
	if(1) {
	  hin[ijh]->Fill(momin);
	  hout[ijh]->Fill(momout);
	  if(momin>xminGen && momin<xmaxGen
	     && momout>xminDet && momout<xmaxDet) {
	    hinout[ijh]->Fill(momout,momin);
	  }
	}
      }	// for(int ijh=0;ijh<2;ijh++) {
    }
  } // for(int ij=0;ij<sentries;ij++) {

  
  f_data->cd();
  
  TTree *MomTreed = (TTree*)f_data->Get("MomTree");
  MomTreed->SetBranchAddress("nhits_finder",&nhits_finder);
  MomTreed->SetBranchAddress("chi2",&chi2);
  MomTreed->SetBranchAddress("trkLen",&trkLen);
  MomTreed->SetBranchAddress("momout",&momout);
  MomTreed->SetBranchAddress("momerr",&momerr);
  MomTreed->SetBranchAddress("theout",&theout);
  MomTreed->SetBranchAddress("phiout",&phiout);
  MomTreed->SetBranchAddress("recosep",&recosep);
  
  int dentries = MomTreed->GetEntries();
  for(int ij=0;ij<dentries;ij++) {
    MomTreed->GetEntry(ij);
    if(ij%1000000==0) {
      cout << " reading data " << dentries-ij << endl;
    }
    if(recosep>0) {hsep->Fill(recosep);}
    
    if(chi2>0) {
      f_hchi2ndf_data->Fill(chi2/(nhits_finder-5.));
      f_hndf_data->Fill(nhits_finder);
      f_theout_data->Fill(theout);
      f_phiout_data->Fill(phiout);
    }
    
    if(1
       && nhits_finder>=minHit
       && chi2/(nhits_finder-5.)<chi2ndf
       ) {
      
      f_hdata->Fill(momout);
      f_hdata1[nhits_finder-minHit]->Fill(momout);
      if(trkLen>0 && trkLen<5.) {
	f_hdata2[int(trkLen/lenBinW)]->Fill(momout);}
      
      for(int ijh=0;ijh<2;ijh++) {
	if(ijh==1) {momout*=-1;}
	hdata[ijh]->Fill(momout);
      }	// for(int ijh=0;ijh<2;ijh++) {
    }
  } // for(int ij=0;ij<dentries;ij++) {
  
  TF1 *fe1 = new TF1("fe1","TMath::Exp([0]+[1]*x)",0.,1.);
  hsep->Fit("fe1","","",0.5e-3,0.1);
  double slope = fabs(fe1->GetParameter(1));
  cout<<" slope "<<slope<<endl;
  double total0 = fe1->Integral(0.,1.)/hsep->GetXaxis()->GetBinWidth(1);
  double total1 = fe1->Integral(0.5e-3,1.)/hsep->GetXaxis()->GetBinWidth(1);
  double totald = fe1->Integral(0.,0.5e-3)/hsep->GetXaxis()->GetBinWidth(1);
  cout<<" total0 "<<total0<<" total1 "<<total1<<endl;
  double totalT = total0/slope;
  double totalD = totald/slope;
  cout<<" totalT "<<totalT/(3600.*24.)<<endl;
  cout<<" totalD "<<totalD/(3600.*24.)<<" "<<totalD/totalT<<endl;
  
  
  TFile *file = new TFile("unfolding/data/coorelation_half.root","recreate");
  file->cd();
  hsep->Write();
  f_thephi_gen->Write();
  f_thephi_trg->Write();
  f_thephi_sr->Write();
  f_hchi2ndf_sim->Write();
  f_hchi2ndf_data->Write();
  f_hndf_sim->Write();
  f_hndf_data->Write();
  f_theout_sim->Write();
  f_phiout_sim->Write();
  f_theout_data->Write();
  f_phiout_data->Write();
  for(int ijh=0;ijh<2;ijh++) {
    hsim[ijh]->Write();
    hin[ijh]->Write();
    hout[ijh]->Write();
    hinout[ijh]->Write();
    hdata[ijh]->Write();
  } // for(int ijh=0;ijh<2;ijh++) {
  f_hsim->Write();
  f_hin->Write();
  f_hout->Write();
  f_hinout->Write();
  f_hdata->Write();
  for(int nl=minHit;nl<=nlayer;nl++) {
    f_hin1[nl-minHit]->Write();
    f_hout1[nl-minHit]->Write();
    f_hinout1[nl-minHit]->Write();
    f_hdata1[nl-minHit]->Write();
  }
  for(int nl=0;nl<lenBinC;nl++) {
    f_hin2[nl]->Write();
    f_hout2[nl]->Write();
    f_hinout2[nl]->Write();
    f_hdata2[nl]->Write();
  }
    
  file->Close();
  f_data->Close();
  f_sim->Close();
  f_simh->Close();

} // unfoldingCorsika_Half() {
