/*
  Simple macro showing how to access branches from the delphes output root file,
  loop over events, and plot simple quantities such as the jet pt and the di-elect
  ron invariant
  mass.

  root -l examples/Example1.C'("delphes_output.root")'
*/


#ifdef __CLING__
//R__ADD_INCLUDE_PATH(/usr/local/src/Delphes-3.3.2)
//R__ADD_INCLUDE_PATH(/usr/local/src/Delphes-3.3.2/external)
R__ADD_INCLUDE_PATH(/usr/local/src/delphes)
R__ADD_INCLUDE_PATH(/usr/local/src/delphes/external)
R__LOAD_LIBRARY(libDelphes)
//#include "classes/DelphesClasses.h"
///#include "external/ExRootAnalysis/ExRootTreeReader.h"
//##include "/usr/local/src/Delphes-3.3.2/external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "THStack.h"
#include "TText.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/SortableObject.h"
//#include "/usr/local/src/Delphes-3.3.2/classes/DelphesClasses.h"
//#include "/usr/local/src/Delphes-3.3.2/external/ExRootAnalysis/ExRootTreeReader.h"
#include "/cvmfs/sft.cern.ch/lcg/releases/delphes/3.4.2pre12-0b66b/x86_64-slc6-gcc8-opt/include/classes/DelphesClasses.h"
#include "/cvmfs/sft.cern.ch/lcg/releases/delphes/3.4.2pre12-0b66b/x86_64-slc6-gcc8-opt/include/ExRootAnalysis/ExRootTreeReader.h"
#include "TF1.h"
#include <iostream>

void histadd(std::vector<TObject *> & hlist,std::vector<TObject *> & hadd);
void showcut(TH1F * plots[2],TString var,TString cut);
void showsb(TH1D * plots[3],TString var,double gam);
void subsb(TH1D * plots[3],TString var,double gam,double mhp);
void makeplots();
double GetCosThetaStar(TLorentzVector & g1,TLorentzVector & gg);
void analyse(const TString inputFiles[3],double gam,double mhp);
TH1D *  mysmooth(TH1D * hist);

bool truth=true;

//------------------------------------------------------------------------------
int main(){
  makeplots();
}

void makeplots() {

  //gStyle->SetOptStat(0);
  
  TString files[3];
  //  int type = 4;
  TString wid="_1";
  double gam = 2.4952;
  double mhp = 91.2;
  //files contains both, signal, bckgrd
  files[0] = "/mercury/data2/ben/SHERPA-MC-2.2.7/bfuller/subPlotM/bothOut.root";
  files[1] = "/mercury/data2/ben/SHERPA-MC-2.2.7/bfuller/subPlotM/bckgrdOut.root";
  files[2] = "/mercury/data2/ben/SHERPA-MC-2.2.7/bfuller/subPlotM/signalOut.root";
  analyse(files,gam,mhp);

}


void analyse(const TString inputFiles[3],double gam,double mhp)
{
  const double lumi = 13200;
  double jetPTCut = 25;
  double jetEtaCut = 4.9;
  double bEtaCut = 2.5;
  TString extn[3] = {"both","bckgrd","signal"};
  double sigBand1 = mhp -2*gam;
  double sigBand2 = mhp +2*gam;
  double lowband1 = mhp -10*gam;
  double lowband2 = mhp -3*gam;
  double highband1 = mhp +4*gam;
  double highband2 = mhp +20*gam;
  double histm1 = mhp - 8*gam;
  double histm2 = mhp + 8*gam;

  // Book histograms
  double x2 = 1000;
  int nb = 50;
  TH2F *h_akt10vakt4 = new TH2F("akt10vakt4", "akt10vakt4", nb,0.,x2, nb,0.,x2);
  TH1F * h_htcpt4[2];
  TH1F * h_htcpt10[2];
  for (int i=0;i<2;++i) {
    TString ext="fail";
    if (i==1) ext="pass";
    h_htcpt4[i] =new TH1F("ht "+ext+" pt4","ht "+ext+" pt4",nb,0.,2*x2);
    h_htcpt10[i] =new TH1F("ht "+ext+" pt10","ht "+ext+" pt10",nb,0.,2*x2);
  }
  TH1D *  h_mass[3];
  TH1D *  h_mass_cths[3];
  TH1D *  h_mass_tot[3];
  TH1D *  h_mass2[3];
  TH1D *  h_mass2_cths[3];
  TH1D *  h_mass2_tot[3];
  TH1D *  h_minbpt[3];
  TH1D *  h_maxbpt[3];
  TH1D *  h_nb[3];
  TH1D *  h_side_cths[3];
  TH1D *  h_sig_cths[3];
  for (int isb=0;isb<3;++isb) {
    h_mass[isb] =new TH1D("mass "+extn[isb],"mass"+extn[isb],100,0,1000);
    h_mass_cths[isb] =new TH1D("mass_cths"+extn[isb],"mass_cths"+extn[isb],100,0,1000);
    h_mass_tot[isb] =new TH1D("mass_tot"+extn[isb],"mass_tot"+extn[isb],100,0,1000);
    h_mass2[isb] =new TH1D("mass2"+extn[isb],"mass2"+extn[isb],100,histm1,histm2);
    std::cout<<"Booking h_mass2_cths["<<isb<<"] with min="<<histm1<<" max="<<histm2<<std::endl;
    h_mass2_cths[isb] =new TH1D("mass2_cths"+extn[isb],"mass2_cths"+extn[isb],100,histm1,histm2);
    std::cout<<" Adress of mass2_cths is "<<h_mass2_cths[isb]<<std::endl;
    h_mass2_tot[isb] =new TH1D("mass2_tot"+extn[isb],"mass2_tot"+extn[isb],50,histm1,histm2);
    h_minbpt[isb] =new TH1D("min b pT"+extn[isb],"min b pT "+extn[isb]+"; b p_{T} [GeV]",100,0,600);
    h_maxbpt[isb] =new TH1D("max b pT"+extn[isb],"max b pT "+extn[isb]+"; b p_{T} [GeV]",100,0,600);
    h_nb[isb] =new TH1D("nb "+extn[isb],"nb "+extn[isb],10,-0.5,9.5);
    h_side_cths[isb] =new TH1D("side cths "+extn[isb],"side cths "+extn[isb]+" ; cos #theta^{*}",20,-1.0,1.0);
    h_sig_cths[isb] =new TH1D("sig cths "+extn[isb],"sig cths "+extn[isb]+" ; cos #theta^{*}",20,-1.0,1.0);
  }
  double CrossSection[3] = {0.14,0.121,0.0217};
  for (int isb=0;isb<3;++isb) {  
    // Create chain of root trees
    TChain chain("Delphes");
    chain.Add(inputFiles[isb]);

    // Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();

    // Get pointers to branches used in this analysis
    TClonesArray *branchEvent = treeReader->UseBranch("Event");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  
    // Loop over all events
    std::cout<<"Will loop through "<<numberOfEntries<<" events for "<<extn[isb]<<std::endl;
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
      std::vector <TLorentzVector> blist;
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
	
      HepMCEvent   *event = (HepMCEvent*) branchEvent->At(0);
      double weight = CrossSection[isb]/numberOfEntries;
      int nchar = 0;

      double ht = 0;
      double maxbPT = 0;
      double minbPT = 1E9;
      int nb = 0;
      for(Int_t ip = 0; ip < branchParticle->GetEntries(); ++ip){
	GenParticle  *par = (GenParticle*) branchParticle->At(ip);
	if (par->Status != 1)  continue;
	if (entry < -5) {
	  std::cout<<"isb="<<isb<<" entry="<<entry<<" ip="<<ip<<" pid="<<par->PID<<"   pt="<<par->PT<<" status="<<par->Status<<std::endl;
	}

	if (fabs(par->Eta) > jetEtaCut) continue;
	ht += par->PT;
	if (fabs(par->Eta) > bEtaCut) continue;
	if (abs(par->PID) == 5) {
          nb++;
	  if (par->PT > maxbPT) maxbPT = par->PT;
	  if (par->PT < minbPT) minbPT = par->PT;
	}
      }
      h_minbpt[isb]->Fill(minbPT,weight);
      h_maxbpt[isb]->Fill(maxbPT,weight);
      h_nb[isb]->Fill(nb,weight);

      // find H+ candidate
      for(Int_t ip = 0; ip < branchParticle->GetEntries(); ++ip){
	GenParticle  *par = (GenParticle*) branchParticle->At(ip);
	if (par->Status != 1)  continue;
	//		std::cout<<"PID="<<par->PID<<std::endl;
	if (abs(par->PID) == 5)	  {
	  //	  std::cout<<"Got 6, PID="<<par->PID<<std::endl;

	  for(Int_t ipb = 0; ipb < branchParticle->GetEntries(); ++ipb){
	    GenParticle  *parb = (GenParticle*) branchParticle->At(ipb);
	    if (abs(parb->PID) == 5)   {
	      if (fabs(parb->Eta) > bEtaCut) continue;
	      if (parb->PT < jetPTCut) continue;

	      // Check we have a b of same charge
	      //   	      std::cout<<"Got 5, p-prod="<<parb->PID*par->PID<<std::endl;
	      TLorentzVector hplus = par->P4() + parb->P4();
	      h_mass_tot[isb]->Fill(hplus.M(),weight); // PLot all tb combinations
	      h_mass2_tot[isb]->Fill(hplus.M(),weight); // PLot all tb combinations
	      if (par->PID * parb->PID < 0) {
		h_mass[isb]->Fill(hplus.M(),weight);
		h_mass2[isb]->Fill(hplus.M(),weight);
		TLorentzVector pb = parb->P4();
		double cosThetaStar = GetCosThetaStar(pb,hplus);
		if ( (hplus.M() > lowband1 && hplus.M() < lowband2)
		     || (hplus.M() > highband1 && hplus.M() < highband2)) {
		  h_side_cths[isb]->Fill(cosThetaStar,weight);
		}
		if ( hplus.M() > sigBand1 && hplus.M() < sigBand2) {
		  h_sig_cths[isb]->Fill(cosThetaStar,weight);
		}
		if (cosThetaStar > -0.6) {
		  h_mass_cths[isb]->Fill(hplus.M(),weight); // passing cths cut
		  h_mass2_cths[isb]->Fill(hplus.M(),weight); // passing cths cut
		}
		
	      }
	    }
	  }

	}
      }

    }
  }
  //=============================================================================

  TString opname="sherpplots_";
  opname+=gam;
  opname+=".root";
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(opname);
  if (!f) {
    f = new TFile(opname,"RECREATE");
  }
  for (int i=0;i<3;++i) {
    h_mass[i]->Write();
    h_mass2[i]->Write();
    h_mass_cths[i]->Write();
    h_mass2_cths[i]->Write();
    h_maxbpt[i]->Write();
    h_side_cths[i]->Write();
    h_sig_cths[i]->Write();
  }
  f->Write();
  f->Close();

  
  showsb(h_mass,"mass",gam);
  showsb(h_mass_cths,"mass_cths",gam);
  showsb(h_mass2,"mass2",gam);
  showsb(h_mass2_cths,"mass2_cths",gam);
  //  showsb(h_mass_tot,"mass_tot");
  //  showsb(h_mass2_tot,"mass2_tot");
  showsb(h_minbpt,"minbpt",gam);
  showsb(h_maxbpt,"maxbpt",gam);
  showsb(h_nb,"nb",gam);
  showsb(h_side_cths,"side cths",gam);
  showsb(h_sig_cths,"sig cths",gam);
  /*
    TH1D * h_mass2_b_raw =   (TH1D*)h_mass2[1]->Clone(); // take a copy
    h_mass2[1] = mysmooth(h_mass2[1]); //smooth it
    TH1D * h_mass2_cths_b_raw =   (TH1D*)h_mass2_cths[1]->Clone(); // take a copy
    h_mass2_cths[1] = mysmooth(h_mass2_cths[1]); //smooth it
  */
  subsb(h_mass2,"mass2",gam,mhp);
  subsb(h_mass2_cths,"mass2_cths",gam,mhp);
  subsb(h_sig_cths,"sig_cths",gam,mhp);
  
  return;
}



void showcut(TH1F * plots[2],TString var,TString cut) {
  TString name=var+"_cut_"+cut;
  
  TCanvas * c = new TCanvas("c_"+name,name,600,800);

  TH1F * tot = (TH1F*)plots[0]->Clone("tot"+name);
  tot->Add(plots[1]);
  TH1F * rat = (TH1F*)plots[1]->Clone("rat"+name);
  rat->Divide(tot);

  TH1F * lum = (TH1F*)plots[1]->Clone("lumi"+name);
  double lumi1 = 0.21;
  double lumi2 = 15;

  int n = lum->GetNbinsX();
  for (int i=1;i<=n;++i) {
    double f2 = rat->GetBinContent(i);
    double f1 = 1-f2;
    double lumeff = lumi1*lumi2/(f1*lumi2 + f2*lumi1);
    if (tot->GetBinContent(i) == 0) lumeff = 0;
    lum->SetBinContent(i,lumeff);
  }

  c->Divide(1,3);
  TPad * pad = (TPad*)c->cd(1);
  pad->SetLogy();
  tot->Draw();
  plots[0]->SetLineColor(2);
  plots[0]->Draw("same");
  plots[1]->SetLineColor(3);
  plots[1]->Draw("same");

  TLegend * leg = new TLegend(0.6,0.75,0.9,0.9);
  leg->AddEntry(plots[0],"Fail cut on "+cut,"l");
  leg->AddEntry(plots[1],"Pass cut on "+cut,"l");
  leg->Draw();

 
  pad = (TPad*)c->cd(2);
  rat->SetMaximum(1.);
  rat->GetXaxis()->SetTitle("p_{T} [GeV]");
  rat->GetYaxis()->SetTitle("Fraction passing");
  rat->GetYaxis()->SetTitleSize(0.07);
  rat->GetYaxis()->SetTitleOffset(0.27);
  rat->SetFillColor(5.);
  rat->SetLineColor(1.);
 
  rat->Draw();

  pad = (TPad*)c->cd(3);
 

  lum->GetXaxis()->SetTitle("p_{T} [GeV]");
  lum->GetYaxis()->SetTitle("Effective Luminosity");
  lum->GetYaxis()->SetTitleSize(0.07);
  lum->GetYaxis()->SetTitleOffset(0.17);

  lum->Draw();

  std::cout<<"Fraction in first plot is "<<plots[0]->GetSumOfWeights()/tot->GetSumOfWeights()<<std::endl;
  c->SaveAs(name+".png");
 
}

void showsb(TH1D * plots[3],TString var,double gam) {
  TString name=var+"_";
  name+=gam;

  TCanvas * c = new TCanvas("c_"+name,name,800,600);

  if (var.Contains("pt")) {
    c->SetLogy();
  } else {
    plots[0]->SetMinimum(0.);
  }
  plots[0]->SetMarkerStyle(20);
  plots[0]->Draw();
  plots[1]->SetMarkerColor(2);
  plots[1]->SetMarkerStyle(20);
  plots[1]->SetLineColor(2);
  plots[1]->Draw("same");
  plots[2]->SetMarkerColor(3);
  plots[2]->SetMarkerStyle(20);
  plots[2]->SetLineColor(3);
  plots[2]->Draw("same");

  TLegend * leg = new TLegend(0.6,0.75,0.9,0.9);
  leg->AddEntry(plots[0],"both","l");
  leg->AddEntry(plots[1],"bckgrd","l");
  leg->AddEntry(plots[2],"signal","l");
  leg->Draw();

  c->SaveAs("sherpa_"+name+".png");
}

void subsb(TH1D * plots[3],TString var,double gam,double mhp) {
  TString name=var+"_";
  name+=gam;
  TCanvas * c = new TCanvas("sub_"+name,name,800,600);

  if (var.Contains("pt")) c->SetLogy();

  TH1D * sb_b = (TH1D*)plots[0]->Clone();
  std::cout<<" Adress of plots[0] is "<<plots[0]<<std::endl;
  std::cout<<" Adress of plots[1] is "<<plots[1]<<std::endl;
  std::cout<<" Adress of plots[2] is "<<plots[2]<<std::endl;
  std::cout << plots[0]->GetName() << std::endl;
  std::cout << "Range "<<plots[0]->GetNbinsX()<<" from "<<plots[0]->GetXaxis()->GetXmin()<<" to "<< plots[0]->GetXaxis()->GetXmax()<< std::endl;
  std::cout << plots[1]->GetName() << std::endl;
  std::cout << "Range "<<plots[1]->GetNbinsX()<<" from "<<plots[1]->GetXaxis()->GetXmin()<<" to "<< plots[1]->GetXaxis()->GetXmax()<< std::endl;
  std::cout << plots[2]->GetName() << std::endl;
  std::cout << "Range "<<plots[2]->GetNbinsX()<<" from "<<plots[2]->GetXaxis()->GetXmin()<<" to "<< plots[2]->GetXaxis()->GetXmax()<< std::endl;

  sb_b->Add(plots[1],-1);
  std::cout << "added 0 to 1" << std::endl;
  TH1D * inter = (TH1D*)sb_b->Clone();
  inter->Add(plots[2],-1);

  std::cout << "added plots" << std::endl;

  double themin = plots[2]->GetMinimum();
  if (sb_b->GetMinimum() < themin) themin = sb_b->GetMinimum();
  if (inter->GetMinimum() < themin) themin = inter->GetMinimum();
  double themax = plots[2]->GetMaximum();
  if (sb_b->GetMaximum() > themax) themax = sb_b->GetMaximum();
  if (inter->GetMaximum() > themax) themax = inter->GetMaximum();

  themin = themin - 0.05*(themax - themin);
  themax = themax + 0.05*(themax - themin);
  plots[2]->SetMinimum(themin);
  plots[2]->SetMaximum(themax);

  bool domass = name.Contains("mass");
  TF1 * fbw1, * fbw2;
  if (domass) {
    // TF1 * fbw1 = new TF1("fbw1"+name,"[0]/((x*x-[1]*[1])*(x*x-[1]*[1])/[2]/[2]+[1]*[2])");
    fbw1 = new TF1("fbw1"+name,"[0]*[2]/((x-[1])*(x-[1]) + [2]*[2]/4)");
    fbw1->SetParameter(0,1.E-3/gam/gam);
    fbw1->FixParameter(1,91.18);
    fbw1->FixParameter(2,2.4952);

    // TF1 * fbw2 = new TF1("fbw2"+name,"[0]/((x*x-[1]*[1])*(x*x-[1]*[1])/[2]/[2]+[1]*[2])");
    fbw2 = new TF1("fbw2"+name,"[0]*[2]/((x-[1])*(x-[1]) + [2]*[2]/4)");
    fbw2->SetParameter(0,1.E-3/gam/gam);
    fbw2->FixParameter(1,91.18);
    fbw2->FixParameter(2,2.4952);
  } else {
    fbw1 = new TF1("fcths1"+name,"[0] + [1]*x");
    fbw1->SetParameter(0,2E-3);
    fbw1->SetParameter(1,-1E-4);

    fbw2 = new TF1("fcths2"+name,"[0] + [1]*x");
    fbw2->SetParameter(0,2E-3);
    fbw2->FixParameter(1,-1E-4);
  }

  fbw1->SetLineColor(2);
  fbw1->SetLineStyle(2);

 
  plots[2]->Fit(fbw1,"","");
  sb_b->Fit(fbw2,"","");

  double hs = fbw1->GetParameter(0);
  double dhs = fbw1->GetParError(0);
  double hsb_b = fbw2->GetParameter(0);
  double dhsb_b = fbw2->GetParError(0);

  double rat = hsb_b/hs;
  double drat = dhsb_b/hsb_b*rat; // ignore error on s
 
  plots[2]->SetMarkerColor(3);
  plots[2]->SetMarkerStyle(20);
  plots[2]->SetLineColor(3);
  plots[2]->Draw("");
  sb_b->SetMarkerStyle(20);
  sb_b->SetLineColor(1);
  sb_b->SetMarkerColor(1);
  sb_b->Draw("same");
  inter->SetMarkerStyle(20);
  inter->SetLineColor(4);
  inter->SetMarkerColor(4);
  inter->Draw("same");

 
  TLegend * leg = new TLegend(0.6,0.75,0.9,0.9);
  TString width="Width ";
  width+=gam;
  width+="GeV";
  leg->SetHeader(width,"C");
  leg->AddEntry(plots[2],"s","p");
  char titsb_b[100];
  sprintf(titsb_b,"sb-b, f=%5.3f+-%5.3f",rat,drat);
 
  leg->AddEntry(sb_b,titsb_b,"p");
  leg->AddEntry(inter,"Interference","p");
  leg->Draw();
  double x1 = sb_b->GetXaxis()->GetXmin();
  double x2 = sb_b->GetXaxis()->GetXmax();
  TLine * l = new TLine(x1,0,x2,0);
  l->SetLineStyle(3);
  l->Draw();

  c->SaveAs("sub_"+name+".png");
 
}

double GetCosThetaStar(TLorentzVector & g1,TLorentzVector & gg) {

  TLorentzVector rotvec  = g1;
  rotvec.Boost(-gg.BoostVector());
  //  std::cout<<"About to get cthstar"<<std::endl;
  TVector3 rot3v=rotvec.Vect();
  TVector3 gg3v=gg.Vect();
  Double_t cthstar = rot3v.Dot(gg3v)/rot3v.Mag()/gg3v.Mag();
  //  std::cout<<"Done cthstar"<<std::endl;

  return cthstar;
}

TH1D *  mysmooth(TH1D * hist) {
  // TF1 * fn = new TF1("fn","pol3");
  //q fn->SetLineStyle(2);
  TString tit=hist->GetTitle();

  TCanvas * c = new TCanvas(tit+" smooth",tit+" smooth",800,600);

  // hist->Fit("pol2","","");

  TF1 * fn = new TF1("fnsmooth"+tit,"[0]+[1]*(x-500) + [2]*(x-500)*(x-500)");
  fn->SetLineColor(2);
  fn->SetLineStyle(2);
  fn->SetParameter(0,0.003);
  fn->SetParameter(1,-2E-5);
  fn->SetParameter(2,8E-8);
  hist->Fit(fn,"",""); 

  int nb = hist->GetNbinsX();
  double x1 = hist->GetXaxis()->GetXmin();
  double x2 = hist->GetXaxis()->GetXmax();
  TH1D * result = new TH1D(tit+"_2",tit+"_s",nb,x1,x2);
  // TF1 * fn = hist->GetFunction("pol2");
  double yv[500];
  double dyv[500];
  for (int i=1;i<nb+1;++i) {
    double x = result->GetBinCenter(i);
    yv[i] = fn->Eval(x);
    dyv[i] = 0;
    //   std::cout<<"i="<<i<<" c="<<x<<" yv="<<yv[i]<<std::endl;
  }
  result->SetContent(yv);
  result->SetError(dyv);
  result->Draw("hist same");
  return result;
  
}
