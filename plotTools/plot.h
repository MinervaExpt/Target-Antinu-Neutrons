#ifndef PLOT_H
#define PLOT_H

#include "TClass.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include "TGaxis.h"
#include "TLatex.h"

#include <iostream>
#include <vector>
#include <math.h>
#include <iostream>
#include <algorithm>

#include "PlotUtils/GridCanvas.h"

using namespace std;
using namespace PlotUtils;

//======================================================================

// Return a unique string on every call, so we can name root hists
// without it clobbering them
TString uniq()
{
  static int i=0;
  return TString::Format("uniq%d", i++);
}


//======================================================================

// Convert a histogram to a graph so we can draw it without ROOT
// drawing extra lines at the end and other such nonsense
TGraphErrors* histToGraph(TH1* h, double multiplier=1, bool includeZeros=true, TString opts="")
{
  TGraphErrors* grE=new TGraphErrors;
  grE->SetLineColor(h->GetLineColor());
  grE->SetLineStyle(h->GetLineStyle());
  grE->SetLineWidth(h->GetLineWidth());
  if(opts.Contains("e3")){
      grE->SetFillColor(h->GetFillColor());
      grE->SetFillStyle(h->GetFillStyle());
  }
  for(int i=0; i<h->GetNbinsX(); ++i){
    int bin=i+1;
    double binVal=h->GetBinContent(bin);
    if(includeZeros || binVal!=0){
      int binNo=grE->GetN();
      grE->SetPoint(binNo, h->GetBinCenter(bin), multiplier*binVal);
      grE->SetPointError(binNo, 0, multiplier*h->GetBinError(bin));
    }
  }
  return grE;
}

//======================================================================

TH1* getPadAxisHist(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    // cout << obj->GetName() << endl;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      // cout << "getPadAxisHist returning hist with name " << obj->GetName() << endl;
      return (TH1*)obj;
    }
  }
  return NULL;
}

//======================================================================
std::vector<TH1*> getPadHists(TPad* pad)
{
  std::vector<TH1*> ret;
  std::vector<string> names;
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    //    cout << obj->GetName() << endl;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      string name = obj->GetName();
      //      cout << "getPadAxisHist returning hist with name " << obj->GetName() <<"t"<< obj->GetTitle() << endl;
      bool isFound = false;
      for(unsigned int i=0;i<names.size();i++){
	if(names[i]==name){
	  isFound = true;
	  break;
	}
      }
      if(!isFound){
	ret.push_back((TH1*)obj);
	names.push_back(name);
      }

    }
  }
  return ret;
}

//======================================================================
TLegend* getPadLegend(TPad* pad)
{
  cout <<"Looking at pad " <<  pad << endl;
  TLegend *ret = NULL;
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    cout << obj->GetName() << endl;
    cout << obj->ClassName() << endl;
    if ( obj->InheritsFrom("TLegend") ) {
      string name = obj->GetName();
      cout << "I found one" << endl;
      ret= (TLegend*)obj;
    }
  }
  cout << "I am done " << ret << endl;
  return ret;
}

//======================================================================
double getPadMax(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  Double_t runningMax=-9e99;//Hparam.ymax;
  while (( obj=next() )) {
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      TH1* curHist=(TH1*)obj;
      const double thisMax=curHist->GetBinContent(curHist->GetMaximumBin());
      if (thisMax > runningMax) {
        runningMax=thisMax;
      }
    }
  }
  return runningMax;
}

//======================================================================
void reMax(TPad* pad, double ymin=0, double headroom=1.1)
{
  TH1* firstHist=getPadAxisHist(pad);
  if(!firstHist){
    std::cerr << "reMax: No histogram in pad " << pad->GetName() << ". Can't reMax" << std::endl;
    return;
  }
  firstHist->GetYaxis()->SetRangeUser(ymin, headroom*getPadMax(pad));
  pad->Update();
}

//======================================================================
void drawBinRange(TH2* h, int axis, int bin, const char* varName, const char* numFormatStr=".2f", int mode=1, int gridx = 0, int gridy = 0)
{

  double varmin=axis==1 ? h->GetXaxis()->GetBinLowEdge(bin) : h->GetYaxis()->GetBinLowEdge(bin);
  double varmax=axis==1 ? h->GetXaxis()->GetBinUpEdge(bin) :  h->GetYaxis()->GetBinUpEdge(bin);
  
  TString formatStr(TString::Format("%%%s < %%s < %%%s", numFormatStr, numFormatStr));

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));

  if(mode==1){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else if(mode==0){
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }
  else if(mode==2){
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  gPad->GetBottomMargin()+0.05,
                  text);
    la->SetTextAlign(31); // bottom right
  }

  la->SetNDC();
  la->SetTextFont(42);
  //Okay math time. The grid will be 150 * ncolumns wide and 150*nrows tall
  //SetTextSize is based off the aspect ratio and then the character height is defined as (size)*pad_width (tall)or (size)*pad_height (wide)
  // This is a problem because the pads for all "frames" are the full canvas with modified margins. That means the size of the text stays constant for any number of "frames". We want smaller text for more frames.
  // A good size for wide plots has been 24. 800 pixel wide * 0.03
  // So, 150*ncolumn = width we want 24 -> 24/(150*ncolumn) or 24/(150*nrow)
  if(gridx==0||gridy==1) la->SetTextSize(0.03);
  else if(gridx<gridy) la->SetTextSize(28/(150.0*gridy));//Tall
  else if(gridx==gridy) la->SetTextSize(10/(150.0*gridx));//Wide
  else la->SetTextSize(18/(150.0*gridx));//Wide
  la->Draw();
}

//======================================================================
void drawVarBinRange(TH2* h, int axis, int minbin, int maxbin, const char* varName, const char* numFormatStr=".2f", bool left=false)
{
  double varmin=axis==1 ? h->GetXaxis()->GetBinLowEdge(minbin) : h->GetYaxis()->GetBinLowEdge(minbin);
  double varmax=axis==1 ? h->GetXaxis()->GetBinUpEdge(maxbin) :  h->GetYaxis()->GetBinUpEdge(maxbin);
  
  TString formatStr(TString::Format("%%%s < %%s < %%%s", numFormatStr, numFormatStr));

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));
  if(left){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  la->SetTextSize(0.03);
  la->Draw();
}

//======================================================================
TH1* rebinpz(TH1* h)
{
  const int nBins=14;
  // Less aggressive version
  // double newBins[nBins+1]={1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 12, 14};
  //1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20,40,60
  double newBins[nBins+1]={1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8,9,10,13,16};

  TH1* ret=new TH1D(uniq(), TString::Format("%s;%s;%s", h->GetTitle(), h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle()),
                    nBins, newBins);

  ret->SetLineColor(h->GetLineColor());
  ret->SetLineStyle(h->GetLineStyle());
  ret->SetLineWidth(h->GetLineWidth());
  ret->SetMarkerStyle(h->GetMarkerStyle());
  ret->SetMarkerSize(h->GetMarkerSize());
  
  // Do under/overflow too
  for(int i=0; i<nBins+2; ++i){
    ret->SetBinContent(i, h->GetBinContent(i));
    ret->SetBinError(i, h->GetBinError(i));
  }
  return ret;
}

TH1* rebinpt(TH1* h)
{
  const int nBins=14;
  // Less aggressive version
  // double newBins[nBins+1]={1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 12, 14};
  //1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20,40,60
  double newBins[nBins+1]={0,0.07,0.15,0.25,0.33,0.4,0.47,0.55,0.70,0.85,1.0,1.25,1.5,2,2.5};

  TH1* ret=new TH1D(uniq(), TString::Format("%s;%s;%s", h->GetTitle(), h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle()),
                    nBins, newBins);

  ret->SetLineColor(h->GetLineColor());
  ret->SetLineStyle(h->GetLineStyle());
  ret->SetLineWidth(h->GetLineWidth());
  ret->SetMarkerStyle(h->GetMarkerStyle());
  ret->SetMarkerSize(h->GetMarkerSize());
  
  // Do under/overflow too
  for(int i=0; i<nBins+2; ++i){
    ret->SetBinContent(i, h->GetBinContent(i));
    ret->SetBinError(i, h->GetBinError(i));
  }
  return ret;
}


//======================================================================
GridCanvas* plotpT1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();


  int grid_x = sqrt(nbins_pz)+1;
  int grid_y = nbins_pz/(grid_x-1);
  //  grid_x = 4;
  //  grid_y = 4;
  cout << "Plotting plotpT1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, "var", ".1f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Muon transverse momentum (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotYAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
			double* multipliers=NULL)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();


  int grid_x = sqrt(nbins_pz)+1;
  int grid_y = nbins_pz/(grid_x-1);
  if(grid_x*grid_y-nbins_pz==grid_x) grid_y-=1;
  cout << nbins_pz - grid_x*grid_y << endl;
  cout << "Plotting plotYAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  //  grid_x = 4;
  //  grid_y = 4;
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 1600, 1000);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle(xaxistitle.c_str());
  //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotYAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int ncolumns, int nrows, int xwidth, int ywidth,
			double* multipliers=NULL)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  //  int nbins_pt = histAndOpts[0].first->GetNbinsY();


  GridCanvas* gc=new GridCanvas(uniq(), ncolumns, nrows, xwidth, ywidth);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }
    cout << pad << "\t" << gPad << "\tI have these" << endl;
    //    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false, ncolumns,nrows);
    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex();
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      /*
      double textsize = 0;
      if(nrows==1) textsize=0.035;
      else if(ncolumns<nrows) textsize=28/(150.0*nrows);
      else if(ncolumns==nrows) textsize=10/(150.0*ncolumns);
      else textsize = 18/(150.0*ncolumns);
      la2->SetTextSize(textsize);
      la2->DrawLatexNDC(1-gPad->GetRightMargin()-0.01, 1-gPad->GetTopMargin()-30/(150.0*nrows),TString::Format("#times %.1f", multipliers[i]));
      */
      la2->SetTextSize(0.035);
      la2->DrawLatexNDC(1-gPad->GetRightMargin()-0.01, 1-gPad->GetTopMargin()-30/(150.0*nrows),TString::Format("#times %.1f", multipliers[i]));
      la2->Draw();
    }
  }

  gc->SetXTitle(xaxistitle.c_str());
  //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotYAxis1DRebinPt(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int ncolumns, int nrows, int xwidth, int ywidth,
			double* multipliers=NULL)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  //  int nbins_pt = histAndOpts[0].first->GetNbinsY();


  GridCanvas* gc=new GridCanvas(uniq(), ncolumns, nrows, xwidth, ywidth);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* projtmp=h2d->ProjectionY(uniq(), i+1, i+1);
      TH1* proj=rebinpt(projtmp);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }
      pad->RedrawAxis();
    }
    
    cout << pad << "\t" << gPad << "\tI have these" << endl;
    //    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false, ncolumns,nrows);
    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex();
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      /*
      double textsize = 0;
      if(nrows==1) textsize=0.035;
      else if(ncolumns<nrows) textsize=28/(150.0*nrows);
      else if(ncolumns==nrows) textsize=10/(150.0*ncolumns);
      else textsize = 18/(150.0*ncolumns);
      la2->SetTextSize(textsize);
      la2->DrawLatexNDC(1-gPad->GetRightMargin()-0.01, 1-gPad->GetTopMargin()-30/(150.0*nrows),TString::Format("#times %.1f", multipliers[i]));
      */
      la2->SetTextSize(0.035);
      la2->DrawLatexNDC(1-gPad->GetRightMargin()-0.01, 1-gPad->GetTopMargin()-30/(150.0*nrows),TString::Format("#times %.1f", multipliers[i]));
      la2->Draw();
    }
  }


  const int nLabels=6;
  const double positions[nLabels]={0.4,1,1.5,2,2.5};
  const char* valueStrings[nLabels]={"0.4","1.0","1.5","2.5","4.5"};
  
  gc->SetManualXLabels(nLabels, positions, valueStrings,0.03);


  gc->SetXTitle(xaxistitle.c_str());
  //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotpT1DAntiNu(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 4, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<11; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, "p_{||}/GeV", ".1f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Muon transverse momentum (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/nucleon)");
  gc->SetTitleSize(22.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotpz1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL, bool is3D=false)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1)-1;
  //  grid_x = 5;
  //  grid_y = 3;
  cout << "Plotting plotpZ1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* projtmp=h2d->ProjectionX(uniq(), i+1, i+1);
      TH1* proj=is3D?(TH1*)projtmp->Clone("clone"):rebinpz(projtmp);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  const int nLabels=3; 
  const double positions[nLabels]={7, 13, 16}; 
    const char* valueStrings[nLabels]={"5", "10", "20"}; 

  if(!is3D){
    gc->SetManualXLabels(nLabels, positions, valueStrings);
    gc->SetXTitle("Muon longitudinal momentum (GeV)");
  }
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotXAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  //  grid_x = 5;
  //  grid_y = 3;
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->SetFillColor(h2d->GetFillColor());
      proj->SetFillStyle(h2d->GetFillStyle());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotXAxis1D_MoveLabels(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
			double* multipliers=NULL, bool moveLabels=false)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  //  grid_x = 5;
  //  grid_y = 3;
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->SetFillColor(h2d->GetFillColor());
      proj->SetFillStyle(h2d->GetFillStyle());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    if(i>=8)drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f",2);
    else drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f",0);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotXAxis1D_IgnoreYBins(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int lowYbin, int upYbin,
				    double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = upYbin-lowYbin+1;

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x)+1;
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  int plot_counter=1;
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    if(i+1 < lowYbin || i>upYbin) continue;
    TPad* pad=(TPad*)gc->cd(plot_counter);
    cout <<i << endl;
    plot_counter+=1;
    cout << "\t" << i << endl;
    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->SetFillColor(h2d->GetFillColor());
      proj->SetFillStyle(h2d->GetFillStyle());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }
  
  for(int i=0;i<grid_x*grid_y;i++){
    TPad *p = (TPad*)gc->cd(i+1);
    TList *l = p->GetListOfPrimitives();
    //    if(l->First()==0) p->SetFillStyle(0);
    p->SetFillStyle(0);
    p->Modified();
  }

  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotXAxis1D_IgnoreYBins_ReduceXRange(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int lowYbin, int upYbin, double lowX, double highX,
				    double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = upYbin-lowYbin+1;

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x)+1;
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  if(lowYbin==upYbin) grid_x=1;
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  vector<int> mybins;
  for(int i=0; i<nbins_pt; ++i){

    bool a = i+1 < lowYbin;
    bool b = i>upYbin;
    cout << "My i check " << i  << "\t" <<a << "\t" << b<< endl;
    //    if(i+1 < lowYbin || i>upYbin) continue;
    if(i+1 >= lowYbin && i<=upYbin) continue;
    cout << "\t" << i << endl;
    mybins.push_back(i+lowYbin);
  }

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  //  for(int i=0; i<nbins_pt; ++i){
  for(int index=0;index<mybins.size();index++){
    cout << index << "\t" << mybins[index] << endl;
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(index+1);
    int i=mybins[index];
    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);
      proj->GetXaxis()->SetRangeUser(lowX,highX);
      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->SetFillColor(h2d->GetFillColor());
      proj->SetFillStyle(h2d->GetFillStyle());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}



//======================================================================
GridCanvas* plotXAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int grid_x, int grid_y, vector<int> panels,
			double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  if((int)panels.size()>grid_x*grid_y){
    cout << "You've picked an incompatible grid layout and number of panels. Exiting" << endl;
    cout << "You have a vector of " << panels.size() << " and chose a gridx*gridy of " << grid_x*grid_y << endl;
    exit(1);
  }
      
  int plot_panel_count = 0;
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    bool plot = true;
    if(panels.size()!=0){
      bool test = false;
      for(unsigned int j=0;j<panels.size();j++){
	if(panels[j]==i+1) test=true;
      }
      plot=test;
    }
    
    if(!plot) continue;
    plot_panel_count+=1;
    TPad* pad=(TPad*)gc->cd(plot_panel_count);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->SetFillColor(h2d->GetFillColor());
      proj->SetFillStyle(h2d->GetFillStyle());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}



//======================================================================
GridCanvas* plotXAxis1D_ReducedXRange(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, double xminval, double xmaxval,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800,500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);
      proj->SetFillStyle(h2d->GetFillStyle());

      proj->GetXaxis()->SetRangeUser(xminval,xmaxval);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotXAxis1D_ReducedXRange(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int grid_x, int grid_y, double xminval, double xmaxval, double xpix, double ypix,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  //  int grid_x = sqrt(nbins_pt)+1;
  //  int grid_y = nbins_pt/(grid_x-1);
  //  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, xpix,ypix);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);
      proj->SetFillStyle(h2d->GetFillStyle());
      //cout << "Setting (Min,max) to: " << xminval << "\t" << xmaxval << endl;
      proj->GetXaxis()->SetRangeUser(xminval,xmaxval);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
    GridCanvas* plotYAxis1D_ReducedXRange(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int grid_x, int grid_y, double xminval, double xmaxval, double xpix, double ypix,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  //  int grid_x = sqrt(nbins_pt)+1;
  //  int grid_y = nbins_pt/(grid_x-1);
  //  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotYAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, xpix,ypix);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);
      proj->SetFillStyle(h2d->GetFillStyle());

      proj->GetXaxis()->SetRangeUser(xminval,xmaxval);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotXAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int ncolumns, int nrows, int xwidth, int ywidth,
                     double* multipliers=NULL)
{

  
  //  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  GridCanvas* gc=new GridCanvas(uniq(), ncolumns, nrows, xwidth, ywidth);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);
    pad->GetFillColor();
    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);
      proj->SetFillColor(h2d->GetFillColor());
      proj->SetFillStyle(h2d->GetFillStyle());


      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    //    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false, ncolumns, nrows);
    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false );

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){

      TLatex* la2=new TLatex();
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      double textsize = 0;
      if(nrows==1) textsize=0.035;
      else if(ncolumns<nrows) textsize=28/(150.0*nrows);
      else if(ncolumns==nrows) textsize=10/(150.0*ncolumns);
      else textsize = 18/(150.0*ncolumns);
      la2->SetTextSize(textsize);
      la2->DrawLatexNDC(1-gPad->GetRightMargin()-0.01, 1-gPad->GetTopMargin()-30/(150.0*nrows),TString::Format("#times %.1f", multipliers[i]));

    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotXAxis1DRebinPz(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, int ncolumns, int nrows, int xwidth, int ywidth,
                     double* multipliers=NULL)
{

  
  //  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  GridCanvas* gc=new GridCanvas(uniq(), ncolumns, nrows, xwidth, ywidth);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();


      TH1* projtmp=h2d->ProjectionX(uniq(), i+1, i+1);
      TH1* proj=rebinpz(projtmp);


      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }
      pad->RedrawAxis();
    }


    

    //    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false,ncolumns,nrows);
    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  const int nLabels=5;
  const double positions[nLabels]={4,10,16,20,24};
  const char* valueStrings[nLabels]={"4","10","20", "40", "60"};
  
  gc->SetManualXLabels(nLabels, positions, valueStrings,0.03);
  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotXAxis1D_Nue(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 2, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<5; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);
      proj->GetXaxis()->SetRangeUser(0,100);
      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotpz1DAntiNu(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 3, 2, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<6; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);
      //      TH1* proj=rebinpz(projtmp);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Muon longitudinal momentum (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/nucleon)");
  gc->SetTitleSize(22.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotvtx1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
		      double* multipliers=NULL, TH2* denominator=NULL, bool track2=false)
{
  cout << "Doing ratios" << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=track2? new GridCanvas(uniq(), 3, 2, 4000, 2000): new GridCanvas(uniq(), 3, 3, 4000, 2000);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();
  
  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them


  const int nbins = track2 ? 4:8;  
  int minb2[4]={1,3,6,9};
  int maxb2[4]={2,5,8,13};

  int minb[8]={1,3,4,5,6,7,8,9};
  int maxb[8]={2,3,4,5,6,7,8,13};

  
  for(int i=0; i<nbins; ++i){
    // Change to the appropriate pad in the canvas
    
    int val = i>1?i+2:i+1;
    if(!track2)val=i+1;
    TPad* pad=(TPad*)gc->cd(val);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();
      TH1* projtmp = track2? h2d->ProjectionX(uniq(), minb2[i],maxb2[i]):h2d->ProjectionX(uniq(), minb[i],maxb[i]);
      TH1* demproj = track2? denominator->ProjectionX(uniq(), minb2[i],maxb2[i]):denominator->ProjectionX(uniq(), minb[i],maxb[i]);
      projtmp->Divide(demproj);

      TH1* proj=(TH1*)projtmp->Clone();
      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
	if(opt.Contains("histe")){
	  proj->SetFillColor(kWhite);
	  proj->SetFillStyle(0);
	  proj->DrawClone("histc");
	  proj->SetLineWidth(2);
	  proj->SetFillColor(kRed);
	  proj->SetFillStyle(3003);
	  proj->Draw("E5SAME");
	}
	else proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
	  
      }

    }

    if(track2) drawVarBinRange(histAndOpts[0].first, 2, minb2[i],maxb2[i], "p_{T}/GeV", ".2f", false);
    else  drawVarBinRange(histAndOpts[0].first, 2, minb[i],maxb[i], "p_{T}/GeV", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Vertex Energy");
  gc->SetYTitle("Events");
  gc->SetTitleSize(100.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotvtxrate1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
			  double* multipliers=NULL,bool track2=false, bool binbybin =false)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=track2? new GridCanvas(uniq(), 3, 2, 4000, 2000): new GridCanvas(uniq(), 3, 3, 4000, 2000);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();
  
  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them


  const int nbins = track2 ? 4:8;  
  int minb2[4]={1,3,6,9};
  int maxb2[4]={2,5,8,13};

  int minb[8]={1,3,4,5,6,7,8,9};
  int maxb[8]={2,3,4,5,6,7,8,13};


  
  for(int i=0; i<nbins; ++i){
    // Change to the appropriate pad in the canvas
    
    int val = i>1?i+2:i+1;
    if(!track2)val=i+1;
    TPad* pad=(TPad*)gc->cd(val);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      //      h2d->SaveAs(Form("%s-%s-%d.root",uniq().Data(),h2d->GetName(),track2));
      opt.ToLower();
      TH1* proj= NULL;
      if(!binbybin) proj = track2? h2d->ProjectionX(uniq(), minb2[i],maxb2[i]):h2d->ProjectionX(uniq(), minb[i],maxb[i]);
      else proj = h2d->ProjectionX(uniq(),i+1,i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);
      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
	if(opt.Contains("histe")){
	  proj->SetFillColor(kWhite);
	  proj->SetFillStyle(0);
	  proj->DrawClone("histc");
	  proj->SetLineWidth(2);
	  proj->SetFillColor(kRed);
	  proj->SetFillStyle(3003);
	  proj->Draw("E5SAME");
	}
	else proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
	  
      }

    }
    if(!binbybin){
      if(track2) drawVarBinRange(histAndOpts[0].first, 2, minb2[i],maxb2[i], "p_{T}/GeV", ".2f", false);
      else  drawVarBinRange(histAndOpts[0].first, 2, minb[i],maxb[i], "p_{T}/GeV", ".2f", false);
    }
    else{
      drawBinRange(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", false);
    }

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Vertex Energy");
  gc->SetYTitle("Events");
  gc->SetTitleSize(100.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

vector<double> GetScales(std::vector<std::pair<TH2*, const char*> >histopts, bool pxProj, double plotMax, double fracMax, bool limitRange = false, double rmin=0, double rmax=100){
  vector<double> tmpvect;
  int nbins = histopts[0].first->GetNbinsX()+2;
  if(pxProj) nbins = histopts[0].first->GetNbinsY()+2;
  for(int i=1;i<nbins-1;i++){
    double maxval = 0;
    for(uint j=0;j<histopts.size();j++){
      TH1D *tmp = pxProj? histopts[j].first->ProjectionX("tmp",i,i): histopts[j].first->ProjectionY("tmp",i,i);
      if(limitRange) tmp->GetXaxis()->SetRangeUser(rmin,rmax);
      int maxbin = tmp->GetMaximumBin();
      double content = tmp->GetBinContent(maxbin);
      if(content>maxval) maxval=content;
    }
    //we want about fracMax of the plotMax. 
    double scale = (plotMax*fracMax)/maxval;
    if(scale>1){
      int tmpscale = floor(scale*10);
      scale = tmpscale/10.0;
    }
    else{
      int tmpscale = ceil(scale*10);
      scale = tmpscale/10.0;
    }
    cout << scale << endl;
    tmpvect.push_back(scale);
  }
    return tmpvect;
}

vector<double> GetScalesRatio(std::vector<std::pair<TH2*, const char*> >histopts, bool pxProj, double plotMax){
  vector<double> tmpvect;
  int nbins = histopts[0].first->GetNbinsX()+2;
  if(pxProj) nbins = histopts[0].first->GetNbinsY()+2;
  for(int i=1;i<nbins-1;i++){
    double maxval = 0;
    for(uint j=0;j<histopts.size();j++){
      TH1D *tmp = pxProj? histopts[j].first->ProjectionX("tmp",i,i): histopts[j].first->ProjectionY("tmp",i,i);
      int maxbin = tmp->GetMaximumBin();
      double content = tmp->GetBinContent(maxbin);
      if(content>maxval) maxval=content;
    }
    //we want about fracMax of the plotMax. 
    //This is a ratio -> we only want easy numbers
    double scale = (plotMax*0.75)/maxval;
    if(scale>1) scale = 1;
    if(scale<1 && scale >0.5) scale = 0.5;
    if(scale<0.5) scale = 0.25;
    tmpvect.push_back(scale);
  }
    return tmpvect;
}


///////////////////Methods for nongrid plotting////////////////////
vector<vector<double>> GetKinValue(int binsx,int binsy,TH1 *rebinnedhisto,TH2 *templatehist) {
  vector<vector <double>> tmpv;
  for (int i = 1; i < binsy-1; i++) {
    vector<double>binsptbin;
    for (int j = 0; j < binsx-1; j++) {
      int binid      = rebinnedhisto->FindBin(i*binsx+j+1);
      int mapped_low = rebinnedhisto->GetBinLowEdge(binid);
      int xindex     = mapped_low%binsx;
      double kin_val = templatehist->GetXaxis()->GetBinLowEdge(int(xindex));
      binsptbin.push_back(kin_val);    
      cout << i << " " << j << "  "<< binid << " " << mapped_low << " " << xindex  << " " << kin_val << endl; 
    } 
    //it = unique(binsptbin.begin(), binsptbin.end());
    //binsptbin.resize(distance(binsptbin.begin(), it));
    
    std::sort(binsptbin.begin(), binsptbin.end()); 
    auto last = std::unique(binsptbin.begin(), binsptbin.end());
    binsptbin.erase(last, binsptbin.end());

    cout << " Bin " << i << " : ";
    for(auto it:binsptbin) {
      cout <<it << " ";
    }
    cout << endl;
    
    tmpv.push_back(binsptbin);
  }
  return tmpv;
}

TH1* makeTH1D(int pt_bin,double pt_width, vector<double> boundaries,TH1* rebinnedhisto,int startindex,const char *prefix,double scale=1.0,bool doRatio=false) {
  
  int n = boundaries.size();
  


  double arr[n];
  cout << " Ptbin " << pt_bin << " " << startindex << " ";

  for (unsigned int i = 0; i < boundaries.size(); i++) {
    arr[i] = boundaries[i];
    cout  << arr[i] << " ";

  }
  
  TH1* retTH1D = new TH1D(Form("%s_pt_bin_%d",prefix,pt_bin),Form("%s_pt_bin_%d", prefix,pt_bin),n-1,arr);
  
  for (unsigned int i = 0; i < boundaries.size(); i++) {
    
    double bv      = doRatio?rebinnedhisto->GetBinContent(startindex+i):rebinnedhisto->GetBinContent(startindex+i)/(pt_width);
    double bv_err  = doRatio?rebinnedhisto->GetBinError(startindex+i):rebinnedhisto->GetBinError(startindex+i)/(pt_width);

    retTH1D->SetBinContent(i+1,bv);
    retTH1D->SetBinError(i+1,bv_err);
  }
  
  //scale by width
  if(!doRatio)retTH1D->Scale(scale,"width");//Now the histograms are width scaled by the area of the bin.
  return retTH1D;
}

//_________________________________
vector <TH1*> MakeTH1Ds(vector<vector<double>> kinVs, TH1 *rebinnedhisto,TH2 *templatehist, const char *prefix,double scale=1.0,bool doRatio=false) {
  vector <TH1*> myTH1Ds;
  int counter = templatehist->GetNbinsX()+4;
  int i = 0;
  cout <<"--------------------------------------------"<<endl;
  for (auto& it : kinVs) {
    
    TH1* tmp_th1D = makeTH1D(i,templatehist->GetYaxis()->GetBinWidth(i+1),it,rebinnedhisto,counter,prefix,scale,doRatio);
    counter+=it.size()+1;
    myTH1Ds.push_back(tmp_th1D);

    cout << " Bin " << i <<" : ";
    for (auto& itt : it) {
      cout << itt << " ";
    }
    cout <<endl;
    i++;
  }
  
  return myTH1Ds;
}

//___________________________________________________________________________
vector <TH1*> PrepareHist(TH2* hisy_temp, TH1* dataMnv, const char *prefix,double scale=1.0,bool doRatio=false){
  int n_binsx_original = hisy_temp->GetNbinsX()+2;
  int n_binsy_original = hisy_temp->GetNbinsY()+2;
  cout << "==================== " << n_binsx_original << "  " << n_binsy_original <<endl;
  vector<vector<double>> kin_value = GetKinValue(n_binsx_original, n_binsy_original, dataMnv, hisy_temp);
  vector <TH1*> mk1Ds = MakeTH1Ds(kin_value, dataMnv, hisy_temp, prefix,scale,doRatio);


  for (unsigned int i = 0; i < mk1Ds.size(); i++){
    mk1Ds[i]->SetLineColor(dataMnv->GetLineColor());
    mk1Ds[i]->SetLineStyle(dataMnv->GetLineStyle());
    mk1Ds[i]->SetLineWidth(dataMnv->GetLineWidth());
    mk1Ds[i]->SetMarkerStyle(dataMnv->GetMarkerStyle());
    mk1Ds[i]->SetMarkerSize(dataMnv->GetMarkerSize());
    mk1Ds[i]->GetXaxis()->SetNdivisions(504);
    mk1Ds[i]->GetYaxis()->SetNdivisions(504);
  }
  return mk1Ds;

}

std::vector<double> GetScales(std::vector <std::pair<std::vector<TH1*> ,const char*> > histopts,double yaxismax,double fraction){
  std::vector<double> tmpvect;
  cout  << " <<<<<<< Calculating the Scale Factor >>>>>>>>>> " << endl;
  int nh = histopts.size();
  int nb = histopts[0].first.size();
  
  double mxv[nb][nh];
  
  cout << nh  << "  " << nb << endl;
  int i = 0;
  for (auto& it : histopts) {
    int bin = 0;
    for (auto& itt : it.first) {
      int maxbin = itt->GetMaximumBin();
      mxv[bin][i] = itt->GetBinContent(maxbin);
      cout << itt->GetName() << " " << i << "  " << bin << " " << mxv[bin][i] << " " << maxbin  << endl;
      bin++;
    }
    i++;
  }
  

  for(int ii = 0; ii < nb; ii++){
    double maxval = 0;
    for(int jj = 0; jj < nh; jj++){
      cout  << ii << " " << " " << jj << " " << mxv[ii][jj] << endl;
      double content =  mxv[ii][jj];
      if (content > maxval) maxval = content;
    }
    cout << maxval << endl;
    double scale = yaxismax/maxval*fraction;
    if(scale>1){
       int tmpscale = floor(scale*10);
       scale = (tmpscale/10.0)-0.5;
       if (scale > 10) scale = scale - 4;
    }
    else{
      int tmpscale = ceil(scale*1000);
      scale = tmpscale/1000.0;
      }
    //scale = scale;
    cout << "Bin " << ii << " Scale Factor " << scale << endl;
    tmpvect.push_back(scale);
  }


  cout << " <<<<<<<<<<< ------------- end scale value calculation  --------- >>>>>>>>> " << endl;
  return tmpvect;
}


//======================================================================
GridCanvas* plotNonGrid1D(std::vector<std::pair<std::vector<TH1*>, const char*> > histAndOpts, string xaxistitle, string celltitle, TH2* mytemplate,
			  double* multipliers=NULL)
{

  
  int nbins_pt = histAndOpts[0].first.size();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  //  grid_x = 5;
  //  grid_y = 3;
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  //  cout << "Plotting plotNonGrid1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH1* proj=histAndOpts[j].first[i];
      TString opt(histAndOpts[j].second[i]);
      opt.ToLower();

      
      proj->GetXaxis()->SetNdivisions(504);
      proj->GetYaxis()->SetNdivisions(504);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros,opt);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(mytemplate, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}



#endif
