#ifndef VARIABLE2D_H
#define VARIABLE2D_H

#include "util/SafeROOTName.h"
#include "PlotUtils/Variable2DBase.h"
#include "util/Categorized.h"

class Variable2D: public PlotUtils::Variable2DBase<CVUniverse>
{
  private:
    typedef PlotUtils::Hist2DWrapper<CVUniverse> Hist;
    bool fAnaVar;
    bool fFillVar;
    TString fDirName;
  public:
    template <class ...ARGS>
    Variable2D(bool isAnalysisVar, ARGS... args): PlotUtils::Variable2DBase<CVUniverse>(args...), fAnaVar(isAnalysisVar), fFillVar(true), fDirName("TwoD")
    {
    }

    template <class ...ARGS>
    Variable2D(bool isAnalysisVar, bool fillVar, ARGS... args): PlotUtils::Variable2DBase<CVUniverse>(args...), fAnaVar(isAnalysisVar), fFillVar(fillVar), fDirName("TwoD")
    {
    }

    template <class ...ARGS>
    Variable2D(TString name, bool isAnalysisVar, ARGS... args): PlotUtils::Variable2DBase<CVUniverse>(args...), fAnaVar(isAnalysisVar), fFillVar(true), fDirName(name+"/TwoD")
    {
    }

    template <class ...ARGS>
    Variable2D(TString name, bool isAnalysisVar, bool fillVar, ARGS... args): PlotUtils::Variable2DBase<CVUniverse>(args...), fAnaVar(isAnalysisVar), fFillVar(fillVar), fDirName(name+"/TwoD")
    {
    }

    //TODO: It's really silly to have to make 2 sets of error bands just because they point to different trees.
    //      I'd rather the physics of the error bands remain the same and just change which tree they point to.
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
                           std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands)
    {

      std::map<int, std::string> BKGLabels = {{1, "1chargePi"},
					      {2, "1neutPi"},
					      {3, "NPi"},
					      {62, "USPlastic"},
					      {63, "DSPlastic"},
					      {44, "Wrong Nucleus"}};  
      
      /*
      std::map<int, std::string> BKGLabels = {{1, "1chargePi"},
                                              {2, "1neutPi"},
                                              {3, "NPi"}};
      */

      //Various breakdowns for signal and background. TODO: Get the mapping correct                                                                                                                                
      std::map<int, std::string> IntTypeLabels = {{1, "QE"},
                                                  {2, "RES"},
                                                  {3, "DIS"},
                                                  {8, "2p2h"}};

      std::map<int, std::string> TargetTypeLabels = {{1, "Plastic"},
						     {2, "USPlastic"},
						     {3, "DSPlastic"},
                                                     {6, "C"},
                                                     {8, "Water"},
                                                     {26, "Fe"},
                                                     {82, "Pb"}};


      std::map<int, std::string> LeadBlobTypeLabels = {{2, "neut"},
                                                       {3, "prot"},
                                                       {4, "pi0"},
                                                       {5, "pip"},
                                                       {6, "pim"},
                                                       {9, "mu"},
                                                       {1, "None"}};
      
      m_backgroundHists = new util::Categorized<Hist, int>(("TwoD_"+GetName() + "_background").c_str(),
							   ("TwoD_" + GetName()).c_str(), BKGLabels,
							   GetBinVecX(), GetBinVecY(), mc_error_bands);

      m_SigIntTypeHists = new util::Categorized<Hist, int>(("TwoD_"+GetName() + "_sig_IntType").c_str(),
                                                           ("TwoD_"+GetName()).c_str(), IntTypeLabels,
                                                           GetBinVecX(), GetBinVecY(), mc_error_bands);

      m_SigTargetTypeHists = new util::Categorized<Hist, int>(("TwoD_"+GetName() + "_sig_TargetType").c_str(),
							      ("TwoD_"+GetName()).c_str(), TargetTypeLabels,
							      GetBinVecX(), GetBinVecY(), mc_error_bands);
      

      m_SigLeadBlobTypeHists = new util::Categorized<Hist, int>(("TwoD_"+GetName() + "_sig_LeadBlobType").c_str(),
								("TwoD_"+GetName()).c_str(), LeadBlobTypeLabels,
								GetBinVecX(), GetBinVecY(), mc_error_bands);

      m_BkgIntTypeHists = new util::Categorized<Hist, int>(("TwoD_"+GetName() + "_bkg_IntType").c_str(),
                                                           ("TwoD_"+GetName()).c_str(), IntTypeLabels,
                                                           GetBinVecX(), GetBinVecY(), mc_error_bands);

      m_BkgTargetTypeHists = new util::Categorized<Hist, int>(("TwoD_"+GetName() + "_bkg_TargetType").c_str(),
							      ("TwoD_"+GetName()).c_str(), TargetTypeLabels,
							      GetBinVecX(), GetBinVecY(), mc_error_bands);

      m_BkgLeadBlobTypeHists = new util::Categorized<Hist, int>(("TwoD_"+GetName() + "_bkg_LeadBlobType").c_str(),
								("TwoD_"+GetName()).c_str(), LeadBlobTypeLabels,
								GetBinVecX(), GetBinVecY(), mc_error_bands);

      efficiencyNumerator = new Hist(("TwoD_" + GetName() + "_efficiency_numerator").c_str(), ("TwoD_"+GetName()).c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
      efficiencyDenominator = new Hist(("TwoD_" + GetName() + "_efficiency_denominator").c_str(), ("TwoD_" + GetName()).c_str(), GetBinVecX(), GetBinVecY(), truth_error_bands);
      selectedSignalReco = new Hist(("TwoD_" + GetName() + "_selected_signal_reco").c_str(),  ("TwoD_"+GetName()).c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
      selectedMCReco = new Hist(("TwoD_" + GetName() + "_selected_mc_reco").c_str(),  ("TwoD_"+GetName()).c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
    }

    //Histograms to be filled
    util::Categorized<Hist, int>* m_backgroundHists;
    util::Categorized<Hist, int>* m_SigIntTypeHists;
    util::Categorized<Hist, int>* m_SigTargetTypeHists;
    util::Categorized<Hist, int>* m_SigLeadBlobTypeHists;
    util::Categorized<Hist, int>* m_BkgIntTypeHists;
    util::Categorized<Hist, int>* m_BkgTargetTypeHists;
    util::Categorized<Hist, int>* m_BkgLeadBlobTypeHists;

    Hist* dataHist;  
    Hist* efficiencyNumerator;
    Hist* efficiencyDenominator;
    Hist* selectedSignalReco; //Effectively "true background subtracted" distribution for warping studies.                                                                                                         
                              //Also useful for a bakground breakdown plot that you'd use to start background subtraction studies.                                                                                 
    Hist* selectedMCReco; //Treat the MC CV just like data for the closure test                                           

    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
      //const char* name = GetName().c_str();
      dataHist = new Hist(("TwoD_" + GetName() + "_data").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), data_error_bands);
      //dataHist = new Hist(Form("_data_%s", name), name, GetBinVecX(), GetBinVecY(), data_error_bands);
    }

    //Repurposing to not write but just set the directory. The writing will be handled by TFile::Write.                                               
    void WriteData(TFile& file)
    {
      TString dirName = (TString)(fDirName);
      TDirectory* dir;
      dir = file.GetDirectory(dirName);
      if (dir == NULL){
        file.mkdir(dirName);
      }
      dir = file.GetDirectory(dirName);

      if (dataHist->hist) {
        dataHist->hist->SetDirectory(dir);
        //dataHist->hist->Write();
      }
    }

    void WriteMC(TFile& file)
    {
      SyncCVHistos();

      TString dirName = (TString)(fDirName);
      TDirectory* dir;
      dir = file.GetDirectory(dirName);
      if (dir == NULL){
        file.mkdir(dirName);
      }
      dir = file.GetDirectory(dirName);

      m_backgroundHists->visit([dir](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(dir);
                                      //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      m_SigIntTypeHists->visit([dir](Hist& categ)
			       {
				 categ.hist->SetDirectory(dir);
				 //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?
			       });

      m_SigTargetTypeHists->visit([dir](Hist& categ)
				  {
				    categ.hist->SetDirectory(dir);
				    //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?
				  });


      m_SigLeadBlobTypeHists->visit([dir](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(dir);
                                      //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?
				    });
      

      m_BkgIntTypeHists->visit([dir](Hist& categ)
			       {
				 categ.hist->SetDirectory(dir);
				 //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?
			       });

      m_BkgTargetTypeHists->visit([dir](Hist& categ)
				  {
				    categ.hist->SetDirectory(dir);
				    //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
				  });

      m_BkgLeadBlobTypeHists->visit([dir](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(dir);
                                      //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?
				    });

      /*
      if (dataHist->hist) {
	dataHist->hist->SetDirectory(dir);
	dataHist->hist->Write();
      }
      */

      if(efficiencyNumerator)
      {
        efficiencyNumerator->hist->SetDirectory(dir); //TODO: Can I get around having to call SetDirectory() this many times somehow?
        //efficiencyNumerator->hist->Write();
      }

      if(efficiencyDenominator)
      {
        efficiencyDenominator->hist->SetDirectory(dir);
        //efficiencyDenominator->hist->Write();
      }
      if(selectedSignalReco)
	{
	  selectedSignalReco->hist->SetDirectory(dir);
	  //selectedSignalReco->hist->Write();
	}

      if(selectedMCReco)
	{
	  selectedMCReco->hist->SetDirectory(dir);
	  selectedMCReco->hist->SetName(("TwoD_"+GetName() + "_data").c_str()); //Make this histogram look just like the data for closure tests
	}
    }

    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
      m_backgroundHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_SigIntTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_SigTargetTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_SigLeadBlobTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_BkgIntTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_BkgTargetTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_BkgLeadBlobTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });

      if(dataHist) dataHist->SyncCVHistos();
      if(efficiencyNumerator) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator) efficiencyDenominator->SyncCVHistos();
      if(selectedSignalReco) selectedSignalReco->SyncCVHistos();
      if(selectedMCReco) selectedMCReco->SyncCVHistos();
    }

    void SetDirectoryName(std::string name){fDirName = name;}
    void SetFillVar(bool fill){fFillVar = fill;}
    bool IsAnaVar(){return fAnaVar;}
    bool IsFill(){return fFillVar;}
    TString GetDirectoryName(){return fDirName;}
};

#endif //VARIABLE2D_H
