	
	/*
	  m_SigLeadBlobTypeHists->visit([dir](Hist& categ)
	  {
	  categ.hist->SetDirectory(dir);
	  //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
	  });
	*/
	
	m_BkgTargetTypeHists->visit([dir](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(dir);
                                      //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });
	/*
	  m_BkgLeadBlobTypeHists->visit([dir](Hist& categ)
	  {
	  categ.hist->SetDirectory(dir);
	  //categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
	  });
	*/
      }

      if(efficiencyNumerator && fAnaVar)
      {
        efficiencyNumerator->hist->SetDirectory(dir); //TODO: Can I get around having to call SetDirectory() this many times somehow?
        //efficiencyNumerator->hist->Write();
      }

      if(efficiencyDenominator && fAnaVar)
      {
        efficiencyDenominator->hist->SetDirectory(dir);
        //efficiencyDenominator->hist->Write();
      }

      if(migration && fAnaVar)
      {
        migration->hist->SetDirectory(dir); 
        //migration->hist->Write();
      }

      if(selectedSignalReco)
      {
        selectedSignalReco->hist->SetDirectory(dir);
        //selectedSignalReco->hist->Write();
      }

      if(selectedMCReco)
      {
        selectedMCReco->hist->SetDirectory(dir);
        selectedMCReco->hist->SetName((GetName() + "_data").c_str()); //Make this histogram look just like the data for closure tests
      }
    }

    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
      m_backgroundHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_BkgIntTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });

      if (fAllBreakdowns){
	m_SigIntTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
	m_SigTargetTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
	//m_SigLeadBlobTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });

	m_BkgTargetTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
	//m_BkgLeadBlobTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      }

      if(dataHist) dataHist->SyncCVHistos();
      if(efficiencyNumerator && fAnaVar) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator && fAnaVar) efficiencyDenominator->SyncCVHistos();
      if(selectedSignalReco) selectedSignalReco->SyncCVHistos();
      if(selectedMCReco) selectedMCReco->SyncCVHistos();
      if(migration && fAnaVar) migration->SyncCVHistos();
    }

    void SetDirectoryName(std::string name){fDirName = name;}
    void SetFillVar(bool fill){fFillVar = fill;}
    bool IsAnaVar(){return fAnaVar;}
    bool IsFill(){return fFillVar;}
    bool IsBroken(){return fAllBreakdowns;}
    TString GetDirectoryName(){return fDirName;}
};

#endif //VARIABLE_H
