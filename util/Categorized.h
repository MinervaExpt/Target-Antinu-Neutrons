//File: Categorized.h
//Brief: HISTs can be Categorized<HIST> according to a NamedCategory<> to make
//       a separate HIST for each NamedCategory<> that can be Fill()ed based on
//       which NamedCategory<> a value belongs to.  Useful for putting together
//       stacked histograms that compare different "channels" with how they
//       contribute to the total histogram for a value.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef UTIL_CATEGORIZED_CPP
#define UTIL_CATEGORIZED_CPP

//Local includes
#include "util/SafeROOTName.h"

//c++ includes
#include <string>
#include <vector>
#ifndef __CINT__
#include <unordered_map>
#endif //__CINT__
#include <set>

namespace util
{
  std::map<int, std::string> TgtList = {{1,"Tgt1"},
                                        {2,"Tgt2"},
                                        {3,"Tgt3"},
                                        {4,"Tgt4"},
                                        {5,"Tgt5"},
                                        {6,"WaterTgt"},
                                        {10,"Plastic_US1"},
                                        {21,"Plastic_US2_DS1"},
                                        {32,"Plastic_US3_DS2"},
                                        {63,"Plastic_USWater_DS3"},
                                        {46,"Plastic_US4_DSWater"},
                                        {54,"Plastic_US5_DS4"},
                                        {0,"Plastic_DS5"}};

  std::map<int, std::map<int, std::string>> TgtCodeList = {{-1,{{1126,"Tgt1_Fe"},
								 {1182,"Tgt1_Pb"},
								 {1199,"Tgt1_Buffer"},
								 {2226,"Tgt2_Fe"},
								 {2282,"Tgt2_Pb"},
								 {2299,"Tgt2_Buffer"},
								 {3326,"Tgt3_Fe"},
								 {3382,"Tgt3_Pb"},
								 {3306,"Tgt3_C"},
								 {3399,"Tgt3_Buffer"},
								 {4482,"Tgt4_Pb"},
								 {5526,"Tgt5_Fe"},
								 {5582,"Tgt5_Pb"},
								 {5599,"Tgt5_Buffer"},
								 {6666,"WaterTgt"}}},
							    {1, {{1126,"Tgt1_Fe"},
								 {1182,"Tgt1_Pb"},
								 {1199,"Tgt1_Buffer"}}},
							    {2, {{2226,"Tgt2_Fe"},
								 {2282,"Tgt2_Pb"},
								 {2299,"Tgt2_Buffer"}}},
							    {3, {{3326,"Tgt3_Fe"},
								 {3382,"Tgt3_Pb"},
								 {3306,"Tgt3_C"},
								 {3399,"Tgt3_Buffer"}}},
							    {4, {{4482,"Tgt4_Pb"}}},
							    {5, {{5526,"Tgt5_Fe"},
								 {5582,"Tgt5_Pb"},
								 {5599,"Tgt5_Buffer"}}},
							    {6, {{6666,"WaterTgt"}}}};
  //Removed the below since I think the interstitial plastic is best handled separately from this list and doesn't need their own directory.
/*,
					    {1000,"Plastic_US1"},
					    {2100,"Plastic_US2_DS1"},
					    {3200,"Plastic_US3_DS2"},
					    {6300,"Plastic_USWater_DS3"},
					    {4600,"Plastic_US4_DSWater"},
					    {5400,"Plastic_US5_DS4"},
					    {7500,"Plastic_DS5"}};
							       */

  std::map<int, std::string> TgtTypeList = {{1, "Plastic"},
					    {200, "USPlastic"},
					    {300, "DSPlastic"},
					    {6, "C"},
					    {8, "Water"},
					    {26, "Fe"},
					    {82, "Pb"}};

  std::map<int, std::string> TgtTypeList2 = {{1, "Plastic"},
					    {6, "C"},
					    {8, "Water"},
					    {26, "Fe"},
					    {82, "Pb"}};

  std::map<int,int> DSTgtMap = {{0,5582},
				{54,4482},
				{46,6666},
				{63,3382},
				{32,2282},
				{21,1182},
				{10,-1}};

  std::map<int,int> USTgtMap = {{54,5582},
				{46,4482},
				{63,6666},
				{32,3382},
				{21,2282},
				{10,1182},
				{0,-1}};

  //Mapping from a set of values to a name.  Helper for constructing a Categorized<>
  template <class value_t>
  struct NamedCategory
  {
    std::vector<value_t> values;
    std::string name;
  };

  //A Categorized holds a total HIST along with a HIST for each category.
  //It works similarly to a Binned<>, but each entry either exactly matches
  //one CATEGORY or is put in the Other CATEGORY.
  //Its Fill() method takes a CATEGORY plus whatever ARGS HIST::Fill() takes.
  template <class HIST, class CATEGORY>
  //HIST is a Fill()able type that takes a c-string as its first constructor argument (like TH1D or Binned<TH1D>)
  //CATEGORY is hashable for std::unordered_map<> and has either a "name" member that can be added to a std::string
  //         or first and second members with first being hashable and second that can be added to a std::string.
  class Categorized
  {
    public:
      #ifndef __CINT__ //Hide "variadic template parameter" c++11 feature from CINT
      //CATEGORIES is an iterable container of objects like NamedCategory<>
      template <class ...HISTARGS>
      Categorized(const std::vector<NamedCategory<CATEGORY>>& categories, const std::string& baseName,
                  const std::string& axes, HISTARGS... args)
      {
        for(const auto& category: categories)
        {
          auto hist = new HIST(SafeROOTName(baseName + "_" + category.name).c_str(), (category.name + ";" + axes).c_str(), args...);
          for(const auto& value: category.values)
          {
            fCatToHist[value] = hist;
          }
        }

        fOther = new HIST((baseName + "_Other").c_str(), ("Other;" + axes).c_str(), args...);
      }

      //CATEGORY is a pointer to an object with a name() member variable
      template <class ...HISTARGS>
      Categorized(const std::vector<CATEGORY>& categories, const std::string& baseName,
                  const std::string& axes, HISTARGS... args)
      {
        for(const auto& catPtr: categories)
        {
          fCatToHist[catPtr] = new HIST(SafeROOTName(baseName + "_" + catPtr->name()).c_str(), (catPtr->name() + ";" + axes).c_str(), args...);
        }

        fOther = new HIST((baseName + "_Other"), ("Other;" + axes).c_str(), args...);
      }

      //CATEGORY is a pointer to an object with a name() member variable
      template <class ...HISTARGS>
      Categorized(const std::string& baseName, const std::string& axes,
                  const std::vector<CATEGORY>& categories, HISTARGS... args)
      {
        for(const auto& catPtr: categories)
        {
          fCatToHist[catPtr] = new HIST(SafeROOTName(baseName + "_" + catPtr->name()).c_str(), (catPtr->name() + ";" + axes).c_str(), args...);
        }

        fOther = new HIST((baseName + "_Other"), ("Other;" + axes).c_str(), args...);
      }

      //Use a std::map<> instead of CATEGORIES
      template <class ...HISTARGS>
      Categorized(const std::string& baseName, const std::string& axes,
                  const std::map<CATEGORY, std::string> categories, HISTARGS... args)
      {
        for(const auto& category: categories)
        {
          auto hist = new HIST(SafeROOTName(baseName + "_" + category.second).c_str(), (category.second + ";" + axes).c_str(), args...);
          fCatToHist[category.first] = hist;
        }

        fOther = new HIST((baseName + "_Other").c_str(), ("Other;" + axes).c_str(), args...);
      }

      //tag is useful for separating the possible overlapping of different Categorizations of the same thing. This exact function also has an extra argument dirName to help navigate the directory saving of variables which currently has its own risks...
      //THIS IS ONLY USEFUL FOR MY VARIABLE CATEGORIZATION RIGHT NOW. DOESN'T WORK AS WELL FOR HISTS...
      template <class ...HISTARGS>
      Categorized(const TString dirName, const TString tag, const bool isAnalysisVariable, const std::string& baseName, const std::string& axes,
                  const std::map<CATEGORY, std::string> categories, HISTARGS... args)
      {
        for(const auto& category: categories)
        {
          auto hist = (dirName != "") ? new HIST(dirName+"/"+tag+"_"+(TString)category.second, isAnalysisVariable, SafeROOTName(baseName + "_" + category.second).c_str(), (category.second + ";" + axes).c_str(), args...) : new HIST(tag+"_"+(TString)category.second, isAnalysisVariable, SafeROOTName(baseName + "_" + category.second).c_str(), (category.second + ";" + axes).c_str(), args...);
          fCatToHist[category.first] = hist;
        }

	auto histOther = (dirName != "") ? new HIST(dirName+"/"+tag+"_Other", isAnalysisVariable, (baseName + "_Other").c_str(), ("Other;" + axes).c_str(), args...) : new HIST(tag+"_Other", isAnalysisVariable, (baseName + "_Other").c_str(), ("Other;" + axes).c_str(), args...);
	fOther = histOther;
      }

      //tag is useful for separating the possible overlapping of different Categorizations of the same thing. This exact function also has an extra argument dirName to help navigate the directory saving of variables which currently has its own risks...
      //THIS IS ONLY USEFUL FOR MY VARIABLE CATEGORIZATION RIGHT NOW. DOESN'T WORK AS WELL FOR HISTS...
      template <class ...HISTARGS>
      Categorized(const TString dirName, const TString tag, const bool isAnalysisVariable, const std::string& baseName,
                  const std::map<CATEGORY, std::string> categories, HISTARGS... args)
      {
        for(const auto& category: categories)
        {
          auto hist = (dirName != "") ? new HIST(dirName+"/"+tag+"_"+(TString)category.second, isAnalysisVariable, SafeROOTName(baseName + "_" + category.second).c_str(), args...) : new HIST(tag+"_"+(TString)category.second, isAnalysisVariable, SafeROOTName(baseName + "_" + category.second).c_str(), args...);
          fCatToHist[category.first] = hist;
        }

	auto histOther = (dirName != "") ? new HIST(dirName+"/"+tag+"_Other", isAnalysisVariable, (baseName + "_Other").c_str(), args...) : new HIST(tag+"_Other", isAnalysisVariable, (baseName + "_Other").c_str(), args...);
	fOther = histOther;
      }
      #endif //__CINT__

      HIST& operator [](const CATEGORY& cat) const
      {
        //Find out whether category is kept track of separately
        const auto found = fCatToHist.find(cat);
        if(found == fCatToHist.end()) return *fOther; //If not, lump this entry in with other uncategorized entries
        return *found->second; //If so, return its category
        
      }

      //Apply a callable object, of type FUNC, to each histogram this object manages.
      //FUNC takes only a reference to the histogram as argument.
      template <class FUNC>
      void visit(FUNC&& func)
      {
        #ifndef __CINT__ //Hide "auto" c++11 feature from CINT
        //Make sure that each histogram is normalized exactly once
        std::set<HIST*> histsNormalized;
        for(auto& category: fCatToHist)
        {
          if(histsNormalized.count(category.second) == 0)
          {
            func(*(category.second));
            histsNormalized.insert(category.second);
          }
        }

        func(*fOther);
        #endif //__CINT__
      }

      //TODO: I think this is needed for nested Categorized<Categorized<HistWrapper<>, >, > because
      //      HistWrapper calls SetDirectory(0) on its MnvH1D.
      /*void SetDirectory(TDirectory* dir)
      {
        visit([dir](auto& hist) { util::detail::set<HIST>::dir(hist, *dir); });
      }*/

    private:
      //All HISTs are referred to as observer pointers for compatability with TH1s created in a TFile.
      //The TFile is responsible for deleting them.
      #ifndef __CINT__ //Hide "std::unordered_map<>" c++11 feature from CINT
      std::unordered_map<CATEGORY, HIST*> fCatToHist;
      #endif //__CINT__
      HIST* fOther; //All entries that don't fit in any other CATEGORY end up in this HIST
  };
}

#endif //UTIL_CATEGORIZED_CPP
