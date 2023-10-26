#include "TROOT.h"
#include "TString.h"

#define PI 3.14159
#define TWO_PI 6.28318

using namespace std;

void getStringFromTag(string tag,std::istringstream &strStream, string inputStr  , string &var)
{
    std::string field;
    if(inputStr.compare(tag)==0){
           getline(strStream, field);
           var =  field;
           cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}


void getBoolFromTag(string tag,std::istringstream &strStream, string inputStr  , Bool_t &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
           getline(strStream, field);
           var =  std::atoi(field.c_str()) > 0 ? 1 : 0;
           cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}

void getFloatFromTag(string tag,std::istringstream &strStream, string inputStr  , Float_t &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
             getline(strStream, field);
             var=std::atof(field.c_str());
             cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}

template <typename T> void getIntFromTag(string tag,std::istringstream &strStream, string inputStr  , T &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
             getline(strStream, field);
             var=std::atof(field.c_str());
             cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}



void getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose)
{
	
    bool cfgModeFlag=false;
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
	string line;
    
    // getting flists
    cfgFile.clear();
	cfgFile.seekp(ios::beg);
    cfgModeFlag=false;
	int nItems(0);
    while(std::getline(cfgFile,line))
	{
	   if(line==beginTag) {cfgModeFlag=true;continue;}
	   if(line==endTag) {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   vecList.push_back(line);
	   nItems++;
    }

    if(verbose)
    {
       std::cout<<" Added "<<nItems<<" between "<<beginTag<<" and "<< endTag<<"\n";
       for( auto &name : vecList)
       {
           std::cout<<"\t"<<name<<"\n";
       }
    }

}

Double_t getDPHI( Double_t phi1, Double_t phi2) {

  Double_t dphi = phi1 - phi2;

  while( dphi > PI)
        dphi-= TWO_PI;
  while( dphi < -PI)
        dphi += TWO_PI;
  //std::cout<<"phi1  : "<<phi1<<" , phi2 : "<<phi2<<" dphi : "<<dphi<<"\n";
  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 "<< dphi << endl;
  }

  return TMath::Abs(dphi);
  //return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2)
{

    return TMath::Abs(eta1 - eta2);
}

Double_t getDR( Double_t eta1, Double_t phi1,Double_t eta2 ,Double_t phi2) {

    Double_t de = getDETA(eta1,eta2);
    Double_t dp = getDPHI(phi1,phi2);
    return sqrt(de*de + dp*dp);
}
