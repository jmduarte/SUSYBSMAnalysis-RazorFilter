#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/BTauReco/interface/JetTag.h"


#include <iostream>
#include <string>

#include <TTree.h>
#include <TFile.h>

#include <fstream>

#include "TVector3.h"
#include "TLorentzVector.h"

class RazorFilter : public edm::EDFilter 
{
public: 
  explicit RazorFilter(const edm::ParameterSet&);
  ~RazorFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void computeHemispheres(std::auto_ptr<std::vector<math::XYZTLorentzVector> >& hlist, const std::vector<math::XYZTLorentzVector>& goodJets);
  static double calcMR(TLorentzVector ja, TLorentzVector jb);
  static double calcRsq(double MR, TLorentzVector ja, TLorentzVector jb, edm::Handle<reco::PFMETCollection> inputMet);
      
  edm::InputTag jetInputTag_;
  edm::InputTag metInputTag_;
  edm::InputTag btagInputTag_;

  PFJetIDSelectionFunctor m_jetIDfunc;

  std::string rootFileName_;
  double minJetPt_;
  double maxJetEta_;
 
  TFile* rootFile_;
  TTree*  razorTree_;

  std::ofstream csvOut_;
  std::string csvFileName_;

  double MR;
  double Rsq;
  double HT;
  double MET;
  int runNum;
  int lumiNum;
  int eventNum;
  int nJets;
  int nBJets;

};

RazorFilter::RazorFilter(const edm::ParameterSet& iConfig)
  : jetInputTag_(iConfig.getParameter<edm::InputTag>("jetInputTag")),
    metInputTag_(iConfig.getParameter<edm::InputTag>("metInputTag")),
    btagInputTag_(iConfig.getParameter<edm::InputTag>("btagInputTag")),
    m_jetIDfunc(PFJetIDSelectionFunctor::FIRSTDATA,PFJetIDSelectionFunctor::LOOSE),
    rootFileName_(iConfig.getParameter<std::string>("rootFileName")),
    minJetPt_(iConfig.getParameter<double>("minJetPt")),
    maxJetEta_(iConfig.getParameter<double>("maxJetEta")),
    csvFileName_(iConfig.getParameter<std::string>("csvFileName"))
{
  rootFile_ = new TFile(rootFileName_.c_str(), "RECREATE");

  
  razorTree_ = new TTree("razorTree",
			 "razorTree");

  csvOut_.open(csvFileName_.c_str());

}


RazorFilter::~RazorFilter()
{}

bool
RazorFilter::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::PFJetCollection> jets;
  event.getByLabel(jetInputTag_, jets);

  edm::Handle<reco::PFMETCollection> met;
  event.getByLabel(metInputTag_, met);

  edm::Handle<reco::JetTagCollection> btags;
  event.getByLabel(btagInputTag_, btags);
  
  if ( ! (jets.isValid() && met.isValid() && btags.isValid()) )
  {
    std::cerr<<"RazorFilter: invalid collection"<<std::endl;
    return false;
  }
  
  nJets = 0;
  nBJets = 0;
  MR = -1;
  Rsq = -1;
  MET = -1;
  HT = 0;
  runNum = -1;
  lumiNum = -1;
  eventNum = -1;
  
  
  // Only examine if there are at least 2 jets
  if ( jets->size() < 2 )
    return false;

  std::vector<math::XYZTLorentzVector> goodJets;
   
  //loop over jets to find those passing pt, eta cuts
  double pt, eta;
  for ( reco::PFJetCollection::const_iterator it = jets->begin(), end = jets->end(); 
        it != end; ++it)
  {
    pt = (*it).pt();
    eta = (*it).eta();
    if (pt > minJetPt_ && fabs(eta) < maxJetEta_){
      pat::strbitset ret = m_jetIDfunc.getBitTemplate();
      ret.set(false);
      bool isGoodJet = m_jetIDfunc((*it),ret);
      if(isGoodJet) {
	goodJets.push_back((*it).p4());
	HT += pt;
	nJets++;
      }
    }
  }

  //loop over btags to find those passing pt, eta cuts
  for ( reco::JetTagCollection::const_iterator it = btags->begin(), end = btags->end(); 
        it != end; ++it)
  {
    pt = (*it).first->pt();
    eta = (*it).first->eta();
    if (pt > minJetPt_ && fabs(eta) < maxJetEta_){
      //if ( (*it).second > 0.898 ) // Tight working point for combined secondary vertex algorithm
      if ( (*it).second > 0.679 ) // Medium working point for combined secondary vertex algorithm
      //if ( (*it).second > 0.244 ) // Loose working point for combined secondary vertex algorithm
	nBJets++;
    }
  }

  // only continue if there are at least 2 jets
  if (nJets < 2)
    return false;
  
  //compute the hemispheres which minimize the invariant masses summed in quadrature
  std::auto_ptr<std::vector<math::XYZTLorentzVector> > hemispheres(new std::vector<math::XYZTLorentzVector> );
  this->computeHemispheres(hemispheres,goodJets);

  if (hemispheres->size() < 2)
    return false;

   TLorentzVector ja(hemispheres->at(0).x(),hemispheres->at(0).y(),hemispheres->at(0).z(),hemispheres->at(0).t());
   TLorentzVector jb(hemispheres->at(1).x(),hemispheres->at(1).y(),hemispheres->at(1).z(),hemispheres->at(1).t());

   // compute the razor variables
   MR = calcMR(ja,jb);
   Rsq  = calcRsq(MR,ja,jb,met);

   // record the MET
   MET = (met->front()).pt();
   
   runNum = event.id().run();
   lumiNum = event.id().luminosityBlock();
   eventNum = event.id().event();
   
   // write out to csv file
   csvOut_<< event.id().run() <<","<< event.id().luminosityBlock() <<","<< event.id().event() <<","
	  << MR << "," << Rsq << "," 
	  << ja.E() << "," << ja.Px() << "," << ja.Py() << "," << ja.Pz() << "," 
	  << jb.E() << "," << jb.Px() << "," << jb.Py() << "," << jb.Pz() << "," 
	  << HT << "," << MET << "," << nJets << "," << nBJets << std::endl;

   // fill the tree
   razorTree_->Fill();
  
  return true;
}

void 
RazorFilter::beginJob()
{
  csvOut_<<"Run,Lumi,Event,MR,Rsq,E1,px1,py1,pz1,E2,px2,py2,pz2,HT,MET,nJets,nBJets"<<std::endl;

  razorTree_->Branch("runNum", &runNum, "runNum/I");
  razorTree_->Branch("lumiNum", &lumiNum, "lumiNum/I");
  razorTree_->Branch("eventNum", &eventNum, "eventNum/I");
  razorTree_->Branch("MR", &MR, "MR/D");
  razorTree_->Branch("Rsq", &Rsq, "Rsq/D");
  razorTree_->Branch("HT", &HT, "HT/D");
  razorTree_->Branch("MET", &MET, "MET/D");
  razorTree_->Branch("nJets", &nJets, "nJets/I");
  razorTree_->Branch("nBJets", &nBJets, "nBJets/I");
}

void 
RazorFilter::endJob() 
{
  rootFile_->cd();
  razorTree_->Write();
  rootFile_->Close();

  csvOut_.close();
}



void
RazorFilter::computeHemispheres(std::auto_ptr<std::vector<math::XYZTLorentzVector> >& hlist, const std::vector<math::XYZTLorentzVector>& goodJets){
  using namespace math;
  using namespace reco;
  XYZTLorentzVector j1R(0.1, 0., 0., 0.1);
  XYZTLorentzVector j2R(0.1, 0., 0., 0.1);
  int nJets = goodJets.size();

  if(nJets<2){ // put empty hemispheres if not enough jets
    hlist->push_back(j1R);
    hlist->push_back(j2R);
    return;
  }
  unsigned int N_comb = pow(2,nJets); // compute the number of combinations of jets possible
  //Make the hemispheres
  double M_minR = 999999999999.0;
  unsigned int j_count;
  for (unsigned int i = 0; i < N_comb; i++) {       
    XYZTLorentzVector j_temp1, j_temp2;
    unsigned int itemp = i;
    j_count = N_comb/2;
    unsigned int count = 0;
    while (j_count > 0) {
      if (itemp/j_count == 1){
	j_temp1 += goodJets.at(count);
      } else {
	j_temp2 += goodJets.at(count);
      }
      itemp -= j_count * (itemp/j_count);
      j_count /= 2;
      count++;
    }
    double M_temp = j_temp1.M2() + j_temp2.M2();
    if (M_temp < M_minR) {
      M_minR = M_temp;
      j1R = j_temp1; 
      j2R = j_temp2; 
    }
  }

  //return the hemispheres in decreasing order of pt
  if(j1R.Pt() > j2R.Pt()){
    hlist->push_back(j1R);
    hlist->push_back(j2R);
  } else {
    hlist->push_back(j2R);
    hlist->push_back(j1R);
  }
  
  return;
}

double 
RazorFilter::calcMR(TLorentzVector ja, TLorentzVector jb){
  if(ja.Pt()<=0.1) return -1;

  ja.SetPtEtaPhiM(ja.Pt(),ja.Eta(),ja.Phi(),0.0);
  jb.SetPtEtaPhiM(jb.Pt(),jb.Eta(),jb.Phi(),0.0);
  
  if(ja.Pt() > jb.Pt()){
    TLorentzVector temp = ja;
    ja = jb;
    jb = temp;
  }
  
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();

  double MR = sqrt((A+B)*(A+B)-(az+bz)*(az+bz));

  return MR;
}

double 
RazorFilter::calcRsq(double MR, TLorentzVector ja, TLorentzVector jb, edm::Handle<reco::PFMETCollection> inputMet){
  //now we can calculate MTR
  TVector3 met;
  met.SetPtEtaPhi((inputMet->front()).pt(),0.0,(inputMet->front()).phi());

  double MTR = sqrt(0.5*(met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect())));
  
  //filter events
  return MTR*MTR/MR/MR; //Rsq
  
}


DEFINE_FWK_MODULE(RazorFilter);
