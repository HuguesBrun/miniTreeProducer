#define miniTreeBuilder_cxx
// The class definition in miniTreeBuilder.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("miniTreeBuilder.C")
// root> T->Process("miniTreeBuilder.C","some options")
// root> T->Process("miniTreeBuilder.C+")
//

#include "miniTreeBuilder.h"
#include <TH2.h>
#include <TStyle.h>

using namespace std;

void miniTreeBuilder::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
    
    myIterator = 0;

}

void miniTreeBuilder::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t miniTreeBuilder::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either miniTreeBuilder::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
   // fChain->GetEntry(entry);
    fChain->GetTree()->GetEntry(entry);
    myIterator++;
   // if (entry%10000==0) cout << "entry=" << entry << endl;
    if (myIterator%100==0) cout << "myIterator=" << myIterator << endl;
    int nbElectrons = T_Elec_Eta->size();
    //cout << "they are " << nbElectrons << " in the event" << endl;
 //   if (nbElectrons<2) return kTRUE;
   // for (int i = 0 ; i < nbElectrons ; i++){
        //cout << "id=" << T_Gen_Elec_PDGid->at(i) << endl;
        //cout << "id=" << T_Elec_Pt->at(i) << endl;
        for (int j = 0 ; j < nbElectrons ; j++){
  //          if (i==j) continue;
            if (fabs(T_Gen_Elec_PDGid->at(j))!=11) continue;
            //if (!(passTrigSelection(j))) continue;
            pt = T_Elec_Pt->at(j);
            eta = T_Elec_Eta->at(j);
            absEta = fabs(T_Elec_Eta->at(j));
            phi = T_Elec_Phi->at(j);
            isFromTau = T_Gen_Elec_fromTAU->at(j);
            int refTrigger = findTriggerMatched(j);
            int refTriggerElec = findTriggerMatchedElec(j);
            saveInTheTree(j, refTrigger, refTriggerElec);
        }
//    }
    

   if (myIterator>100000) Abort("the end !!!");
   return kTRUE;
}

bool miniTreeBuilder::passNewTrigSelection(int theIte){
    if (T_Elec_TKiso->at(theIte)/T_Elec_Pt->at(theIte)>0.19) return false;
    
    if (fabs(T_Elec_SC_Eta->at(theIte))<1.479){
        if (T_Elec_noZSsee->at(theIte)>0.0125) return false;
        if (fabs(T_Elec_deltaEtaIn->at(theIte))>0.009) return false;
        if (fabs(T_Elec_deltaPhiIn->at(theIte))>0.06) return false;
        if (T_Elec_HtoE->at(theIte)>0.14) return false;
        if ((T_Elec_ECALiso->at(theIte))/T_Elec_Pt->at(theIte)>0.14) return false;
        if (T_Elec_HCALiso->at(theIte)/T_Elec_Pt->at(theIte)>0.45) return false;
  //      if (T_Elec_TKiso->at(theIte)/T_Elec_Pt->at(theIte)>0.19) return false;


    }
    else  {
        if (T_Elec_noZSsee->at(theIte)>0.034) return false;
        if (fabs(T_Elec_deltaEtaIn->at(theIte))>0.014) return false;
        if (fabs(T_Elec_deltaPhiIn->at(theIte))>0.09) return false;
        if (T_Elec_HtoE->at(theIte)>0.35) return false;
        if ((T_Elec_ECALiso->at(theIte))/T_Elec_Pt->at(theIte)>0.73) return false;
        if (T_Elec_HCALiso->at(theIte)/T_Elec_Pt->at(theIte)>0.45) return false;

    }
    return true;
}

/*bool miniTreeBuilder::passNewTrigSelection(int theIte){
    if (T_Elec_TKiso->at(theIte)/T_Elec_Pt->at(theIte)>10) return false;
    
    if (fabs(T_Elec_SC_Eta->at(theIte))<1.479){
        if (T_Elec_noZSsee->at(theIte)>0.0185) return false;
        if (fabs(T_Elec_deltaEtaIn->at(theIte))>0.035) return false;
        if (fabs(T_Elec_deltaPhiIn->at(theIte))>0.026) return false;
        if (T_Elec_HtoE->at(theIte)>0.07) return false;
        if ((T_Elec_ECALiso->at(theIte))/T_Elec_Pt->at(theIte)>0.13) return false;
        if (T_Elec_HCALiso->at(theIte)/T_Elec_Pt->at(theIte)>0.34) return false;
        if (T_Elec_TKiso->at(theIte)/T_Elec_Pt->at(theIte)>10) return false;
        
        
    }
    else  {
        if (T_Elec_noZSsee->at(theIte)>0.039) return false;
        if (fabs(T_Elec_deltaEtaIn->at(theIte))>0.035) return false;
        if (fabs(T_Elec_deltaPhiIn->at(theIte))>0.063) return false;
        if (T_Elec_HtoE->at(theIte)>0.075) return false;
        if ((T_Elec_ECALiso->at(theIte))/T_Elec_Pt->at(theIte)>0.73) return false;
        if (T_Elec_HCALiso->at(theIte)/T_Elec_Pt->at(theIte)>0.2) return false;
        
    }
    return true;
}*/

bool miniTreeBuilder::passTrigSelection(int theIte){
    if (T_Elec_HCALiso->at(theIte)/T_Elec_Pt->at(theIte)>0.2) return false;
    if (T_Elec_TKiso->at(theIte)/T_Elec_Pt->at(theIte)>0.2) return false;
    
    if (fabs(T_Elec_SC_Eta->at(theIte))<1.479){
        if (T_Elec_noZSsee->at(theIte)>0.01) return false;
        if (fabs(T_Elec_deltaEtaIn->at(theIte))>0.007) return false;
        if (fabs(T_Elec_deltaPhiIn->at(theIte))>0.15) return false;
        if (T_Elec_HtoE->at(theIte)>0.12) return false;
        if ((T_Elec_ECALiso->at(theIte)-1)/T_Elec_Pt->at(theIte)>0.2) return false;
    }
    else  {
        if (T_Elec_noZSsee->at(theIte)>0.03) return false;
        if (fabs(T_Elec_deltaEtaIn->at(theIte))>0.009) return false;
        if (fabs(T_Elec_deltaPhiIn->at(theIte))>0.1) return false;
        if (T_Elec_HtoE->at(theIte)>0.1) return false;
        if ((T_Elec_ECALiso->at(theIte))/T_Elec_Pt->at(theIte)>0.2) return false;
    }
    return true;
}

void miniTreeBuilder::saveInTheTree(int theIte, int triggerIte, int triggerIteElec){
    rho = T_Event_Rho->at(0);
 energy = T_Elec_Energy->at(theIte);
 nLost = T_Elec_nLost->at(theIte);
 //nHits = T_Elec_nHits->at(theIte);
 kfchi2 = T_Elec_kfchi2->at(theIte);
 kfhits = T_Elec_kfhits->at(theIte);
 gsfhits = T_Elec_gsfhits->at(theIte);
 gsfChi2 = T_Elec_gsfchi2->at(theIte);
 fBrem = T_Elec_fbrem->at(theIte);
 Dist = T_Elec_Dist->at(theIte);
 Dcot = T_Elec_Dcot->at(theIte);
 D0 = T_Elec_D0->at(theIte);
 DZ =T_Elec_Dz->at(theIte);
 ip3d = T_Elec_ip3d->at(theIte);
 ip3ds = T_Elec_ip3ds->at(theIte);
 detacalo = T_Elec_detacalo->at(theIte);
 eledeta = T_Elec_eledeta->at(theIte);
 eledphi = T_Elec_eledphi->at(theIte);
 dphicalo = T_Elec_dphicalo->at(theIte);
 deltaPhiIn = fabs(T_Elec_deltaPhiIn->at(theIte));
 deltaEtaIn = fabs(T_Elec_deltaEtaIn->at(theIte));
 EoP = T_Elec_EoP->at(theIte);
 ESeedoP = T_Elec_ESeedoP->at(theIte);
 ESeedoPout = T_Elec_ESeedoPout->at(theIte);
 EEleoPout = T_Elec_EEleoPout->at(theIte);
 IoEmIoP = T_Elec_IoEmIoP->at(theIte);
 SC_Et = T_Elec_SC_Et->at(theIte);
 SC_Eta = T_Elec_SC_Eta->at(theIte);
 SC_Phi = T_Elec_SC_Phi->at(theIte);
 SC_RawEnergy = T_Elec_SC_RawEnergy->at(theIte);
 EcalEnergy = T_Elec_EcalEnergy->at(theIte);
 EsEnergy = T_Elec_EsEnergy->at(theIte);
 PreShowerOverRaw = T_Elec_PreShowerOverRaw->at(theIte);
 NClusters = T_Elec_NClusters->at(theIte);
 EtaSeed = T_Elec_EtaSeed->at(theIte);
 PhiSeed = T_Elec_PhiSeed->at(theIte);
 ESeed = T_Elec_ESeed->at(theIte);
 HtoE = T_Elec_HtoE->at(theIte);
 EmaxSeed = T_Elec_EmaxSeed->at(theIte);
 EtopSeed = T_Elec_EtopSeed->at(theIte);
 EbottomSeed = T_Elec_EbottomSeed->at(theIte);
 EleftSeed = T_Elec_EleftSeed->at(theIte);
 ErightSeed = T_Elec_ErightSeed->at(theIte);
 E2ndSeed = T_Elec_E2ndSeed->at(theIte);
 E2x5RightSeed = T_Elec_E2x5RightSeed->at(theIte);
 E2x5LeftSeed = T_Elec_E2x5LeftSeed->at(theIte);
 E2x5TopSeed = T_Elec_E2x5TopSeed->at(theIte);
 E2x5BottomSeed = T_Elec_E2x5BottomSeed->at(theIte);
 E2x5MaxSeed =T_Elec_E2x5MaxSeed->at(theIte);
// E1x5Seed =T_Elec_E1x5Seed->at(theIte);
// E2x2Seed =T_Elec_E2x2Seed->at(theIte);
 E3x3Seed =T_Elec_E3x3Seed->at(theIte);
 E5x5Seed =T_Elec_E5x5Seed->at(theIte);
 see =T_Elec_see->at(theIte);
 spp = T_Elec_spp->at(theIte);
 sep = T_Elec_sep->at(theIte);
 etawidth = T_Elec_etawidth->at(theIte);
 phiwidth = T_Elec_phiwidth->at(theIte);
 e1x5e5x5 = T_Elec_e1x5e5x5->at(theIte);
 s9e25 = T_Elec_s9e25->at(theIte);
 R9 = T_Elec_R9->at(theIte);
 MatchConv = T_Elec_MatchConv->at(theIte);
 EcalDriven = T_Elec_EcalDriven->at(theIte);
 noZSsee = T_Elec_noZSsee->at(theIte);
 noZSspp = T_Elec_noZSspp->at(theIte);
 noZSsep = T_Elec_noZSsep->at(theIte);
 noZSr9 = T_Elec_noZSr9->at(theIte);
 noZSe1x5 = T_Elec_noZSe1x5->at(theIte);
 noZSe2x5MaxSeed = T_Elec_noZSe2x5MaxSeed->at(theIte);
 noZSe5x5 = T_Elec_noZSe5x5->at(theIte);
 noZSe1x5e5x5  = 1-1.*(noZSe1x5/noZSe5x5);
 ECALiso = T_Elec_ECALiso->at(theIte);
 HCALiso = T_Elec_HCALiso->at(theIte);
 TKiso = T_Elec_TKiso->at(theIte);
    relatECALiso = (fabs(eta)<1.479) ? (ECALiso-rho*0.083)/pt: (ECALiso-rho*0.046)/pt;
    relatHCALiso = (fabs(eta)<1.479) ? (HCALiso-rho*0.030)/pt: (HCALiso-rho*0.076)/pt;
    relatTKiso = TKiso/pt;
    passVeto = T_Elec_Veto->at(theIte);
    passLoose = T_Elec_Loose->at(theIte);
    passMedium = T_Elec_Medium->at(theIte);
    passTight = T_Elec_Tight->at(theIte);
 nbBC = T_Elec_nbBC->at(theIte);
 oldTrigering = passTrigSelection(theIte);
    newTrigering = passNewTrigSelection(theIte);
 PFisoChargedHadron = T_Elec_chargedHadronIso->at(theIte);
 PFisoChargedHadronAll = T_Elec_allChargedHadronIso->at(theIte);
 PFisoNeutralHadron = T_Elec_neutralHadronIso->at(theIte);
 PFisoPhoton = T_Elec_photonIso->at(theIte);
 matchedToTrigger = isPassingBit(1,T_Elec_TriggerLeg->at(theIte))||isPassingBit(2,T_Elec_TriggerLeg->at(theIte));
   // cout << "ite1=" << triggerIte << " ite2=" << triggerIteElec << endl;
    if (triggerIte == -1){
     matchedToTriggerFilter =0;
     trig_eta=-1;
     trig_phi=-1;
     trig_pt=-1;
     trig_sigEta=-1;
     trig_isoEcal=-1;
     trig_isoHcal=-1;
     trig_HoE=-1;
     /*trig_dPhi=-1;
     trig_dEta=-1;
     trig_isoTracker=-1;*/
    }
    else {
        matchedToTriggerFilter = 1;
        trig_eta=T_Trig_Eta->at(triggerIte);
        trig_phi=T_Trig_Phi->at(triggerIte);
        trig_pt=T_Trig_Pt->at(triggerIte);
        trig_sigEta=T_Trig_sigEta->at(triggerIte);
        trig_isoEcal=T_Trig_isoECAL->at(triggerIte);
        trig_isoHcal=T_Trig_isoHCAL->at(triggerIte);
        trig_HoE=T_Trig_HoE->at(triggerIte)/trig_pt;
       /* trig_dPhi=T_Trig_dPhi->at(triggerIte);
        trig_dEta=T_Trig_dEta->at(triggerIte);
        trig_isoTracker=T_Trig_isoTracker->at(triggerIte);*/
    }

    if (triggerIteElec == -1){
        matchedToTriggerFilter =0;
        trig_eta=-1;
        trig_phi=-1;
        trig_pt=-1;
       /* trig_sigEta=-1;
        trig_isoEcal=-1;
        trig_isoHcal=-1;
        trig_HoE=-1;*/
        trig_dPhi=-1;
        trig_dEta=-1;
        trig_isoTracker=-1;
    }
    else {
        matchedToTriggerFilter = 1;
        trig_eta=T_Trig_Eta->at(triggerIteElec);
        trig_phi=T_Trig_Phi->at(triggerIteElec);
        trig_pt=T_Trig_Pt->at(triggerIteElec);
       /* trig_sigEta=T_Trig_sigEta->at(triggerIte);
        trig_isoEcal=T_Trig_isoECAL->at(triggerIte);
        trig_isoHcal=T_Trig_isoHCAL->at(triggerIte);
        trig_HoE=T_Trig_HoE->at(triggerIte);*/
        trig_dPhi= (T_Trig_dPhi->at(triggerIteElec)<1 ? T_Trig_dPhi->at(triggerIteElec): 1);
        trig_dEta=(T_Trig_dEta->at(triggerIteElec)<1 ? T_Trig_dEta->at(triggerIteElec): 1) ;
        trig_isoTracker=T_Trig_isoTracker->at(triggerIteElec);
    }
    
    trig_isoRelatEcal = (trig_isoEcal/trig_pt<0.9 ? trig_isoEcal/trig_pt : 0.9);
    trig_isoRelatHcal = (trig_isoHcal/trig_pt<0.9 ? trig_isoHcal/trig_pt: 0.9);
    trig_isoRelatTracker = (trig_isoTracker/trig_pt<0.9 ? trig_isoTracker/trig_pt : 0.9);
    
    elec_passHLT = willPassTheHLT(trig_eta,  trig_pt, trig_sigEta, trig_isoEcal,  trig_HoE,  trig_isoHcal,  trig_dPhi,  trig_dEta,  trig_isoTracker);
myTree->Fill();

}

void miniTreeBuilder::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void miniTreeBuilder::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
    myFile->Write();
    myFile->Close();
}

int miniTreeBuilder::findTriggerMatched(int theIte){
    int nbTriggerObject = T_Trig_Eta->size();
    float elecPhi = T_Elec_Phi->at(theIte);
    float elecEta = T_Elec_Eta->at(theIte);
    
    
    for (int m=0 ; m < nbTriggerObject ; m++){
        float deltaR = sqrt(pow(elecEta-T_Trig_Eta->at(m),2)+ pow(acos(cos((elecPhi-T_Trig_Phi->at(m)))),2));
       // cout << "deltaR=" << deltaR << endl;
        if ((deltaR<0.1)&&((T_Trig_Leg->at(m)==0)||(T_Trig_Leg->at(m)==1))) {
            //cout << "electron= eta=" << elecEta << " phi=" << elecPhi << " pt=" << T_Elec_Pt->at(theIte) << endl;
            //cout << "trigger= eta=" << T_Trig_Eta->at(m) << " phi=" << T_Trig_Phi->at(m) << " pt=" << T_Trig_Pt->at(m) << endl;
            return m;
        }
    }
    return -1;
    
}

int miniTreeBuilder::findTriggerMatchedElec(int theIte){
    int nbTriggerObject = T_Trig_Eta->size();
    float elecPhi = T_Elec_Phi->at(theIte);
    float elecEta = T_Elec_Eta->at(theIte);
    
    
    for (int m=0 ; m < nbTriggerObject ; m++){
        float deltaR = sqrt(pow(elecEta-T_Trig_Eta->at(m),2)+ pow(acos(cos((elecPhi-T_Trig_Phi->at(m)))),2));
        // cout << "deltaR=" << deltaR << endl;
        if ((deltaR<0.1)&&((T_Trig_Leg->at(m)==2)||(T_Trig_Leg->at(m)==3))) {
            //cout << "electron= eta=" << elecEta << " phi=" << elecPhi << " pt=" << T_Elec_Pt->at(theIte) << endl;
            //cout << "trigger= eta=" << T_Trig_Eta->at(m) << " phi=" << T_Trig_Phi->at(m) << " pt=" << T_Trig_Pt->at(m) << endl;
            //cout << "m2=" << m << endl;
            return m;
        }
    }
    return -1;
    
}

bool miniTreeBuilder::willPassTheHLT(float eta, float pt, float sigIeta, float isoECAL, float HoE, float isoHCAL, float dPhi, float dEta, float isoTracker){
    if (fabs(eta)<1.479){
        if (sigIeta>0.013) return false;
        if (HoE>0.15) return false;
        if (isoECAL/pt>0.5) return false;
        if (isoHCAL/pt>0.3) return false;
        if (fabs(dPhi)>0.07) return false;
        if (fabs(dEta)>0.01) return false;
        if (isoTracker/pt>0.2) return false;
    }
    else{
        if (sigIeta>0.035) return false;
        if (HoE>0.4) return false;
        if (isoECAL/pt>0.5) return false;
        if (isoHCAL/pt>0.4) return false;
        if (fabs(dPhi)>0.1) return false;
        if (fabs(dEta)>0.015) return false;
        if ((isoTracker/pt>0.2)&&(isoTracker!=-1)) return false;
    }
    return true;
    
}


bool miniTreeBuilder::isPassingBit(int i, int triggerInfo){
    int bitVal = pow(2,i);
    int afterMask = bitVal&triggerInfo;
    
    if (afterMask>0) return 1;
    else return 0;
}