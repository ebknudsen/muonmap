/*given a point in the sample and an initial spin
 * pick a random number t*/
/*then by precession around the local field
 * muon will decay an emit a positron along
 * preferentially along the spin direction*/
/*trace the positron to the detector pixel, or rather store the ntuple xyzt in either BFLR*/

/*positron flight time?*/
/*positron absorption?*/

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TTree.h"

using namespace ROOT::Math;

#ifndef __CINT__
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/Rotation3D.h"
#include "Math/EulerAngles.h"
#include "Math/AxisAngle.h"
#include "Math/Quaternion.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/RotationZYX.h"
#include "Math/LorentzRotation.h"
#include "Math/Boost.h"
#include "Math/BoostX.h"
#include "Math/BoostY.h"
#include "Math/BoostZ.h"
#include "Math/Transform3D.h"
#include "Math/Plane3D.h"
#include "Math/VectorUtil.h"
#endif

class Muon{
  public:
    Double_t t,p; /*time and weight*/
    XYZVector S;/*spin*/
    XYZPoint r;/*position*/
    Double_t E;
    const Double_t gamma=135.5e6; /*Gyromagnetic ratio Hz/T*/
    const Double_t tau=2.197034e-6; /*mean lifetime*/

    Muon(Double_t x,Double_t y, Double_t z){
      r.SetCoordinates(x,y,z);
      S.SetCoordinates(0,0,0);
    }
    
    void precess(Double_t dt, XYZVector B){
      Double_t omega = 2*TMath::Pi() * gamma*t*B.R();
      AxisAngle rot(Bf.unit(),omega);
      this.S=rot*this.S;
    }

    XYZVector decay_to_positron() {
      /*set up a mersenne twister*/ 
      TRandom3 *generator = new TRandom3(); 

      /*Now add the physical asymmetry*/
      Double_t a=0.3333333333333333333333333333333;
      TF1 *f1 = new TF1("f1","1+0.3333333333333333333333333333*cos(x)",-TMath::Pi(),TMath::Pi());
    }
};
  

void muon_sim(int ncount=1000000, double Bx, double By, double Bz) {
  const Double_t gamma_mu=135.5e6; /*Hz/T gyromagnetic ratio of muon*/ 
  const Double_t tau_mu=2.197034e-6; /*mean lifetime*/

  gSystem->Load("libMathCore");
  gSystem->Load("libGenVector");
  
  TTree *back = new TTree("Back","Back muon detector");
  TTree *front = new TTree("Front","Front muon detector");
  TTree *left = new TTree("Right","Right muon detector");
  TTree *right = new TTree("Left","Left muon detector");

  TCanvas *c1 = new TCanvas("c1","Detected positrons",200,10,600,400);

  
  XYZVector Smu(0.0, 0.0, -1.0);
  
  /*create a tree to hold the data*/
  Double_t B[4],F[4],R[4],L[4], *data;

  back->Branch("b",B,"x/D:y/D:z/D:t/D:p/D"); 
  front->Branch("f",F,"x/D:y/D:z/D:t/D:p/D"); 
  left->Branch("l",L,"x/D:y/D:z/D:t/D:p/D"); 
  right->Branch("r",R,"x/D:y/D:z/D:t/D:p/D"); 

  /*set up a mersenne twister*/ 
  TRandom3 *generator = new TRandom3(); 
  
  /*Now add the physical asymmetry*/
  Double_t a=0.3333333333333333333333333333333;

  TF1 *f1 = new TF1("f1","1+0.3333333333333333333333333333*cos(x)",-TMath::Pi(),TMath::Pi());
    
  XYZVector Bf(Bx,By,Bz);
  XYZPoint P(0.0,0.0,0.0); 

  for ( Int_t n=0; n<ncount; n++ ){
    XYZVector Smup;
    Double_t t=generator->Exp(tau_mu);//-lambda * log( rand() * lambda );
    Double_t omega = 2*TMath::Pi() * gamma_mu*t*Bf.R();
    //t=omega/(gamma_mu*2*Bf.R());

    Double_t theta = f1->GetRandom();
    AxisAngle rot(Bf.unit(),omega);
    Smup=rot*Smu;
    if (n && n%(ncount/10) ==0){
      cout << n/(ncount/100) << "\%" <<endl;
    }

    AxisAngle rot_symmetry_violation(Bf.unit(),theta);
    Double_t eta = generator->Uniform(2*TMath::Pi());
    AxisAngle rot_symmetry_violation2(Smup,eta);

    //Smup=rot_symmetry_violation*(rot_symmetry_violation2*Smup);
    Smup=rot_symmetry_violation*Smup;

    TBranch *branchptr;
    Double_t dz;
    Int_t rc,lc,fc,bc;

    if( fabs(Smup.x())<fabs(Smup.z())){
      /*hit back or front*/
      if ( Smup.z()>0){
        dz=(1-P.z())/Smup.unit().z();
        treeptr=front;
        data=F;
        fc++;
      }else{
        dz=(-1-P.z())/Smup.unit().z();
        treeptr=back;
        data=B;
        bc++;
      }
      data[3] = t;
    }else {
      /*we hit left or right first*/
      if (Smup.x()>0){
        dz=(1-P.x())/Smup.unit().x();
        treeptr=left;
        data=L;
        lc++;
      }else{
        dz=(-1-P.x())/Smup.unit().x();
        treeptr=right;
        data=R;
        rc++;
      }
      data[3] = t;
    }
    data[0] =P.X()+ dz*Smup.unit().X();
    data[1] =P.Y()+ dz*Smup.unit().Y();
    data[2] =P.Z()+ dz*Smup.unit().Z();
    data[4]=1;
    //cout << data[0] << " " << data[1] <<" " << data[2] << " " << data[3] << endl;
    //cout << B[0] << " " << B[1] <<" " << B[2] << " " << B[3] << endl;
    treeptr->Fill();
    //cout << treeptr->GetEntries() << endl;
  }
  cout << "done!\n" <<endl;
  printf("%d %d %d %d %d\n",bc,lc,fc,rc,bc+fc+lc+rc); 
  /*plot the detector histograms*/


  TH1D *hb=new TH1D("hb","Back  Detector I vs. t",500,0,24e-6);
  hb->SetLineWidth(2);
  TH1D *hr=new TH1D("hr","Right Detector I vs. t",500,0,24e-6);
  hr->SetLineColor(2);
  hr->SetLineWidth(2);
  TH1D *hf=new TH1D("hf","Front Detector I vs. t",500,0,24e-6);
  hf->SetLineColor(3);
  hf->SetLineWidth(2);
  TH1D *hl=new TH1D("hl","Left  Detector I vs. t",500,0,24e-6);
  hl->SetLineColor(6);
  hl->SetLineWidth(2);

  /*plot times to histograms*/
  back->Draw("t>>hb");
  right->Draw("t>>hr");
  front->Draw("t>>hf");
  left->Draw("t>>hl");

  /*now plot hists to canvas*/
  hb->Draw();
  //hr->Draw("same");
  hf->Draw("same");
  //hl->Draw("same");
}
