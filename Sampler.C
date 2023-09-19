#include <complex>
#include <tuple>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>


Double_t Biv_Gauss(Double_t *val, Double_t *par){

    double pi = acos(-1.0);
               
    double x = val[0];
    double y = val[1];
    double u1 = par[0];
    double u2 = par[1];
    double s1 = par[2];
    double s2 = par[3];
    double p = par[4];

    double a = 1/(2*pi*s1*s2*sqrt(1-p*p));
   
    double b = -1/(2*(1-p*p));

    double c = pow((val[0] - par[0]),2.0)/(s1*s1) + pow((val[1] - par[1]),2.0)/(s2*s2) - 2*p*(val[0] - par[0])*(val[1] - par[1])/(s1*s2);

    return a*exp(b*c);
}

Double_t EPL_Gauss(Double_t *val, Double_t *par){

    double x = val[0];
    double y = val[1];
    //double x = 0.0697043;
    //double y = 0.258077;

    double a1 = par[0];
    double a2 = par[1];
    double p = par[2];
    
    double f1 = 2*a1*x*pow((1-x*x),a1-1);
    double f2 = 2*a2*y*pow((1-y*y),a2-1);

    double phi1 = ROOT::Math::gaussian_quantile((1-pow((1-x*x),a1)),1);
    double phi2 = ROOT::Math::gaussian_quantile((1-pow((1-y*y),a2)),1);
   
    double z = -0.5*(1/(1 - p*p) - 1)*(phi1*phi1 + phi2*phi2) + p/(1-p*p)*(phi1*phi2);

    double r = (1/sqrt(1-p*p))*exp(z)*f1*f2;
        if(r != r){
            r =0;
        }
    return r;
}

Double_t fEPL(Double_t *val, Double_t *par){

    double x = val[0];

    double e = par[0];
    double a = par[1];

    double t1 = 2*x*a*pow((1-x*x),a-1);
    double t2 = pow((1-e*e),a+0.5);
    double exp = -1-2*a;
    double t3 = pow((1-e*x),exp);
    double t4 = ROOT::Math::hyperg(0.5,1+2*a,1,(2*x*e)/(e*x-1));

    return t1*t2*t3*t4;
}

Double_t fEPL2(Double_t *val, Double_t *par){
    double pi = 3.14159265358979323846;
    double x = val[0];
    double e = par[0];
    double a = par[1];

    double t1 = 2/pi*x*a*pow((1-(x*x)),a-1);
    double t2 = pow((1-e*e),a+0.5);
    double phi = 0;
    double vem = 0;
    
    double t4 = 0;
    for(int i = 1; i<101; i++){
        phi =((2*i*pi-0.5)/200);
        vem = (pi/100)*pow((1-e*x*cos(phi)),-2*a-1);
        t4 = t4+ vem;

    }

    return t1*t2*t4;
        
}

Double_t fEPL3(Double_t *val, Double_t *par){
    double pi = 3.14159265358979323846;
    double x = val[0];
    double e = par[0];
    double a = par[1];
    //double s = par[2];

    double t1 = 2/pi*x*a*pow((1-(x*x)),a-1);
    double t2 = pow((1-e*e),a+0.5);
    double phi = 0;
    
    double t4 = 0;
    TF1*f = new TF1("f","pow((1-[0]*[1]*cos(x)),-2*[2]-1)");
    f->SetParameter(0,e);
    f->SetParameter(1,val[0]);
    f->SetParameter(2,a);
    t4 = f->Integral(0.000001,pi);

    return t2*t1*t4;
      
}


Double_t EPL_Hyper(Double_t *val, Double_t *par){

    double x = val[0];
    double y = val[1];
    //double x = 0.0697043;
    //double y = 0.258077;

    if((x<0)||(y<0)){
        return 0;
    }

    double e1 = par[0];
    double e2 = par[1];

    double a1 = par[2];
    double a2 = par[3];
    double p = par[4];
    

    auto f_1 = new TF1("f",fEPL,0.0,1.0,2);
    f_1->SetParameter(0,e1);
    f_1->SetParameter(1,a1);

    auto f_2 = new TF1("g",fEPL,0.0,1.0,2);
    f_2->SetParameter(0,e2);
    f_2->SetParameter(1,a2);
    


    double f1 = f_1->Eval(x);
    double f2 = f_2->Eval(y);
    double C1 = f_1->Integral(0.0,x);
    double C2 = f_2->Integral(0.0,y);
    

    double phi1 = ROOT::Math::gaussian_quantile(C1,1);
    double phi2 = ROOT::Math::gaussian_quantile(C2,1);
   
    double z = -0.5*(1/(1 - p*p) - 1)*(phi1*phi1 + phi2*phi2) + p/(1-p*p)*(phi1*phi2);

    double r = (1/sqrt(1-p*p))*exp(z)*f1*f2;
        if(r != r){
            r =0;
        }
    return r;

    cout<< "iter" <<endl;



    return 0;
}

void Gauss_Sample(std::vector<double>&x, std::vector<double>&y,double n,double mu1, double mu2, double s1, double s2, double p){


    TRandom *rand = new TRandom2();
    rand->SetSeed();

    //std::vector<double> nvec;
    //std::vector<double> mvec;



    TF2*Gauss = new TF2("Biv",Biv_Gauss,-0.5,1,-0.5,1,5);

    Gauss->SetNpx(300);
    Gauss->SetNpy(300);
    Gauss->SetParameter(0,mu1);
    Gauss->SetParameter(1,mu2);
    Gauss->SetParameter(2,s1);
    Gauss->SetParameter(3,s2);
    Gauss->SetParameter(4,p);

   
    for(int i = 0; i<n; i++){
     double_t a,b;

        Gauss->GetRandom2(a,b,rand);
        
        x[i] = a;
        y[i] = b;
        
    }

}

void EPL_Sample(std::vector<double>&x, std::vector<double>&y,double n,double a1, double a2, double p){


    TRandom *rand = new TRandom2();
    rand->SetSeed();

    //std::vector<double> nvec;
    //std::vector<double> mvec;



    auto EPL = new TF2("EPL", EPL_Gauss,0,1.0,0,1.0,3);

    EPL->SetNpx(300);
    EPL->SetNpy(300);
    EPL->SetParameter(0,a1);
    EPL->SetParameter(1,a2);
    EPL->SetParameter(2,p);

    double_t a,b;
    for(int i = 0; i<n; i++){
    
        EPL->GetRandom2(a,b,rand);
        
        x[i] = a;
        y[i] = b;
        
    }

}

void fEPL_Sample(std::vector<double>&x, std::vector<double>&y,double n,double e1, double e2,double a1, double a2, double p){
    TRandom *rand = new TRandom2();
    rand->SetSeed();

    //std::vector<double> nvec;
    //std::vector<double> mvec;



    TF2*fEPL = new TF2("EPL", EPL_Hyper,0,0.99999,0,0.999999,5);

    fEPL->SetNpx(150);
    fEPL->SetNpy(150);
    fEPL->SetParameter(0,e1);
    fEPL->SetParameter(1,e2);
    fEPL->SetParameter(2,a1);
    fEPL->SetParameter(3,a2);
    fEPL->SetParameter(4,p);

    double_t a,b;
    for(int i = 0; i<n; i++){
    
        fEPL->GetRandom2(a,b,rand);
        
        
        x[i] = a;
        y[i] = b;
        
    }
    cout << "subiter" <<endl;

}