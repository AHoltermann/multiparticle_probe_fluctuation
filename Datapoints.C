#include <complex>
#include <tuple>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "Quantities.C"
#include "Sampler.C"

void Datapoints(const char* outfile = "EPL_Line_iter3_raw.root",int isEPL = 1,int n = 20000,int n_samples = 40){

    //Gauss Params 
    double mu_1 = 0.07;
    double sig_1 = 0.55*mu_1;
    double mu_2 = 0.05;
    double sig_2;
    //std::vector<double> sig_2s{0.25*(mu_2/mu_1)*sig_1,0.50*(mu_2/mu_1)*sig_1,0.75*(mu_2/mu_1)*sig_1,1*(mu_2/mu_1)*sig_1,1.25*(mu_2/mu_1)*sig_1,1.50*(mu_2/mu_1)*sig_1,2.00*(mu_2/mu_1)*sig_1};
    std::vector<double> sig_2s{0.80*(mu_2/mu_1)*sig_1,1.00*(mu_2/mu_1)*sig_1,1.2*(mu_2/mu_1)*sig_1,1.4*(mu_2/mu_1)*sig_1,1.6*(mu_2/mu_1)*sig_1,1.80*(mu_2/mu_1)*sig_1,2.00*(mu_2/mu_1)*sig_1};
    
    /*std::vector<double> sig_2s{0.25*(mu_2/mu_1)*sig_1,0.30*(mu_2/mu_1)*sig_1,0.35*(mu_2/mu_1)*sig_1,0.40*(mu_2/mu_1)*sig_1,0.45*(mu_2/mu_1)*sig_1,0.50*(mu_2/mu_1)*sig_1,0.55*(mu_2/mu_1)*sig_1,
0.60*(mu_2/mu_1)*sig_1,0.65*(mu_2/mu_1)*sig_1,0.70*(mu_2/mu_1)*sig_1,0.75*(mu_2/mu_1)*sig_1,0.80*(mu_2/mu_1)*sig_1,0.85*(mu_2/mu_1)*sig_1,0.90*(mu_2/mu_1)*sig_1,
0.95*(mu_2/mu_1)*sig_1,1.00*(mu_2/mu_1)*sig_1,1.05*(mu_2/mu_1)*sig_1,1.10*(mu_2/mu_1)*sig_1,1.15*(mu_2/mu_1)*sig_1,1.20*(mu_2/mu_1)*sig_1,1.25*(mu_2/mu_1)*sig_1,1.30*(mu_2/mu_1)*sig_1,
1.35*(mu_2/mu_1)*sig_1,1.40*(mu_2/mu_1)*sig_1,1.45*(mu_2/mu_1)*sig_1,1.50*(mu_2/mu_1)*sig_1,1.55*(mu_2/mu_1)*sig_1,1.60*(mu_2/mu_1)*sig_1,1.65*(mu_2/mu_1)*sig_1,
1.70*(mu_2/mu_1)*sig_1,1.75*(mu_2/mu_1)*sig_1,1.80*(mu_2/mu_1)*sig_1,1.85*(mu_2/mu_1)*sig_1,1.90*(mu_2/mu_1)*sig_1,1.95*(mu_2/mu_1)*sig_1,1.00*(mu_2/mu_1)*sig_1,
1.05*(mu_2/mu_1)*sig_1,1.10*(mu_2/mu_1)*sig_1,1.15*(mu_2/mu_1)*sig_1,1.20*(mu_2/mu_1)*sig_1,1.25*(mu_2/mu_1)*sig_1,1.30*(mu_2/mu_1)*sig_1,1.35*(mu_2/mu_1)*sig_1,1.40*(mu_2/mu_1)*sig_1,1.45*(mu_2/mu_1)*sig_1,
1.50*(mu_2/mu_1)*sig_1,1.55*(mu_2/mu_1)*sig_1,1.60*(mu_2/mu_1)*sig_1,1.65*(mu_2/mu_1)*sig_1,1.70*(mu_2/mu_1)*sig_1,1.75*(mu_2/mu_1)*sig_1,1.80*(mu_2/mu_1)*sig_1,
1.85*(mu_2/mu_1)*sig_1,1.90*(mu_2/mu_1)*sig_1,1.95*(mu_2/mu_1)*sig_1,2.00*(mu_2/mu_1)*sig_1};
*/   

    //EPD Params
    double e1 = 0.05;
    double a1 = 120;
    double e2 = 0.03;
    double a2;
    std::vector<int> a_2s{2,5,10,30,70,120,500};
    //std::vector<int> a_2s{2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100,150,200,250,300,400,500};
    double rho;
    std::vector<double> rhos{0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.999};
    /*std::vector<double> rhos{0.001,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.999};*/

    cout << "starting" << endl;

    TFile *out = new TFile(outfile,"RECREATE");
    TTree*T = new TTree("T","T");

    double c_mu_v,c_mu_vp,c_sig_v,c_sig_vp,c_RS,c_rho,c_Q_n,c_scalarprod;

    double m_X,m_Y,m_Z;
    double m_XX,m_XY,m_YY,m_XZ;
    double m_XXX,m_XXY,m_XYY,m_XXZ;

    double cm_XX,cm_XY,cm_YY,cm_XZ;
    double cm_XXX,cm_XXY,cm_XYY,cm_XXZ;

    double g_XY_XX,g_YY_XX,g_XZ_XX;
    double g_XXY_XXX,g_XYY_XXX,g_XXZ_XXX;

    double z_XY_XX,z_YY_XX,z_XZ_XX;
    double z_XXY_XXX,z_XYY_XXX,z_XXZ_XXX;
        
    cout << "vars + files init" << endl;

    T->Branch("mu_v",&c_mu_v);
    T->Branch("mu_vp",&c_mu_vp);
    T->Branch("sig_v",&c_sig_v);
    T->Branch("sig_vp",&c_sig_vp);
    T->Branch("RS",&c_RS);
    T->Branch("rho",&c_rho);
    T->Branch("Qn",&c_Q_n);
    T->Branch("scalarprod",&c_scalarprod);

    T->Branch("m_X",&m_X);
    T->Branch("m_Y",&m_Y);
    T->Branch("m_Z",&m_Z);

    T->Branch("m_XX",&m_XX);
    T->Branch("m_XY",&m_XY);
    T->Branch("m_YY",&m_YY);
    T->Branch("m_XZ",&m_XZ);

    T->Branch("m_XXX",&m_XXX);
    T->Branch("m_XXY",&m_XXY);
    T->Branch("m_XYY",&m_XYY);
    T->Branch("m_XXZ",&m_XXZ);

    T->Branch("cm_XX",&cm_XX);
    T->Branch("cm_XY",&cm_XY);
    T->Branch("cm_YY",&cm_YY);
    T->Branch("cm_XZ",&cm_XZ);

    T->Branch("cm_XXX",&cm_XXX);
    T->Branch("cm_XXY",&cm_XXY);
    T->Branch("cm_XYY",&cm_XYY);
    T->Branch("cm_XXZ",&cm_XXZ);

    T->Branch("g_XY_XX",&g_XY_XX);
    T->Branch("g_YY_XX",&g_YY_XX);
    T->Branch("g_XZ_XX",&g_XZ_XX);

    T->Branch("z_XY_XX",&z_XY_XX);
    T->Branch("z_YY_XX",&z_YY_XX);
    T->Branch("z_XZ_XX",&z_XZ_XX);

    T->Branch("g_XXY_XXX",&g_XXY_XXX);
    T->Branch("g_XYY_XXX",&g_XYY_XXX);
    T->Branch("g_XXZ_XXX",&g_XXZ_XXX);

    T->Branch("z_XXY_XXX",&z_XXY_XXX);
    T->Branch("z_XYY_XXX",&z_XYY_XXX);
    T->Branch("z_XXZ_XXX",&z_XXZ_XXX);

    cout << "branch init" << endl;

    for(int l = 0; l<sig_2s.size(); l++){
        for(int j = 0; j<rhos.size(); j++){
            cout << sig_2s.size()*rhos.size() << endl;
            cout << "iter: " << j+l*rhos.size() << endl;
            
            sig_2 = sig_2s[l];
            rho = rhos[j];
            a2 = a_2s[l];

            std::vector<double> x(n*n_samples);
            std::vector<double> y(n*n_samples);

            std::vector<double> x_1(n);
            std::vector<double> y_1(n);
            
            //cout << rho << endl;

            if(isEPL==0){
            Gauss_Sample(x,y,n*n_samples,mu_1,mu_2,sig_1,sig_2,rho);
            }
            else{
            fEPL_Sample(x,y,n*n_samples,e1,e2,a1,a2,rho);  
            }

            for(int i = 0; i< n_samples; i++){
                //cout << "sampling" << endl;
                for(int k = 0; k<n; k++){

                    //if((i*n+k)%100 == 0){
                    //    cout << (i*n+k)/(1.0*n*n_samples) << endl;
                    //}
                    x_1[k] = x[i*n+k];
                    y_1[k] = y[i*n+k];
                }


                    std::vector<double> X = prod(x_1,y_1,2,0);
                    std::vector<double> Y = prod(x_1,y_1,1,1);
                    std::vector<double> Z = prod(x_1,y_1,0,2);
                        
                        c_mu_v = moment(x_1,y_1,x_1,1,0,0);
                        c_mu_vp = moment(x_1,y_1,x_1,0,1,0);
                        c_sig_v = sqrt(variance(x_1,x_1));
                        c_sig_vp = sqrt(variance(y_1,y_1));
                        c_RS = RS(x_1,y_1);
                        c_rho = variance(x_1,y_1)/(c_sig_v*c_sig_vp);
                        //cout << c_rho << endl;
                        c_Q_n = Q_n(x_1,y_1);
                        c_scalarprod = Scalarprod(x_1,y_1);

                        m_X = moment(X,Y,Z,1,0,0);
                        m_Y = moment(X,Y,Z,0,1,0);
                        m_Z = moment(X,Y,Z,0,0,1);

                        m_XX = moment(X,Y,Z,2,0,0);
                        m_XY = moment(X,Y,Z,1,1,0);
                        m_YY = moment(X,Y,Z,0,2,0);
                        m_XZ = moment(X,Y,Z,1,0,1);

                        m_XXX = moment(X,Y,Z,3,0,0);
                        m_XXY = moment(X,Y,Z,2,1,0);
                        m_XYY = moment(X,Y,Z,1,2,0);
                        m_XXZ = moment(X,Y,Z,2,0,1);

                        cm_XX = variance(X,X);
                        cm_XY = variance(X,Y);
                        cm_YY = variance(Y,Y);
                        cm_XZ = variance(X,Z);

                        cm_XXX = skewness(X,X,X);
                        cm_XXY = skewness(X,X,Y);
                        cm_XYY = skewness(X,Y,Y);
                        cm_XXZ = skewness(X,X,Z);

                        g_XY_XX = gamma4(X,Y,X,X);
                        g_YY_XX = gamma4(Y,Y,X,X);
                        g_XZ_XX = gamma4(X,Z,X,X);

                        g_XXY_XXX = gamma6(X,X,Y,X,X,X);
                        g_XYY_XXX = gamma6(X,Y,Y,X,X,X);
                        g_XXZ_XXX = gamma6(X,X,Z,X,X,X);

                        z_XY_XX = gamma4(X,Y,X,X);
                        z_YY_XX = gamma4(Y,Y,X,X);
                        z_XZ_XX = gamma4(X,Z,X,X);

                        z_XXY_XXX = zeta6(X,X,Y,X,X,X);
                        z_XYY_XXX = zeta6(X,Y,Y,X,X,X);
                        z_XXZ_XXX = zeta6(X,X,Z,X,X,X);

                        T->Fill();
                    
                
            }
        }

    }
         
    out->WriteObject(T,"T");
    out->Close();
         

}

