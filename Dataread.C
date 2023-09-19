#include <complex>
#include <tuple>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

void Dataread(const char* infile = "EPL_Line_iter3_raw.root",const char*outfile = "EPL3.root",int iters = 77,int n_samples = 40){

    TFile*f = new TFile(infile);
    TTree *T = (TTree*)f -> Get("T");
    TFile*g = new TFile(outfile,"RECREATE");
    TTree *To = new TTree("T","T");
    
    double c_mu_v,c_mu_vp,c_sig_v,c_sig_vp,c_RS,c_rho,c_Q_n,scalarprod;
    
    double m_X,m_Y,m_Z;
    double m_XX,m_XY,m_YY,m_XZ;
    double m_XXX,m_XXY,m_XYY,m_XXZ;
    
    double cm_XX,cm_XY,cm_YY,cm_XZ;
    double cm_XXX,cm_XXY,cm_XYY,cm_XXZ;

    double g_XY_XX,g_YY_XX,g_XZ_XX;
    double g_XXY_XXX,g_XYY_XXX,g_XXZ_XXX;
    
    double z_XY_XX,z_YY_XX,z_XZ_XX;
    double z_XXY_XXX,z_XYY_XXX,z_XXZ_XXX;

    double dc_mu_v,dc_mu_vp,dc_sig_v,dc_sig_vp,dc_RS,dc_rho,dc_Q_n,dscalarprod;
    
    double dm_X,dm_Y,dm_Z;
    double dm_XX,dm_XY,dm_YY,dm_XZ;
    double dm_XXX,dm_XXY,dm_XYY,dm_XXZ;

    double dcm_XX,dcm_XY,dcm_YY,dcm_XZ;
    double dcm_XXX,dcm_XXY,dcm_XYY,dcm_XXZ;
    
    double dg_XY_XX,dg_YY_XX,dg_XZ_XX;
    double dg_XXY_XXX,dg_XYY_XXX,dg_XXZ_XXX;
    
    double dz_XY_XX,dz_YY_XX,dz_XZ_XX;
    double dz_XXY_XXX,dz_XYY_XXX,dz_XXZ_XXX;

    double ec_mu_v,ec_mu_vp,ec_sig_v,ec_sig_vp,ec_RS,ec_rho,ec_Q_n,escalarprod;
    
    double em_X,em_Y,em_Z;
    double em_XX,em_XY,em_YY,em_XZ;
    double em_XXX,em_XXY,em_XYY,em_XXZ;

    double ecm_XX,ecm_XY,ecm_YY,ecm_XZ;
    double ecm_XXX,ecm_XXY,ecm_XYY,ecm_XXZ;

    double eg_XY_XX,eg_YY_XX,eg_XZ_XX;
    double eg_XXY_XXX,eg_XYY_XXX,eg_XXZ_XXX;
    
    double ez_XY_XX,ez_YY_XX,ez_XZ_XX;
    double ez_XXY_XXX,ez_XYY_XXX,ez_XXZ_XXX;

    TH1D*hc_mu_v = new TH1D("hc_mu_v","hc_mu_v",1000,-1000,1000);
    TH1D*hc_mu_vp = new TH1D("hc_mu_vp","hc_mu_vp",1000,-1000,1000);
    TH1D*hc_sig_v = new TH1D("hc_sig_v","hc_sig_v",1000,-1000,1000);
    TH1D*hc_sig_vp = new TH1D("hc_sig_vp","hc_sig_vp",1000,-1000,1000);
    TH1D*hc_RS = new TH1D("hc_RS","hc_RS",1000,-1000,1000);
    TH1D*hc_rho = new TH1D("hc_rho","hc_rho",1000,-1000,1000);
    TH1D*hc_Q_n = new TH1D("hc_Q_n","hc_Q_n",1000,-1000,1000);
    TH1D*hscalarprod= new TH1D("hc_scalarprod","hc_scalarpord",1000,-1000,1000);
    
    TH1D*hm_X = new TH1D("hm_X","hm_X",1000,-1000,1000);
    TH1D*hm_Y = new TH1D("hm_Y","hm_Y",1000,-1000,1000);
    TH1D*hm_Z = new TH1D("hm_Z","hm_Z",1000,-1000,1000);

    TH1D*hm_XX = new TH1D("hm_XX","hm_XX",1000,-1000,1000);
    TH1D*hm_XY = new TH1D("hm_XY","hm_XY",1000,-1000,1000);
    TH1D*hm_YY = new TH1D("hm_YY","hm_YY",1000,-1000,1000);
    TH1D*hm_XZ = new TH1D("hm_XZ","hm_XZ",1000,-1000,1000);

    TH1D*hm_XXX = new TH1D("hm_XXX","hm_XXX",1000,-1000,1000);
    TH1D*hm_XXY = new TH1D("hm_XXY","hm_XXY",1000,-1000,1000);
    TH1D*hm_XYY = new TH1D("hm_XYY","hm_XYY",1000,-1000,1000);
    TH1D*hm_XXZ = new TH1D("hm_XXZ","hm_XXZ",1000,-1000,1000);

    TH1D*hcm_XX = new TH1D("hcm_XX","hcm_XX",1000,-1000,1000);
    TH1D*hcm_XY = new TH1D("hcm_XY","hcm_XY",1000,-1000,1000);
    TH1D*hcm_YY = new TH1D("hcm_YY","hcm_YY",1000,-1000,1000);
    TH1D*hcm_XZ = new TH1D("hcm_XZ","hcm_XZ",1000,-1000,1000);

    TH1D*hcm_XXX = new TH1D("hcm_XXX","hcm_XXX",1000,-1000,1000);
    TH1D*hcm_XXY = new TH1D("hcm_XXY","hcm_XXY",1000,-1000,1000);
    TH1D*hcm_XYY = new TH1D("hcm_XYY","hcm_XYY",1000,-1000,1000);
    TH1D*hcm_XXZ = new TH1D("hcm_XXZ","hcm_XXZ",1000,-1000,1000);

    TH1D*hg_XY_XX = new TH1D("hg_XY_XX","hg_XY_XX",1000,-1000,1000);
    TH1D*hg_YY_XX = new TH1D("hg_YY_XX","hg_YY_XX",1000,-1000,1000);
    TH1D*hg_XZ_XX = new TH1D("hg_XZ_XX","hg_XZ_XX",1000,-1000,1000);

    TH1D*hg_XXY_XXX = new TH1D("hg_XXY_XXX","hg_XXY_XXX",1000,-1000,1000);
    TH1D*hg_XYY_XXX = new TH1D("hg_XYY_XXX","hg_XYY_XXX",1000,-1000,1000);
    TH1D*hg_XXZ_XXX = new TH1D("hg_XXZ_XXX","hg_XXZ_XXX",1000,-1000,1000);

    TH1D*hz_XY_XX = new TH1D("hz_XY_XX","hz_XY_XX",1000,-1000,1000);
    TH1D*hz_YY_XX = new TH1D("hz_YY_XX","hz_YY_XX",1000,-1000,1000);
    TH1D*hz_XZ_XX = new TH1D("hz_XZ_XX","hz_XZ_XX",1000,-1000,1000);

    TH1D*hz_XXY_XXX = new TH1D("hz_XXY_XXX","hz_XXY_XXX",1000,-1000,1000);
    TH1D*hz_XYY_XXX = new TH1D("hz_XYY_XXX","hz_XYY_XXX",1000,-1000,1000);
    TH1D*hz_XXZ_XXX = new TH1D("hz_XXZ_XXX","hz_XXZ_XXX",1000,-1000,1000);
        
    T->SetBranchAddress("mu_v",&c_mu_v);
    T->SetBranchAddress("mu_vp",&c_mu_vp);
    T->SetBranchAddress("sig_v",&c_sig_v);
    T->SetBranchAddress("sig_vp",&c_sig_vp);
    T->SetBranchAddress("RS",&c_RS);
    T->SetBranchAddress("rho",&c_rho);
    T->SetBranchAddress("Qn",&c_Q_n);
    T->SetBranchAddress("scalarprod",&scalarprod);

    T->SetBranchAddress("m_X",&m_X);
    T->SetBranchAddress("m_Y",&m_Y);
    T->SetBranchAddress("m_Z",&m_Z);

    T->SetBranchAddress("m_XX",&m_XX);
    T->SetBranchAddress("m_XY",&m_XY);
    T->SetBranchAddress("m_YY",&m_YY);
    T->SetBranchAddress("m_XZ",&m_XZ);

    T->SetBranchAddress("m_XXX",&m_XXX);
    T->SetBranchAddress("m_XXY",&m_XXY);
    T->SetBranchAddress("m_XYY",&m_XYY);
    T->SetBranchAddress("m_XXZ",&m_XXZ);

    T->SetBranchAddress("cm_XX",&cm_XX);
    T->SetBranchAddress("cm_XY",&cm_XY);
    T->SetBranchAddress("cm_YY",&cm_YY);
    T->SetBranchAddress("cm_XZ",&cm_XZ);

    T->SetBranchAddress("cm_XXX",&cm_XXX);
    T->SetBranchAddress("cm_XXY",&cm_XXY);
    T->SetBranchAddress("cm_XYY",&cm_XYY);
    T->SetBranchAddress("cm_XXZ",&cm_XXZ);

    T->SetBranchAddress("g_XY_XX",&g_XY_XX);
    T->SetBranchAddress("g_YY_XX",&g_YY_XX);
    T->SetBranchAddress("g_XZ_XX",&g_XZ_XX);

    T->SetBranchAddress("g_XXY_XXX",&g_XXY_XXX);
    T->SetBranchAddress("g_XYY_XXX",&g_XYY_XXX);
    T->SetBranchAddress("g_XXZ_XXX",&g_XXZ_XXX);

    T->SetBranchAddress("z_XY_XX",&z_XY_XX);
    T->SetBranchAddress("z_YY_XX",&z_YY_XX);
    T->SetBranchAddress("z_XZ_XX",&z_XZ_XX);

    T->SetBranchAddress("z_XXY_XXX",&z_XXY_XXX);
    T->SetBranchAddress("z_XYY_XXX",&z_XYY_XXX);
    T->SetBranchAddress("z_XXZ_XXX",&z_XXZ_XXX);

    To->Branch("mu_v",&dc_mu_v);
    To->Branch("mu_vp",&dc_mu_vp);
    To->Branch("sig_v",&dc_sig_v);
    To->Branch("sig_vp",&dc_sig_vp);
    To->Branch("RS",&dc_RS);
    To->Branch("rho",&dc_rho);
    To->Branch("Qn",&dc_Q_n);
    To->Branch("scalarprod",&dscalarprod);

    To->Branch("m_X",&dm_X);
    To->Branch("m_Y",&dm_Y);
    To->Branch("m_Z",&dm_Z);

    To->Branch("m_XX",&dm_XX);
    To->Branch("m_XY",&dm_XY);
    To->Branch("m_YY",&dm_YY);
    To->Branch("m_XZ",&dm_XZ);

    To->Branch("m_XXX",&dm_XXX);
    To->Branch("m_XXY",&dm_XXY);
    To->Branch("m_XYY",&dm_XYY);
    To->Branch("m_XXZ",&dm_XXZ);

    To->Branch("cm_XX",&dcm_XX);
    To->Branch("cm_XY",&dcm_XY);
    To->Branch("cm_YY",&dcm_YY);
    To->Branch("cm_XZ",&dcm_XZ);

    To->Branch("cm_XXX",&dcm_XXX);
    To->Branch("cm_XXY",&dcm_XXY);
    To->Branch("cm_XYY",&dcm_XYY);
    To->Branch("cm_XXZ",&dcm_XXZ);

    To->Branch("g_XY_XX",&dg_XY_XX);
    To->Branch("g_YY_XX",&dg_YY_XX);
    To->Branch("g_XZ_XX",&dg_XZ_XX);

    To->Branch("g_XXY_XXX",&dg_XXY_XXX);
    To->Branch("g_XYY_XXX",&dg_XYY_XXX);
    To->Branch("g_XXZ_XXX",&dg_XXZ_XXX);

    To->Branch("z_XY_XX",&dz_XY_XX);
    To->Branch("z_YY_XX",&dz_YY_XX);
    To->Branch("z_XZ_XX",&dz_XZ_XX);

    To->Branch("z_XXY_XXX",&dz_XXY_XXX);
    To->Branch("z_XYY_XXX",&dz_XYY_XXX);
    To->Branch("z_XXZ_XXX",&dz_XXZ_XXX);

    To->Branch("emu_v",&ec_mu_v);
    To->Branch("emu_vp",&ec_mu_vp);
    To->Branch("esig_v",&ec_sig_v);
    To->Branch("esig_vp",&ec_sig_vp);
    To->Branch("eRS",&ec_RS);
    To->Branch("erho",&ec_rho);
    To->Branch("eQn",&ec_Q_n);
    To->Branch("escalarprod",&escalarprod);

    To->Branch("em_X",&em_X);
    To->Branch("em_Y",&em_Y);
    To->Branch("em_Z",&em_Z);

    To->Branch("em_XX",&em_XX);
    To->Branch("em_XY",&em_XY);
    To->Branch("em_YY",&em_YY);
    To->Branch("em_XZ",&em_XZ);

    To->Branch("em_XXX",&em_XXX);
    To->Branch("em_XXY",&em_XXY);
    To->Branch("em_XYY",&em_XYY);
    To->Branch("em_XXZ",&em_XXZ);

    To->Branch("ecm_XX",&ecm_XX);
    To->Branch("ecm_XY",&ecm_XY);
    To->Branch("ecm_YY",&ecm_YY);
    To->Branch("ecm_XZ",&ecm_XZ);

    To->Branch("ecm_XXX",&ecm_XXX);
    To->Branch("ecm_XXY",&ecm_XXY);
    To->Branch("ecm_XYY",&ecm_XYY);
    To->Branch("ecm_XXZ",&ecm_XXZ);

    To->Branch("eg_XY_XX",&eg_XY_XX);
    To->Branch("eg_YY_XX",&eg_YY_XX);
    To->Branch("eg_XZ_XX",&eg_XZ_XX);

    To->Branch("eg_XXY_XXX",&eg_XXY_XXX);
    To->Branch("eg_XYY_XXX",&eg_XYY_XXX);
    To->Branch("eg_XXZ_XXX",&eg_XXZ_XXX);

    To->Branch("ez_XY_XX",&ez_XY_XX);
    To->Branch("ez_YY_XX",&ez_YY_XX);
    To->Branch("ez_XZ_XX",&ez_XZ_XX);

    To->Branch("ez_XXY_XXX",&ez_XXY_XXX);
    To->Branch("ez_XYY_XXX",&ez_XYY_XXX);
    To->Branch("ez_XXZ_XXX",&ez_XXZ_XXX);

    
    for( int i = 0; i< iters; i++){
        for(int j = 0; j< n_samples; j++){
            T->GetEntry(i*n_samples+j);
            
            hc_mu_v->Fill(c_mu_v);
            hc_mu_vp->Fill(c_mu_vp);
            hc_sig_v->Fill(c_sig_v);
            hc_sig_vp->Fill(c_sig_vp);
            hc_RS->Fill(c_RS);
            hc_rho->Fill(c_rho);
            hc_Q_n->Fill(c_Q_n);
            hscalarprod->Fill(scalarprod);

            hm_X->Fill(m_X);
            hm_Y->Fill(m_Y);
            hm_Z->Fill(m_Z);

            hm_XX->Fill(m_XX);
            hm_XY->Fill(m_XY);
            hm_YY->Fill(m_YY);
            hm_XZ->Fill(m_XZ);

            hm_XXX->Fill(m_XXX);
            hm_XXY->Fill(m_XXY);
            hm_XYY->Fill(m_XYY);
            hm_XXZ->Fill(m_XXZ);

            hcm_XX->Fill(cm_XX);
            hcm_XY->Fill(cm_XY);
            hcm_YY->Fill(cm_YY);
            hcm_XZ->Fill(cm_XZ);

            hcm_XXX->Fill(cm_XXX);
            hcm_XXY->Fill(cm_XXY);
            hcm_XYY->Fill(cm_XYY);
            hcm_XXZ->Fill(cm_XXZ);

            hg_XY_XX->Fill(g_XY_XX);
            hg_YY_XX->Fill(g_YY_XX);
            hg_XZ_XX->Fill(g_XZ_XX);

            hg_XXY_XXX->Fill(g_XXY_XXX);
            hg_XYY_XXX->Fill(g_XYY_XXX);
            hg_XXZ_XXX->Fill(g_XXZ_XXX);

            hz_XY_XX->Fill(z_XY_XX);
            hz_YY_XX->Fill(z_YY_XX);
            hz_XZ_XX->Fill(z_XZ_XX);

            hz_XXY_XXX->Fill(z_XXY_XXX);
            hz_XYY_XXX->Fill(z_XYY_XXX);
            hz_XXZ_XXX->Fill(z_XXZ_XXX);

        }

            dc_mu_v = hc_mu_v->GetMean();
           dc_mu_vp = hc_mu_vp->GetMean();
           dc_sig_v = hc_sig_v->GetMean();
           dc_sig_vp = hc_sig_vp->GetMean();
           dc_RS = hc_RS->GetMean();
           dc_rho = hc_rho->GetMean();
           dc_Q_n = hc_Q_n->GetMean();
           dscalarprod = hscalarprod->GetMean();

           dm_X = hm_X->GetMean();
           dm_Y = hm_Y->GetMean();
           dm_Z = hm_Z->GetMean();

           dm_XX = hm_XX->GetMean();
           dm_XY = hm_XY->GetMean();
           dm_YY = hm_YY->GetMean();
           dm_XZ = hm_XZ->GetMean();

           dm_XXX = hm_XXX->GetMean();
           dm_XXY = hm_XXY->GetMean();
           dm_XYY = hm_XYY->GetMean();
           dm_XXZ = hm_XXZ->GetMean();

           dcm_XX = hcm_XX->GetMean();
           dcm_XY = hcm_XY->GetMean();
           dcm_YY = hcm_YY->GetMean();
           dcm_XZ = hcm_XZ->GetMean();

           dcm_XXX = hcm_XXX->GetMean();
           dcm_XXY = hcm_XXY->GetMean();
           dcm_XYY = hcm_XYY->GetMean();
           dcm_XXZ = hcm_XXZ->GetMean();

           dg_XY_XX = hg_XY_XX->GetMean();
           dg_YY_XX = hg_YY_XX->GetMean();
           dg_XZ_XX = hg_XZ_XX->GetMean();

           dg_XXY_XXX = hg_XXY_XXX->GetMean();
           dg_XYY_XXX = hg_XYY_XXX->GetMean();
           dg_XXZ_XXX = hg_XXZ_XXX->GetMean();

           dz_XY_XX = hz_XY_XX->GetMean();
           dz_YY_XX = hz_YY_XX->GetMean();
           dz_XZ_XX = hz_XZ_XX->GetMean();

           dz_XXY_XXX = hz_XXY_XXX->GetMean();
           dz_XYY_XXX = hz_XYY_XXX->GetMean();
           dz_XXZ_XXX = hz_XXZ_XXX->GetMean();

           ec_mu_v = hc_mu_v->GetStdDev()/sqrt(n_samples);
           ec_mu_vp = hc_mu_vp->GetStdDev()/sqrt(n_samples);
           ec_sig_v = hc_sig_v->GetStdDev()/sqrt(n_samples);
           ec_sig_vp = hc_sig_vp->GetStdDev()/sqrt(n_samples);
           ec_RS = hc_RS->GetStdDev()/sqrt(n_samples);
           ec_rho = hc_rho->GetStdDev()/sqrt(n_samples);
           ec_Q_n = hc_Q_n->GetStdDev()/sqrt(n_samples);
           escalarprod = hscalarprod->GetStdDev()/sqrt(n_samples);

           em_X = hm_X->GetStdDev()/sqrt(n_samples);
           em_Y = hm_Y->GetStdDev()/sqrt(n_samples);
           em_Z = hm_Z->GetStdDev()/sqrt(n_samples);

           em_XX = hm_XX->GetStdDev()/sqrt(n_samples);
           em_XY = hm_XY->GetStdDev()/sqrt(n_samples);
           em_YY = hm_YY->GetStdDev()/sqrt(n_samples);
           em_XZ = hm_XZ->GetStdDev()/sqrt(n_samples);

           em_XXX = hm_XXX->GetStdDev()/sqrt(n_samples);
           em_XXY = hm_XXY->GetStdDev()/sqrt(n_samples);
           em_XYY = hm_XYY->GetStdDev()/sqrt(n_samples);
           em_XXZ = hm_XXZ->GetStdDev()/sqrt(n_samples);

           ecm_XX = hcm_XX->GetStdDev()/sqrt(n_samples);
           ecm_XY = hcm_XY->GetStdDev()/sqrt(n_samples);
           ecm_YY = hcm_YY->GetStdDev()/sqrt(n_samples);
           ecm_XZ = hcm_XZ->GetStdDev()/sqrt(n_samples);

           ecm_XXX = hcm_XXX->GetStdDev()/sqrt(n_samples);
           ecm_XXY = hcm_XXY->GetStdDev()/sqrt(n_samples);
           ecm_XYY = hcm_XYY->GetStdDev()/sqrt(n_samples);
           ecm_XXZ = hcm_XXZ->GetStdDev()/sqrt(n_samples);

           eg_XY_XX = hg_XY_XX->GetStdDev()/sqrt(n_samples);
           eg_YY_XX = hg_YY_XX->GetStdDev()/sqrt(n_samples);
           eg_XZ_XX = hg_XZ_XX->GetStdDev()/sqrt(n_samples);

           eg_XXY_XXX = hg_XXY_XXX->GetStdDev()/sqrt(n_samples);
           eg_XYY_XXX = hg_XYY_XXX->GetStdDev()/sqrt(n_samples);
           eg_XXZ_XXX = hg_XXZ_XXX->GetStdDev()/sqrt(n_samples);

           ez_XY_XX = hz_XY_XX->GetStdDev()/sqrt(n_samples);
           ez_YY_XX = hz_YY_XX->GetStdDev()/sqrt(n_samples);
           ez_XZ_XX = hz_XZ_XX->GetStdDev()/sqrt(n_samples);

           ez_XXY_XXX = hz_XXY_XXX->GetStdDev()/sqrt(n_samples);
           ez_XYY_XXX = hz_XYY_XXX->GetStdDev()/sqrt(n_samples);
           ez_XXZ_XXX = hz_XXZ_XXX->GetStdDev()/sqrt(n_samples);

            hc_mu_v->Reset();
            hc_mu_vp->Reset();
            hc_sig_v->Reset();
            hc_sig_vp->Reset();
            hc_RS->Reset();
            hc_rho->Reset();
            hc_Q_n->Reset();
            hscalarprod->Reset();

            hm_XXX->Reset();
            hm_XXY->Reset();
            hm_XYY->Reset();
            hm_XXZ->Reset();

            hcm_XXX->Reset();
            hcm_XXY->Reset();
            hcm_XYY->Reset();
            hcm_XXZ->Reset();

            hm_XX->Reset();
            hm_XY->Reset();
            hm_YY->Reset();
            hm_XZ->Reset();

            hcm_XX->Reset();
            hcm_XY->Reset();
            hcm_YY->Reset();
            hcm_XZ->Reset();

            hm_X->Reset();
            hm_Y->Reset();
            hm_Z->Reset();

            hg_XY_XX->Reset();
            hg_YY_XX->Reset();
            hg_XZ_XX->Reset();

            hg_XXY_XXX->Reset();
            hg_XYY_XXX->Reset();
            hg_XXZ_XXX->Reset();

            hz_XY_XX->Reset();
            hz_YY_XX->Reset();
            hz_XZ_XX->Reset();

            hz_XXY_XXX->Reset();
            hz_XYY_XXX->Reset();
            hz_XXZ_XXX->Reset();


        To->Fill();


    }

    g->cd();
    To->Write();
    g->Close();       


}