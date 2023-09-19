#include <complex>
#include <tuple>
#include <sstream>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <cmath>
#include <vector>


//Where most mathematically calculated quantities reside. Note a few cross check quantities such as Scalar Product

std::vector<double> prod(const std::vector<double>& v, const std::vector<double>& vp, int nv, int nvp) {
    std::vector<double> out;
    for (size_t i = 0; i < v.size(); ++i) {
        double result = std::pow(v[i], nv) * std::pow(vp[i], nvp);
        out.push_back(result);
    }

    return out;
}


double moment(std::vector<double>X1,std::vector<double>X2,std::vector<double>X3,int n_1, int n_2, int n_3){
    
    int L = X1.size();
    std::vector<double>varray(L);

    for(int i = 0; i< L; i++){
        varray[i] = pow(X1[i],n_1)*pow(X2[i],n_2)*pow(X3[i],n_3);
    }
      
    return std::accumulate(varray.begin(), varray.end(),0.0) / L;

}

//Central moments

double variance(std::vector<double>X1,std::vector<double>X2){

    double m1 = moment(X1,X2,X1,1,0,0);
    double m2 = moment(X1,X2,X1,0,1,0);
    double m12 = moment(X1,X2,X1,1,1,0);

    return m12 - m1*m2;

}

double skewness(std::vector<double>X1,std::vector<double>X2,std::vector<double>X3){

    double m1 = moment(X1,X2,X3,1,0,0);
    double m2 = moment(X1,X2,X3,0,1,0);
    double m3 = moment(X1,X2,X3,0,0,1);

    double m12 = moment(X1,X2,X3,1,1,0);
    double m13 = moment(X1,X2,X3,1,0,1);
    double m23 = moment(X1,X2,X3,0,1,1);

    double m123 = moment(X1,X2,X3,1,1,1);

    return m123 - m12*m3 - m13*m2 - m23*m1 + 2*m1*m2*m3;

}

//Normalized Quantities

double Nvariance(std::vector<double>X1,std::vector<double>X2){

    double v12 = variance(X1,X2);
    double m1 = moment(X1,X2,X1,1,0,0);
    double m2 = moment(X1,X2,X1,0,1,0);

    
    return v12/(m1*m2);

}

double Nskewness(std::vector<double>X1,std::vector<double>X2,std::vector<double>X3){

    double s123 = skewness(X1,X2,X3);
    double m1 = moment(X1,X2,X3,1,0,0);
    double m2 = moment(X1,X2,X3,0,1,0);
    double m3 = moment(X1,X2,X3,0,0,1);
    
    return s123/(m1*m2*m3);

}

//Gamma + Zeta

double gamma4(std::vector<double>X1,std::vector<double>X2,std::vector<double>X3,std::vector<double>X4){
    
    double m12 = moment(X1,X2,X1,1,1,0);
    double m1 = moment(X1,X2,X1,1,0,0);
    double m2 = moment(X1,X2,X1,0,1,0);

    double m34 = moment(X3,X4,X3,1,1,0);
    double m3 = moment(X3,X2,X1,1,0,0);
    double m4 = moment(X4,X2,X1,1,0,0);

    return m12/(m1*m2) - m34/(m3*m4);

}

double gamma6(std::vector<double>X1,std::vector<double>X2,std::vector<double>X3,std::vector<double>X4,std::vector<double>X5,std::vector<double>X6){
    
    double m123 = moment(X1,X2,X1,1,1,1);
    double m1 = moment(X1,X2,X3,1,0,0);
    double m2 = moment(X1,X2,X3,0,1,0);
    double m3 = moment(X1,X2,X3,0,0,1);

    double m456 = moment(X4,X5,X6,1,1,1);
    double m4 = moment(X4,X5,X6,1,0,0);
    double m5 = moment(X4,X5,X6,0,1,0);
    double m6 = moment(X4,X5,X6,0,0,1);

    return m123/(m1*m2*m3) - m456/(m4*m5*m6);

    }

double zeta4(std::vector<double>X1,std::vector<double>X2,std::vector<double>X3,std::vector<double>X4){
    
    double m12 = moment(X1,X2,X1,1,1,0);
    double m1 = moment(X1,X2,X1,1,0,0);
    double m2 = moment(X1,X2,X1,0,1,0);

    double m34 = moment(X3,X4,X3,1,1,0);
    double m3 = moment(X3,X2,X1,1,0,0);
    double m4 = moment(X4,X2,X1,1,0,0);

    return (m12/(m1*m2))/(m34/(m3*m4));


}

double zeta6(std::vector<double>X1,std::vector<double>X2,std::vector<double>X3,std::vector<double>X4,std::vector<double>X5,std::vector<double>X6){
    

    double m123 = moment(X1,X2,X1,1,1,1);
    double m1 = moment(X1,X2,X3,1,0,0);
    double m2 = moment(X1,X2,X3,0,1,0);
    double m3 = moment(X1,X2,X3,0,0,1);

    double m456 = moment(X4,X5,X6,1,1,1);
    double m4 = moment(X4,X5,X6,1,0,0);
    double m5 = moment(X4,X5,X6,0,1,0);
    double m6 = moment(X4,X5,X6,0,0,1);

    return (m123/(m1*m2*m3))/(m456/(m4*m5*m6));
    
}

// Other Quantities of interest

double Scalarprod(std::vector<double>v,std::vector<double>vp){

    double num = moment(v,vp,v,1,1,0);
    double denom = sqrt(moment(v,vp,v,2,0,0));
    return num/denom;
}

double Q_n(std::vector<double>v,std::vector<double>vp){
    
    double num = moment(v,vp,v,1,1,0);
    double denom1 = moment(v,vp,v,2,0,0);
    double denom2 = moment(v,vp,v,0,2,0);
    

    return num/(sqrt(denom1 * denom2));

}

double RS(std::vector<double>v,std::vector<double>vp){
    double vvar = variance(v,v);
    double vm = moment(v,v,v,1,0,0);
    double vpvar = variance(vp,vp);
    double vpm = moment(vp,vp,vp,1,0,0);

    return (sqrt(vpvar)/vpm)/(sqrt(vvar)/vm);
}

double Rdiff(std::vector<double>v,std::vector<double>vp){
    double vvar = variance(v,v);
    double vm = moment(v,v,v,1,0,0);
    double vpvar = variance(vp,vp);
    double vpm = moment(vp,vp,vp,1,0,0);

    return (sqrt(vpvar)/vpm)-(sqrt(vvar)/vm);
}