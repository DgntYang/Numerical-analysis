#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>

using namespace std;
class Romberg{
public:
    double fun(double x){
        //return 1.0/(1+100*pow(x, 2));
        //return sin(x)/x;
        return 1.0/x;
    }

    double Tn(double a, double b, int n){
        double h = (b-a) / double(n);
        double x1, x2, t = 0;
        for (int i=0; i<n; i++){
            x1 = a + i*h;
            x2 = a + (i+1)*h;
            t += 0.5*h*(fun(x1)+fun(x2));
        }
        return t;
    }

    double Sn(double t2n, double tn){
        return 4.0/3*t2n - 1.0/3*tn;
    }

    double Cn(double s2n, double sn){
        return 16.0/15*s2n - 1.0/15*sn;
    }

    double Rn(double c2n, double cn){
        return 64.0/63*c2n - 1.0/63*cn;
    }

};

int main(){
    int n=1, a=2, b=8;
    vector<double> Tn, Sn, Cn, Rn;
    double tn, sn, cn, rn;
    double diff=100, epsilon=0.5*pow(10,-5);

    tn = Romberg().Tn(a, b, pow(2, 0));
    Tn.push_back(tn);

    while (diff>epsilon){
        n++;
        tn = Romberg().Tn(a, b, pow(2, n-1));
        Tn.push_back(tn);

        if(n>=2){
            sn = Romberg().Sn(Tn.at(Tn.size()-1), Tn.at(Tn.size()-2));
            Sn.push_back(sn);
        }

        if (n>=3){
            cn = Romberg().Cn(Sn.at(Sn.size()-1), Sn.at(Sn.size()-2));
            Cn.push_back(cn);
        }

        if (n>=4){
            rn = Romberg().Rn(Cn.at(Cn.size()-1), Cn.at(Cn.size()-2));
            Rn.push_back(cn);
            cout << "Epoch " << n-3 << ",   R_2n = " << setprecision(8) << rn << endl;
        }

        if (Rn.size()>=2){
            double r1=Rn.at(Rn.size()-1), r2=Rn.at(Rn.size()-2);
            diff = abs(r1 - r2)/255;
        }
    }
    for (auto num:Cn){
        cout << num << endl;
    }


    return 0;
}