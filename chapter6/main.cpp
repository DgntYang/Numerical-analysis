//C++
//Created by yangtan on 2022/12/7

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <windows.h>


using namespace std;
class Solution{
public:
    void print_all(vector<double> V1, vector<double> V2,vector<double> V3,vector<double> V4,vector<double> V5,
                   double x0, double h){
        cout << "                      RK4                   AB4                  AB4-AM4             "
                "Improved-AB4-AM4             精确解" << endl;
        for (int i=0; i<V1.size(); i++){
            if (i==0){
                cout << "x = " << setprecision(2) << x0+i*h
                     << ";  \t  y = "<<setprecision(10) << V1[i] << " "
                     << ";  \t\t\t\t  y = "<<setprecision(10) << V2[i] << " "
                     << ";  \t\t\t\t  y = "<<setprecision(10) << V3[i] << " "
                     << ";  \t\t\t\t  y = "<<setprecision(10) << V4[i] << " "
                     << ";  \t\t\t\t  y = "<<setprecision(10) << V5[i] << " " << endl;
            }
            else{
                cout << "x = " << setprecision(2) << x0+i*h
                     << ";  \t  y = "<<setprecision(10) << V1[i] << " "
                     << ";  \t  y = "<<setprecision(10) << V2[i] << " "
                     << ";  \t  y = "<<setprecision(10) << V3[i] << " "
                     << ";  \t  y = "<<setprecision(10) << V4[i] << " "
                     << ";  \t  y = "<<setprecision(10) << V5[i] << " " << endl;
            }
        }
    }

    void print_vector(vector<double> V, double x0, double h){
        for (int i=0; i<V.size(); i++){
            cout << "x = " << setprecision(2) << x0+i*h
                 << ";  \t  y = "<<setprecision(10) << V[i] << "  " << endl;
        }
    }

    double fun(double x, double u){
        //return 2.0/x*u + pow(x,2)*exp(x);
        return -pow(x,2)*pow(u,2);
    }

    double f(double x){
        return 3.0/(1+pow(x,3));
    }

    vector<double> rk4(double x0, double u0, double h, double xup){
        vector<double> u={u0};
        double k1,k2,k3,k4;
        while(x0<=xup){
            k1 = fun(x0, u0);
            k2 = fun(x0+0.5*h, u0+0.5*h*k1);
            k3 = fun(x0+0.5*h, u0+0.5*h*k2);
            k4 = fun(x0+h, u0+h*k3);
            u0 = u0 + h/6*(k1+2*k2+2*k3+k4);
            x0 += h;
            u.push_back(u0);
        }
        return u;
    }

    vector<double> ab4(vector<double> X, vector<double> U, double h, double xup){
        vector<double> u={U[0], U[1], U[2], U[3]};
        double tu3;
        while (X[3]<=xup){
            tu3 = U[3] + h/24*(55*fun(X[3],U[3]) - 59*fun(X[2],U[2]) + 37*fun(X[1],U[1]) - 9*fun(X[0],U[0]));
            u.push_back(tu3);
            X[0]+=h; X[1]+=h; X[2]+=h; X[3]+=h;
            U[0]=U[1]; U[1]=U[2]; U[2]=U[3]; U[3]=tu3;
        }
        return u;
    }

    vector<double> ab4_am4(vector<double> X, vector<double> U, double h, double xup){
        vector<double> u={U[0], U[1], U[2], U[3]};
        double y_p;
        while (X[3]<=xup){
            y_p = U[3] + h/24*(55*fun(X[3],U[3]) - 59*fun(X[2],U[2]) + 37*fun(X[1],U[1]) - 9*fun(X[0],U[0]));
            y_p = U[3] + h/24*(9*fun(X[3]+h,y_p) + 19*fun(X[3],U[3]) - 5*fun(X[2],U[2]) + fun(X[1],U[1]));
            u.push_back(y_p);
            X[0]+=h; X[1]+=h; X[2]+=h; X[3]+=h;
            U[0]=U[1]; U[1]=U[2]; U[2]=U[3]; U[3]=y_p;
        }
        return u;
    }

    vector<double> improved_ab4am4(vector<double> X, vector<double> U, double h, double xup){
        vector<double> u={U[0], U[1], U[2], U[3]};
        double y_p, y_c, y;
        while (X[3]<=xup){
            y_p = U[3] + h/24*(55*fun(X[3],U[3]) - 59*fun(X[2],U[2]) + 37*fun(X[1],U[1]) - 9*fun(X[0],U[0]));
            y_c = U[3] + h/24*(9*fun(X[3]+h,y_p) + 19*fun(X[3],U[3]) - 5*fun(X[2],U[2]) + fun(X[1],U[1]));
            y = 251.0/270*y_c + 19.0/270*y_p;
            u.push_back(y);
            X[0]+=h; X[1]+=h; X[2]+=h; X[3]+=h;
            U[0]=U[1]; U[1]=U[2]; U[2]=U[3]; U[3]=y;
        }
        return u;
    }

    vector<double> compute(double x0, double h, double xup){
        vector<double> ut;
        for (double x=x0; x<xup+h; x+=h){
            ut.push_back(f(x));
        }
        return ut;
    }

};


int main() {
    SetConsoleOutputCP(65001); //更改cmd编码为utf8
    double u0=3;
    double h=0.1;
    double xup=1.5;
    double x0=0.0;

    vector<double> x0_ab4 = {x0, x0+h, x0+2*h, x0+3*h};
    vector<double> x0_ab4_am4 = {x0, x0+h, x0+2*h, x0+3*h};
    vector<double> x0_ipv = {x0, x0+h, x0+2*h, x0+3*h};
    vector<double> u_rk4, u_ab4, u_ab4am4, u_ipv, ut;

    cout << "------------RK4算法求解结果-------------" << endl;
    u_rk4 = Solution().rk4(x0, u0, h, xup);
    Solution().print_vector(u_rk4, x0, h);

    cout << "\n------------AB4算法求解结果-------------" << endl;
    vector<double> u_f4(u_rk4.begin(), u_rk4.begin()+4);
    u_ab4 = Solution().ab4(x0_ab4, u_f4, h, xup);
    Solution().print_vector(u_ab4, x0, h);

    cout << "\n----------AB4-AM4算法求解结果-----------" << endl;
    u_ab4am4 = Solution().ab4_am4(x0_ab4_am4, u_f4, h, xup);
    Solution().print_vector(u_ab4am4, x0, h);

    cout << "\n---------加速AB4-AM4算法求解结果---------" << endl;
    u_ipv = Solution().improved_ab4am4(x0_ipv, u_f4, h, xup);
    Solution().print_vector(u_ipv, x0, h);

    cout << "\n-----------精确解------------" << endl;
    ut = Solution().compute(x0, h, xup);
    Solution().print_vector(ut, x0, h);

    cout << "\n                          -------------------------综上，上述算法求解结果-------------------------" << endl;
    Solution().print_all(u_rk4, u_ab4, u_ab4am4, u_ipv, ut, x0, h);

    return 0;
}
