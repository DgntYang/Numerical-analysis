#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>

using namespace std;
class Solution{
public:
    double fun(double x){
        return pow(x, 3)/3 - x;
        // return sin(x);
    }

    double autograd(double x, double step){
        return (fun(x) - fun(x-step))/step;
    }

    double grad(double x, double step){
        return (pow(x, 2) - 1);
    }

    vector<double> find_ranges(int start, int end, double step){
        vector <double> ranges;
        for (double i=start; i<=end; i+=step){
            if (fun(i) * fun(i+step) < 0){
                ranges.push_back(i);
                ranges.push_back(i+step);
            }
        }
        return ranges;
    }

    double iterate(double x0, int maxiter, double epsilon, double grad_step=0.02){
        double x1;
        for (int i=0; i<= maxiter; i++){
            x1 = x0 - fun(x0) / grad(x0, grad_step);
            if (i%maxiter == 0) {
                //cout << "Epoch " << i + 1 << ",  x: "  << setprecision(10) << x1 << endl;
            }
            if (abs(x1 - x0) <= epsilon){
                break;
            }
            x0 = x1;
        }
        return x1;
    }
};

int main() {
    // Init the solution
    double x, x0, search_step = 0.02, grad_step = 0.02;
    int maxiter = 50;
    int start = -10, end = 10;
    double epsilon = pow(10, -8);
    double delta;
    double r, r_;
    vector<double> ans;

    //零点定理寻找零点所在区间
    vector<double> ranges = Solution().find_ranges(start, end, search_step);
    int N = ranges.size()/2; //区间的总对数

    // Start solving the equation
    for (int i=0; i<N; i++){
        cout << "-----------------------------" << endl;
        cout << "Start solving "<< i << "-th solution." << endl;
        x0 = (ranges[2*i] + ranges[2*i+1])/2;
        x = Solution().iterate(x0, maxiter, epsilon, grad_step);
        ans.push_back(x);
    }

    // Output the results
    cout << "-------------------------\n" << "Solutions for equation are \n" ;
    for (auto j=ans.begin(); j != ans.end(); j++){
        cout << setprecision(10) << *j << endl;
    }

    // compute the delta

    cout << "-------------------------\n" << "Starting computing delta ....\n" ;
    double delta_start = -0.8, delta_end = 0.8;
    int count = 0;
    for (double x=delta_start; x<=delta_end; x+=0.0001){
        count ++;
        r = Solution().iterate(x, maxiter, epsilon);

        //if (abs(r-r_) >= 1.7 && count >= 2){
        //    cout << "-----------------------------" << endl;
        //    cout << "delta: "<< x-0.0001 << "\t  r: " << r_ << endl;
        cout << "delta: "<< x << "\t  r: " << r << endl;
        //}
        r_ = r;
        //if (abs(r) >= pow(10,-4)){
        //    break;
        //}
        //if (abs(r) <= pow(10,-4)){
        //    break;
        //}
    }

    cout << "Completed." << endl;



    return 0;
}
