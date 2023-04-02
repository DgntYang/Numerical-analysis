//
// Created by yangtan on 2022/9/8.
//
#include "stdio.h"
#include "math.h"
#include "iostream"
#include "iomanip"

using namespace std;

class solution{
public:
    void compute(int N){
        //compute baseline
        float res = 0;
        res = 0.5 * (1.5 - 1.0/N - 1.0/(N+1));
        cout << "Baseline: \t" << fixed << setprecision(7)  << res << endl;

        // compute in ascent
        res = 0;
        for (int i=2; i<= N; i++){
            res += 1.0 / (pow(i, 2)-1);
        }
        cout << "Ascent: \t" << fixed << setprecision(7) <<res << endl;

        // compute in descent
        res = 0;
        for (int i=N; i >= 2; i--){
            res += 1.0 / (pow(i, 2)-1);
        }
        cout << "Descent: \t" << fixed << setprecision(7) << res << "\n" <<endl;
    }
};


int main(){
    int N[3] = {100, 10000,1000000};
    for (int j = 0; j<=2; j++){
        cout << "start" << " epoch " << j+1 << endl;
        solution().compute(N[j]);
    }
    return 0;
}









