#include <bits/stdc++.h>
#include <filesystem>
#include "tu_algo.hpp"
extern "C" {
  #include <sleef.h>
  #include <crlibm.h>
}

static inline uint64_t tobits(double x){ uint64_t u; std::memcpy(&u,&x,sizeof u); return u; }
static inline double ulp_dist(double a, double b){
    if (std::isnan(a) || std::isnan(b)) return std::numeric_limits<double>::infinity();
    if (std::isinf(a) || std::isinf(b)) return (a==b)?0.0:std::numeric_limits<double>::infinity();
    int64_t ia = (int64_t)tobits(a); if (ia<0) ia = 0x8000000000000000LL - ia;
    int64_t ib = (int64_t)tobits(b); if (ib<0) ib = 0x8000000000000000LL - ib;
    return (double)std::llabs(ia-ib);
}

int main(){
    std::mt19937_64 rng(321);
    std::uniform_real_distribution<double> U01(0.0,1.0);

    const size_t N = 1<<20;
    std::vector<double> A(N), Nth(N);
    for(size_t i=0;i<N;i++){
        double u=U01(rng), v=U01(rng);
        double e = -500 + u*(1000);
        A[i] = std::ldexp(1.0,(int)e) * (1.0 + 1e-12*v);
        double choices[5]={0.5,2.0,3.0,10.0,100.0};
        Nth[i]=choices[(int)std::floor(v*5)%5];
        if (A[i]<0.0 && Nth[i]!=2.0) { A[i]=std::abs(A[i]); }
    }

    auto run = [&](const char* name, auto F){
        using clk=std::chrono::high_resolution_clock;
        std::vector<double> ulps; ulps.reserve(N);
        double ulpmax=0; size_t bad=0;
        auto t0=clk::now();
        for(size_t i=0;i<N;i++){
            double got = F(A[i], Nth[i]);
            double ref = std::pow(A[i], 1.0/Nth[i]);
            double u = ulp_dist(ref, got);//MARCARA DIFERENCIA DECIMALESLa función ulp_dist calcula la distancia en unidades de último lugar (ULP) entre dos números de punto flotante de tipo double. Maneja casos especiales como valores NaN o infinitos, devolviendo infinito si alguno de los valores es NaN o si los infinitos son diferentes, y utiliza la representación binaria de los números para calcular la distancia.
            ulps.push_back(u); ulpmax = std::max(ulpmax,u);
            bad += std::isnan(got);
        }
        auto t1=clk::now();
        std::sort(ulps.begin(), ulps.end());
        double ulp95 = ulps[(size_t)(0.95*N)];
        double ns = std::chrono::duration<double,std::nano>(t1-t0).count()/N;
        std::filesystem::create_directories("results");
        std::ofstream f(std::string("results/")+name+".csv");
        f<<"name,time_ns,ulp95,ulpmax,bad\n"<<name<<","<<ns<<","<<ulp95<<","<<ulpmax<<","<<bad<<"\n";
        std::cerr<<name<<"  "<<ns<<" ns/call  ULP95="<<ulp95<<"  ULPmax="<<ulpmax<<"  bad="<<bad<<"\n";
    };

    run("ROOT_TU", [&](double a,double n){ return tu_root(a,n); });
    run("ROOT_LIBM",    [&](double a,double n){ return std::pow(a, 1.0/n); });
    run("ROOT_SLEEF",   [&](double a,double n){ return std::exp(Sleef_log_u10(a)/n);} );
    run("ROOT_CRLIBM",  [&](double a,double n){ return std::exp(log_rn(a)/n);} );
    return 0;
}
