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
    std::mt19937_64 rng(123);
    std::uniform_real_distribution<double> U01(0.0,1.0);

    const size_t N = 1<<20;
    std::vector<double> A(N), B(N);
    for(size_t i=0;i<N;i++){
        double u=U01(rng), v=U01(rng);
        double e = -500 + u*(1000);
        A[i] = std::ldexp(1.0,(int)e) * (1.0 + 1e-12*v);
        B[i] = 1.0 + v*15.0;
    }

    auto run = [&](const char* name, auto F){
        using clk=std::chrono::high_resolution_clock;
        double ulp95=0, ulpmax=0; size_t bad=0; double sum=0;
        std::vector<double> ulps; ulps.reserve(N);
        auto t0=clk::now();
        for(size_t i=0;i<N;i++){
            double a=A[i], b=B[i];
            double got = F(a,b);
            double ref = std::log(a)/std::log(b);
            double u = ulp_dist(got, ref);
            ulps.push_back(u); ulpmax = std::max(ulpmax,u);
            bad += std::isnan(got);
            sum += got;
        }
        auto t1=clk::now();
        std::sort(ulps.begin(), ulps.end());
        ulp95 = ulps[(size_t)(0.95*N)];
        double ns = std::chrono::duration<double,std::nano>(t1-t0).count()/N;
        std::filesystem::create_directories("results");
        std::ofstream f(std::string("results/")+name+".csv");
        f<<"name,time_ns,ulp95,ulpmax,bad\n"<<name<<","<<ns<<","<<ulp95<<","<<ulpmax<<","<<bad<<"\n";
        std::cerr<<name<<"  "<<ns<<" ns/call  ULP95="<<ulp95<<"  ULPmax="<<ulpmax<<"  bad="<<bad<<"\n";
    };

    run("TU_ALGO", [&](double a,double b){ return tu_log_b_a(a,b); });
    run("LIBM",    [&](double a,double b){ return std::log(a)/std::log(b); });
    run("CRLIBM",  [&](double a,double b){ return log_rn(a)/log_rn(b); });
    run("SLEEF",   [&](double a,double b){ return Sleef_log_u10(a)/Sleef_log_u10(b);} );
    return 0;
}
