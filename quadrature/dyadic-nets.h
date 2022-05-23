#pragma once

#include <array>
#include <random>
#include "range.h"
#include "integrate.h"
#include "integrate-bins-stepper.h"
#include "vector-dimensions.h"

#include <filesystem>
#include <iostream>
#include <math.h>

#if (__cplusplus < 201703L)
namespace std {
    template< class T >
    inline constexpr bool is_integral_v = std::is_integral<T>::value;
}
#endif

using namespace std;

namespace viltrum {

template<typename RNG>
class StepperMonteCarloDyadicUniform {
    mutable RNG rng;

    mutable vector<array<float,2>> dyadic_net;
    mutable vector<int> dyadicIndex;
    mutable vector<array<float,2>> dyadicDims;

    template<typename Result>
    struct Samples {
        Result sumatory;
        unsigned long counter;
        Samples() : sumatory(0),counter(0) { }
    };
public:
    template<typename F, typename Float, std::size_t DIM>
    auto init(const F& f, const Range<Float,DIM>& range) const {
        //First element of tuple is sum, second element is number of samples
        return Samples<decltype(f(range.min()))>();
    }

    template<typename F, typename Float, std::size_t DIM, typename Result>
    void step(const F& f, const Range<Float,DIM>& range, Samples<Result>& samples) const {
    	std::array<Float,DIM> sample;
        int dimInd = 0;
	    for (std::size_t i=0;i<DIM;++i) {
            std::uniform_real_distribution<Float> dis(range.min(i),range.max(i));
            if(dyadicDims[dimInd][0] == i){
                sample[i] = dyadic_net[samples.counter][0];      //Get the dyadic-nets values
                sample[++i] = dyadic_net[samples.counter][1];
                if(dimInd<dyadicDims.size()-1) dimInd++;
            }
            else{
                
                sample[i] = dis(rng);
            }
	    }
        samples.sumatory += range.volume()*f(sample);
        ++samples.counter;
    }

    template<typename F, typename Float, std::size_t DIM, typename Result>
    Result integral(const F& f, const Range<Float,DIM>& range, const Samples<Result>& samples) const {
        return (samples.counter==0)?decltype(samples.sumatory)(0):(samples.sumatory/double(samples.counter));
    }

    StepperMonteCarloDyadicUniform(RNG&& r, vector<array<float,2>> dyadicDims_,int spp) :
        rng(std::forward<RNG>(r)), dyadicDims(dyadicDims_) {
            dyadicIndex.resize(dyadicDims_.size());
            int i;
            for(i=0; i<=12; i++){
                if(spp == pow(2,i)) break;
            }
            if(i>12) {
                std::cout<<"Incorrect number of samples: " << spp <<"\nPlease, use a power of 2 between 0 and 12."<<std::endl;
                exit(0);
            }
            string path = "../external/viltrum/utils/Blue-Nets/";
            if(i<=3){
                path = path + "000" + to_string((int) pow(2,i)) + ".txt";
            }
            else if(i<=6){
                path = path + "00" + to_string((int) pow(2,i)) + ".txt";
            }
            else if(i<=9){
                path = path + "0" + to_string((int) pow(2,i)) + ".txt";
            }
            else{
                path = path + to_string((int) pow(2,i)) + ".txt";
            }
            //0001.txt";
            ifstream filein(path);
            int k = 1;
            cout<<path<<endl;
            if (filein.is_open()) {
                string line;
                getline(filein, line);
                while (getline(filein, line)) {
                    //cout<<line<<endl;
                    std::size_t pos = line.find(" ");      // position of "live" in str
                    std::string f2 = line.substr(pos);
                    //cout<<stof(line)<<" "<<stof(f2)<<endl;
                    dyadic_net.push_back({stof(line),stof(f2)});
                }
            }
            
        }

    StepperMonteCarloDyadicUniform(RNG&& r) :
    rng(std::forward<RNG>(r)) { }
};

template<typename RNG>
auto stepper_monte_carlo_dyadic_uniform(RNG&& rng) {
    return StepperMonteCarloDyadicUniform<RNG>(std::forward<RNG>(rng));
}

template<typename RNG>
auto stepper_monte_carlo_dyadic_uniform(RNG&& rng,vector<array<float,2>> dyadicDims, int spp) {
    return StepperMonteCarloDyadicUniform<RNG>(std::forward<RNG>(rng),dyadicDims,spp);
}

auto stepper_monte_carlo_dyadic_uniform(std::size_t seed = std::random_device()()) {
    return stepper_monte_carlo_dyadic_uniform(std::mt19937_64(seed));
}

auto stepper_monte_carlo_dyadic_uniform(vector<array<float,2>> dyadicDims, int spp, std::size_t seed = std::random_device()()) {
    return stepper_monte_carlo_dyadic_uniform(std::mt19937_64(seed),dyadicDims,spp);
}

template<typename RNG>
auto integrator_monte_carlo_dyadic_uniform(RNG&& rng, unsigned long samples, 
    std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return integrator_stepper(stepper_monte_carlo_dyadic_uniform(std::forward<RNG>(rng)),samples);
}

auto integrator_monte_carlo_dyadic_uniform(unsigned long samples, std::size_t seed = std::random_device()()) {
    int i;
    for(i=0; i<=12; i++){
        if(samples == pow(2,i)) break;
    }
    if(i>12) {
        std::cout<<"Incorrect number of samples, please, use a power of 2 between 0 and 12."<<std::endl;
        exit(0);
    }
    else {
        std::cout<<"Using a "<<pow(2,i)<<" net."<<std::endl;
    }
    return integrator_stepper(stepper_monte_carlo_dyadic_uniform(seed),samples);
}


}

