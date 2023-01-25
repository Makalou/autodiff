//
// Created by 王泽远 on 2023/1/22.
//

//
// Created by 王泽远 on 2023/1/21.
//
#include "gtest/gtest.h"
#include "../auto_diff/auto_diff.h"
#include <random>
#include <array>

using Double = ad::differentiable_var<ad::differential_mode::FORWARD>;

using Vector100 = std::array<Double ,100>;
const int N = 3;
using VectorN = std::array<Double,N>;

Double f1(VectorN k,const std::array<double,k.size()>& x){
    Double ret{};
    for(int i = 0;i<k.size();++i){
        ret = ret + x[i]*k[i];
    }
    return ret;
}

Double loss(VectorN k,const std::array<std::pair<double,double>,10>& samples){
    std::array<double,k.size()> x{};
    Double lse{};

    for(const auto & sample : samples){
        double y = sample.second;
        x[0] = 1.;
        x[1] = sample.first;
        for(int i = 2; i<k.size(); ++i){
            x[i] = x[i-1]*x[1];
        }
        lse = lse + pow(f1(k,x)-y,2); // Least Square Error
    }

    return lse;
}

Double sphere100D(Vector100 v){
    Double ret;
    for(int i =0;i<100;++i){
        ret = ret + pow(v[i] - i,2);
    }
    return ret + 1.0;
}

template<int n ,bool silence = true>
std::pair<std::pair<double,std::array<double,n>>,std::array<double,n>>
gradient_descent_optimizer(const std::function<Double(std::array<Double ,n>)>& f,
                           std::array<double,n> init_val,
                           double threshold, double step_length){

    auto v  = init_val;

    unsigned int descent_axis = 0;
    auto result = ad::value_and_gradient_at<n>(f,v);
    auto output = result.first;
    auto grad = result.second;

    while(std::any_of(std::begin(grad),std::end(grad),[threshold](double g){return abs(g) > threshold;})){

        v[descent_axis] -= step_length * grad[descent_axis];

        descent_axis = (descent_axis+1)%n;

        result = ad::value_and_gradient_at<n>(f,v);

        output = result.first;
        grad = result.second;

        if constexpr (!silence){
            std::cout << output << std::endl;
        }
    }

    return {{output,v},grad};
}

TEST(ScalarFunctionnDOptimization, BasicAssertions) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-10.0, 10.0);
    std::uniform_real_distribution<> dis2(-0.1, 0.1);

    std::array<std::pair<double,double>,10> samples;

    for(auto & sample : samples){
        double x = dis(gen);

        double a = dis2(gen)+3.0;
        double b = dis2(gen)+2.0;
        double c = dis2(gen)+1.0;

        sample = std::make_pair(x,a*x*x+b*x+c);
    }

    std::array<double,N> init_v{};

    auto [output,grad] = gradient_descent_optimizer<N>([samples](auto && PH1) { return loss(std::forward<decltype(PH1)>(PH1), samples); },
                                                         init_v,1e-6,0.00001);
    //EXPECT_LT(abs(1-output.first),1e-5);
    std::cout<<output.first<<std::endl;
    for(int i = 0;i<N;++i){
        std::cout<<output.second[i]<<",";
    }
    std::cout<<std::endl;
}
