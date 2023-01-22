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

using Vector100 = std::array<Double,100>;

Double f1(Vector100 k,std::array<double,100> x){

    auto ret = k[0];

    for(int i = 1;i<100;++i){
        ret = ret + x[i-1]*k[i];
    }

    return ret;
}

Double loss(Vector100 v){
    Double lse;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    std::array<double,100> x{};

    for(int n =0;n<1000;++n){
        x[0] = dis(gen);
        double y = sin(x[0]);
        for(int i =1;i<100;++i){
            x[i] = x[i-1]*x[0];
        }
        lse = lse + pow(f1(v,x) - y,2);
    }

    return lse;
}

Double sphere100D(Vector100 v){
    Double ret;
    for(int i =0;i<v.size();++i){
        ret = ret + pow(v[i] - i,2);
    }
    return ret + 1;
}

template<int n ,bool silence = true>
std::pair<std::pair<double,std::array<double,n>>,std::array<double,n>>
gradient_descent_optimizer(const std::function<Double(std::array<Double ,n>)>& f,
                           std::array<double,n> init_val,
                           double threshold, double step_length){

    auto v  = init_val;

    unsigned int descent_axis = 0;
    auto result = ad::gradient_at<n>(f,v);
    auto output = result.first;
    auto grad = result.second;

    while(std::any_of(std::begin(grad),std::end(grad),[threshold](double g){return abs(g) > threshold;})){

        v[descent_axis] -= step_length * grad[descent_axis];

        descent_axis = (descent_axis+1)%n;

        result = ad::gradient_at<n>(f,v);

        output = result.first;
        grad = result.second;

        if constexpr (!silence){
            //std::cout << output <<" , "<<"["<<grad.first<<","<<grad.second<<"]"<<" , "<<"["<<x<<","<<y<<"]"<< std::endl;
        }
    }

    return {{output,v},grad};
}

TEST(ScalarFunctionnDOptimization, BasicAssertions) {
    std::array<double,100> init_v{0};
    auto [output,grad] = gradient_descent_optimizer<100>(sphere100D,init_v,1e-6,0.01);
    EXPECT_LT(abs(1-output.first),1e-5);
}
