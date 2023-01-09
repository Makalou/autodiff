//
// Created by 王泽远 on 2023/1/9.
//

#include "gtest/gtest.h"
#include "../auto_diff/auto_diff.h"
#include <random>

using Double = ad::differentiable_var<ad::differential_mode::FORWARD>;

Double g(Double x){
    return 0.01*pow(x,4)+pow(3,x)+0.001*pow(x,2);
}

template<bool silence = true>
std::pair<double,double> gradient_descent_optimizer(const std::function<Double(Double)>& f, double init, double threshold, double step_length){
    auto x = init;

    auto [output,grad] = ad::gradient_at(f,x);

    while(std::abs(grad)>threshold){
        x -= step_length * grad;
        auto new_result = ad::gradient_at(f,x);

        output = new_result.first;
        grad = new_result.second;
        if constexpr (!silence){
            std::cout << output <<" , "<<grad<<" , "<<x<< std::endl;
        }
    }

    return {output,x};
}

TEST(ScalarFunction1DOptimization, BasicAssertions) {

    auto [minimum,xx] = gradient_descent_optimizer(g,3,1e-10,0.0001);

    EXPECT_LT(abs(0.24-minimum),1e-4);
    EXPECT_LT(abs(-1.6443-xx),1e-2);
}