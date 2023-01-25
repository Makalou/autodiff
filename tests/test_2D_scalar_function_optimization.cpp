//
// Created by 王泽远 on 2023/1/9.
//
#include "gtest/gtest.h"
#include "../auto_diff/auto_diff.h"
#include <random>

using Double = ad::differentiable_var<ad::differential_mode::FORWARD>;

inline Double target2D(Double x, Double y){
    return pow(x-2,2) + pow(y-1,4)+10.0;
}

template<bool silence = true>
std::pair<double,std::pair<double,double>>
gradient_descent_optimizer(const std::function<Double(Double, Double)>& f, double init_x, double init_y, double threshold, double step_length){
    auto x = init_x;
    auto y = init_y;

    unsigned int descent_axis = 0;
    auto [output,grad] = ad::value_and_gradient_at(f,x,y);

    while(std::abs(grad.first)>threshold||std::abs(grad.second)>threshold){

        if(descent_axis == 0){
            x -= step_length * grad.first;
        }else{
            y -= step_length * grad.second;
        }

        descent_axis = (descent_axis+1)%2;

        auto new_result = ad::value_and_gradient_at(f,x,y);

        output = new_result.first;
        grad = new_result.second;

        if constexpr (!silence){
            std::cout << output <<" , "<<"["<<grad.first<<","<<grad.second<<"]"<<" , "<<"["<<x<<","<<y<<"]"<< std::endl;
        }
    }

    return {output,{x,y}};
}

TEST(ScalarFunction2DOptimization, BasicAssertions) {

    auto [minimum,xxyy] = gradient_descent_optimizer(target2D,3,3,1e-6,0.01);

    EXPECT_LT(abs(10-minimum),1e-5);
    EXPECT_LT(abs(2-xxyy.first),1e-2);
    EXPECT_LT(abs(1-xxyy.second),1e-2);
}
