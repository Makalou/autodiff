//
// Created by 王泽远 on 2023/1/9.
//

#include "gtest/gtest.h"
#include "../auto_diff/auto_diff.h"
#include <random>

using Double = ad::differentiable_var<ad::differential_mode::FORWARD>;

Double func(Double x){
    return 4.0*ad::pow(x,2.0)+5.0;
}

Double func1(Double x){
    return 3*ad::pow(func(x),2.0)+1;
}

double func1_truth(double x){
    auto f = 4*x*x+5;
    return 3*(f)*(f)+1;
}

double func1_dx_truth(double x){
    return 192*pow(x,3) + 240*x;
}

Double xpowx(Double x){
    return ad::pow(x,x);
}

double xpowx_truth(double x){
    return std::pow(x,x);
}

double xpowx_dx_truth(double x){
    return pow(x,x)*(log(x)+1);
}

TEST(ScalarFunction1D, BasicAssertions) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    for (int n = 0; n < 100; ++n) {
        double x = dis(gen);
        auto res = ad::value_and_gradient_at(func1,x);
        EXPECT_DOUBLE_EQ(res.first,func1_truth(x));
        EXPECT_DOUBLE_EQ(res.second, func1_dx_truth(x));
    }

    dis = std::uniform_real_distribution<double>(0.00001, 100.0);

    for (int n = 0; n < 100; ++n) {
        double x = dis(gen);
        auto res = ad::value_and_gradient_at(xpowx,x);
        EXPECT_DOUBLE_EQ(res.first, xpowx_truth(x));
        EXPECT_DOUBLE_EQ(res.second, xpowx_dx_truth(x));
    }

}
