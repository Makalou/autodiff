//
// Created by 王泽远 on 2023/1/9.
//
#include "gtest/gtest.h"
#include "../auto_diff/auto_diff.h"
#include <random>

using Double = ad::differentiable_var<ad::differential_mode::FORWARD>;

Double func2D(Double x, Double y){
    return pow(x-2,2) + pow(y-1,4)+10.0;
}

double func2D_truth(double x,double y){
    return pow(x-2,2) + pow(y-1,4)+10.0;
}

std::pair<double,double> func2D_dxdy_truth(double x, double y){
    return {2*(x-2),4* pow(y-1,3)};
}

TEST(ScalarFunction2D, BasicAssertions) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    for (int n = 0; n < 100; ++n) {
        double x = dis(gen);
        double y = dis(gen);

        auto res = ad::gradient_at(func2D,x,y);

        EXPECT_DOUBLE_EQ(res.first, func2D_truth(x,y));

        auto [dfdx,dfdy] = func2D_dxdy_truth(x,y);
        EXPECT_DOUBLE_EQ(res.second.first, dfdx);
        EXPECT_DOUBLE_EQ(res.second.second, dfdy);
    }

}