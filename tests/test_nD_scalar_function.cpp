//
// Created by 王泽远 on 2023/1/21.
//
#include "gtest/gtest.h"
#include "../auto_diff/auto_diff.h"
#include <random>
#include <array>

using Double = ad::differentiable_var<ad::differential_mode::REVERSE>;

using Vector = std::array<Double,4>;

Double funcND(Vector v){
    return v[0] + pow(v[1],2) + pow(v[2],3) + pow(v[3],4);
}

double funcND_truth(double x,double y,double z,double w){
    return x + pow(y,2) + pow(z,3) + pow(w,4);
}

std::vector<double> funcND_dxdy_truth(double x, double y,double z,double w){
    return {1,2*y,3*pow(z,2),4*pow(w,3)};
}

TEST(ScalarFunctionnD, BasicAssertions) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    for (int n = 0; n < 100; ++n) {
        double x = dis(gen);
        double y = dis(gen);
        double z = dis(gen);
        double w = dis(gen);

        auto res = ad::gradient_at<4>(funcND,{x,y,z,w});

        EXPECT_DOUBLE_EQ(res.first, funcND_truth(x,y,z,w));

        auto grad_truth = funcND_dxdy_truth(x,y,z,w);

        EXPECT_DOUBLE_EQ(res.second[0], grad_truth[0]);
        EXPECT_DOUBLE_EQ(res.second[1], grad_truth[1]);
        EXPECT_DOUBLE_EQ(res.second[2], grad_truth[2]);
        EXPECT_DOUBLE_EQ(res.second[3], grad_truth[3]);
    }

}
