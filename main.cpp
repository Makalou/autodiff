#include <iostream>
#include "auto_diff/auto_diff.h"

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
    return 192*x*x*x + 240*x;
}

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

Double func2D(Double x, Double y){
    return pow(x-2,2) + pow(y-1,4)+10.0;
}

template<bool silence = true>
std::pair<double,std::pair<double,double>> gradient_descent_optimizer(const std::function<Double(Double, Double)>& f, double init_x, double init_y, double threshold, double step_length){
    auto x = init_x;
    auto y = init_y;

    uint descent_axis = 0;
    auto [output,grad] = ad::gradient_at(f,x,y);

    while(std::abs(grad.first)>threshold||std::abs(grad.second)>threshold){

        if(descent_axis == 0){
            x -= step_length * grad.first;
        }else{
            y -= step_length * grad.second;
        }

        descent_axis = (descent_axis+1)%2;

        auto new_result = ad::gradient_at(f,x,y);

        output = new_result.first;
        grad = new_result.second;

        if constexpr (!silence){
            std::cout << output <<" , "<<"["<<grad.first<<","<<grad.second<<"]"<<" , "<<"["<<x<<","<<y<<"]"<< std::endl;
        }
    }

    return {output,{x,y}};
}

int main() {

    double x = 3;

    auto [output, dx] = ad::gradient_at(func1,x);

    assert(output == func1_truth(x));
    assert(dx == func1_dx_truth(x));

    std::cout << output <<","<<dx<< std::endl;

    auto [minimum,xx] = gradient_descent_optimizer<true>(g,3,1e-10,0.0001);

    std::cout << minimum <<" at "<<xx<< std::endl;

    auto[output1,grad] = ad::gradient_at(func2D,2,2);

    std::cout << output1 <<","<<"["<<grad.first<<","<<grad.second<<"]"<< std::endl;

    auto [minimum1,xxyy] = gradient_descent_optimizer<true>(func2D,3,3,1e-6,0.01);

    std::cout<<minimum1<<" at "<<"["<<xxyy.first<<","<<xxyy.second<<"]"<< std::endl;

    return 0;
}
