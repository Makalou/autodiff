#include <iostream>
#include "auto_diff/auto_diff.h"

using dual = ad::dual_number;

dual func(dual x){
    return 4.0*ad::pow(x,2.0)+5.0;
}

dual func1(dual x){
    return 3*ad::pow(func(x),2.0)+1;
}

double func1_truth(double x){
    auto f = 4*x*x+5;
    return 3*(f)*(f)+1;
}

double func1_dx_truth(double x){
    return 192*x*x*x + 240*x;
}

dual g(dual x){
    return 0.01*pow(x,4)+pow(3,x)+0.001*pow(x,2);
}

template<bool silence = true>
std::pair<double,double> gradient_descent_optimizer(const std::function<dual(dual)>& f,double init,double threshold,double step_length){
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

dual func2D(dual x, dual y){
    return x*x + y*y*y;
}

int main() {

    double x = 3;

    auto [output, dx] = ad::gradient_at(func1,x);

    assert(output == func1_truth(x));
    assert(dx == func1_dx_truth(x));

    std::cout << output <<","<<dx<< std::endl;

    auto [minimum,xx] = gradient_descent_optimizer<false>(g,3,1e-10,0.0001);

    std::cout << minimum <<" at "<<xx<< std::endl;

    auto[output1,grad] = ad::gradient_at(func2D,2,2);

    std::cout << output1 <<","<<"["<<grad.first<<","<<grad.second<<"]"<< std::endl;

    return 0;
}
