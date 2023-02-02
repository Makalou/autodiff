//
// Created by 王泽远 on 2023/2/2.
//

#ifndef AUTODIFF_EVALUATOR_H
#define AUTODIFF_EVALUATOR_H

#include "auto_diff.h"

namespace ad{
    template<differential_mode mod = differential_mode::FORWARD>
    struct evaluator{};

    template<>
    struct evaluator<differential_mode::FORWARD>
    {
        evaluator(){

        }
    };

    template<>
    struct evaluator<differential_mode::REVERSE>
    {
        evaluator(){

        }

        float value_at(double x) {

        }

        float gradient_at(double x){

        }

        std::pair<double,double> value_and_gradient_at(double x){

        }

        differential_node z;
    };
}

#endif //AUTODIFF_EVALUATOR_H
