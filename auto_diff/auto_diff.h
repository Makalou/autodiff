//
// Created by 王泽远 on 2023/1/8.
//

#ifndef AUTODIFF_AUTO_DIFF_H
#define AUTODIFF_AUTO_DIFF_H

#include "forward/dual.h"
#include "reverse/node.h"

namespace ad{
    enum class differential_mode{
        FORWARD,
        REVERSE
    };

    template<differential_mode mod = differential_mode::FORWARD>
    struct differentiable_var{};

    template<>
    struct differentiable_var<differential_mode::FORWARD>{
        dual_number _dual;
    };

    template<>
    struct differentiable_var<differential_mode::REVERSE>{
        std::shared_ptr<differential_node> _node;
    };

    namespace detail{
        using dvf = differentiable_var<differential_mode::FORWARD>;
        using dvr = differentiable_var<differential_mode::REVERSE>;
    }

    inline detail::dvf operator+(const  detail::dvf &l, const  detail::dvf &r) {
        return detail::dvf{l._dual + r._dual};
    }

    inline detail::dvf operator-(const  detail::dvf &l, const  detail::dvf &r) {
        return detail::dvf{l._dual-r._dual};
    }

    inline detail::dvf operator*(const  detail::dvf &l, const  detail::dvf &r) {
        return detail::dvf{l._dual*r._dual};
    }

    inline detail::dvf operator/(const  detail::dvf &l, const  detail::dvf &r) {
        return detail::dvf{l._dual/r._dual};
    }

    inline detail::dvf operator+(double l, const  detail::dvf &r) {
        return detail::dvf{l+r._dual};
    }

    inline detail::dvf operator-(double l, const  detail::dvf &r) {
        return detail::dvf{l-r._dual};
    }

    inline detail::dvf operator*(double l, const  detail::dvf &r) {
        return detail::dvf{l*r._dual};
    }

    inline detail::dvf operator/(double l, const  detail::dvf &r) {
        return detail::dvf{l/r._dual};
    }

    inline detail::dvf operator+(const  detail::dvf &l, double r){
        return detail::dvf{l._dual+r};
    }

    inline detail::dvf operator-(const  detail::dvf &l, double r){
        return detail::dvf{l._dual-r};
    }

    inline detail::dvf operator*(const  detail::dvf &l, double r){
        return detail::dvf{l._dual*r};
    }

    inline detail::dvf operator/(const  detail::dvf &l, double r){
        return detail::dvf{l._dual/r};
    }

    inline detail::dvf sin( detail::dvf d) {
        return detail::dvf{sin(d._dual)};
    }

    inline detail::dvf cos( detail::dvf d) {
        return detail::dvf{cos(d._dual)};
    }

    inline detail::dvf tan( detail::dvf d) {
        return detail::dvf{tan(d._dual)};
    }

    inline detail::dvf pow( detail::dvf d, double n){
        return detail::dvf{pow(d._dual,n)};
    }

    inline detail::dvf pow(double a,  detail::dvf d){
        return detail::dvf{pow(a,d._dual)};
    }

    inline detail::dvf sqrt( detail::dvf d){
        return detail::dvf{sqrt(d._dual)};
    }

    inline detail::dvf exp( detail::dvf d){
        return detail::dvf{exp(d._dual)};
    }

    inline detail::dvf ln( detail::dvf d){
        return detail::dvf{ln(d._dual)};
    }

    inline detail::dvf asin( detail::dvf d){
        return detail::dvf{asin(d._dual)};
    }

    inline detail::dvf acos( detail::dvf d){
        return detail::dvf{acos(d._dual)};
    }

    inline detail::dvf atan( detail::dvf d){
        return detail::dvf{atan(d._dual)};
    }

    inline detail::dvr operator+(const  detail::dvr &l, const  detail::dvr &r){
        return detail::dvr{l._node+r._node};
    }

    static inline std::pair<double,double> gradient_at(const std::function<detail::dvf( detail::dvf)>& f,double x){
        auto u_dual = detail::dvf{ad::dual_number{x,1}};
        auto res = f(u_dual)._dual;
        return {res._real_part,res._dual_part};
    }

    static inline std::pair<double,std::pair<double,double>> gradient_at(const std::function< detail::dvf( detail::dvf, detail::dvf)>& f,double x,double y){

        auto dfdx = f(detail::dvf{dual_number{x,1}},detail::dvf{dual_number{y,0}})._dual;
        auto dfdy = f(detail::dvf{dual_number{x,0}},detail::dvf{dual_number{y,1}})._dual;

        assert(dfdx._real_part == dfdy._real_part);
        return {dfdx._real_part,{dfdx._dual_part,dfdy._dual_part}};
    }
}
#endif //AUTODIFF_AUTO_DIFF_H
