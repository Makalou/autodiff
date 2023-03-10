//
// Created by 王泽远 on 2023/1/8.
//

#ifndef AUTODIFF_AUTO_DIFF_H
#define AUTODIFF_AUTO_DIFF_H

#include <functional>
#include <cassert>

#include "forward/dual.h"
#include "reverse/node_variant.h"

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
        detail::nsp _node;

        differentiable_var(){
            _node = std::make_shared<differential_node>(constant_node{.0});
        }

        differentiable_var(double v,uint idx){
            _node = std::make_shared<differential_node>(variable_node{v,idx});
        }

        explicit differentiable_var<differential_mode::REVERSE>(const detail::nsp& node){
            _node = node;
        }

    };

    namespace detail{
        using dvf = differentiable_var<differential_mode::FORWARD>;
        using dvr = differentiable_var<differential_mode::REVERSE>;
    }

    //Begin Forward Differential Variable Operator Overloading

    inline detail::dvf operator+(const detail::dvf &l, const detail::dvf &r) {
        return detail::dvf{l._dual + r._dual};
    }

    inline detail::dvf operator-(const detail::dvf &l, const detail::dvf &r) {
        return detail::dvf{l._dual - r._dual};
    }

    inline detail::dvf operator*(const detail::dvf &l, const detail::dvf &r) {
        return detail::dvf{l._dual * r._dual};
    }

    inline detail::dvf operator/(const detail::dvf &l, const detail::dvf &r) {
        return detail::dvf{l._dual / r._dual};
    }

    inline detail::dvf operator+(double l, const detail::dvf &r) {
        return detail::dvf{l + r._dual};
    }

    inline detail::dvf operator-(double l, const detail::dvf &r) {
        return detail::dvf{l - r._dual};
    }

    inline detail::dvf operator*(double l, const detail::dvf &r) {
        return detail::dvf{l * r._dual};
    }

    inline detail::dvf operator/(double l, const detail::dvf &r) {
        return detail::dvf{l / r._dual};
    }

    inline detail::dvf operator+(const detail::dvf &l, double r) {
        return detail::dvf{l._dual + r};
    }

    inline detail::dvf operator-(const detail::dvf &l, double r) {
        return detail::dvf{l._dual - r};
    }

    inline detail::dvf operator*(const detail::dvf &l, double r) {
        return detail::dvf{l._dual * r};
    }

    inline detail::dvf operator/(const detail::dvf &l, double r) {
        return detail::dvf{l._dual / r};
    }

    inline detail::dvf sin(detail::dvf d) {
        return detail::dvf{sin(d._dual)};
    }

    inline detail::dvf cos(detail::dvf d) {
        return detail::dvf{cos(d._dual)};
    }

    inline detail::dvf tan(detail::dvf d) {
        return detail::dvf{tan(d._dual)};
    }

    inline detail::dvf pow(detail::dvf d, double n) {
        return detail::dvf{pow(d._dual, n)};
    }

    inline detail::dvf pow(double a, detail::dvf d) {
        return detail::dvf{pow(a, d._dual)};
    }

    inline detail::dvf pow(detail::dvf d1, detail::dvf d2){
        return detail::dvf{pow(d1._dual, d2._dual)};
    }

    inline detail::dvf sqrt(detail::dvf d) {
        return detail::dvf{sqrt(d._dual)};
    }

    inline detail::dvf exp(detail::dvf d) {
        return detail::dvf{exp(d._dual)};
    }

    inline detail::dvf ln(detail::dvf d) {
        return detail::dvf{ln(d._dual)};
    }

    inline detail::dvf asin(detail::dvf d) {
        return detail::dvf{asin(d._dual)};
    }

    inline detail::dvf acos(detail::dvf d) {
        return detail::dvf{acos(d._dual)};
    }

    inline detail::dvf atan(detail::dvf d) {
        return detail::dvf{atan(d._dual)};
    }

    //END Forward Differential Variable Operator Overloading

    //Begin Reverse Differential Variable Operator Overloading
    inline detail::dvr operator+(const  detail::dvr &l, const  detail::dvr &r){
        return detail::dvr{l._node+r._node};
    }

    inline detail::dvr operator-(const detail::dvr &l, const detail::dvr &r) {
        return detail::dvr{l._node - r._node};
    }

    inline detail::dvr operator*(const detail::dvr &l, const detail::dvr &r) {
        return detail::dvr{l._node * r._node};
    }

    inline detail::dvr operator/(const detail::dvr &l, const detail::dvr &r) {
        return detail::dvr{l._node / r._node};
    }

    inline detail::dvr operator+(double l, const detail::dvr &r) {
        return detail::dvr{l + r._node};
    }

    inline detail::dvr operator-(double l, const detail::dvr &r) {
        return detail::dvr{l - r._node};
    }

    inline detail::dvr operator*(double l, const detail::dvr &r) {
        return detail::dvr{l * r._node};
    }

    inline detail::dvr operator/(double l, const detail::dvr &r) {
        return detail::dvr{l / r._node};
    }

    inline detail::dvr operator+(const detail::dvr &l, double r) {
        return detail::dvr{l._node + r};
    }

    inline detail::dvr operator-(const detail::dvr &l, double r) {
        return detail::dvr{l._node - r};
    }

    inline detail::dvr operator*(const detail::dvr &l, double r) {
        return detail::dvr{l._node * r};
    }

    inline detail::dvr operator/(const detail::dvr &l, double r) {
        return detail::dvr{l._node / r};
    }

    inline detail::dvr sin(const detail::dvr& d) {
        return detail::dvr{sin(d._node)};
    }

    inline detail::dvr cos(const detail::dvr& d) {
        return detail::dvr{cos(d._node)};
    }

    inline detail::dvr tan(const detail::dvr& d) {
        return detail::dvr{tan(d._node)};
    }

    inline detail::dvr pow(const detail::dvr& d, double n) {
        return detail::dvr{pow(d._node, n)};
    }

    inline detail::dvr pow(double a, const detail::dvr& d) {
        return detail::dvr{pow(a, d._node)};
    }

    inline detail::dvr sqrt(const detail::dvr& d) {
        return detail::dvr{sqrt(d._node)};
    }

    inline detail::dvr exp(const detail::dvr& d) {
        return detail::dvr{exp(d._node)};
    }

    inline detail::dvr ln(const detail::dvr& d) {
        return detail::dvr{ln(d._node)};
    }

    inline detail::dvr asin(const detail::dvr& d) {
        return detail::dvr{asin(d._node)};
    }

    inline detail::dvr acos(const detail::dvr& d) {
        return detail::dvr{acos(d._node)};
    }

    inline detail::dvr atan(const detail::dvr& d) {
        return detail::dvr{atan(d._node)};
    }

    //End Reverse Differential Variable Operator Overloading

    /*
     * FORWARD MODE OVERLOADING
     */
    static double value_at(const std::function<detail::dvf( detail::dvf)>& f,double x){
        auto u_dual = detail::dvf{ad::dual_number{x,1.0}};
        auto res = f(u_dual)._dual;
        return res._real_part;
    }

    static double gradient_at(const std::function<detail::dvf( detail::dvf)>& f,double x){
        auto u_dual = detail::dvf{ad::dual_number{x,1.0}};
        auto res = f(u_dual)._dual;
        return res._dual_part;
    }

    static std::pair<double,double> value_and_gradient_at(const std::function<detail::dvf( detail::dvf)>& f,double x){
        auto u_dual = detail::dvf{ad::dual_number{x,1.0}};
        auto res = f(u_dual)._dual;
        return {res._real_part,res._dual_part};
    }

    static double value_at(const std::function<detail::dvf(detail::dvf, detail::dvf)>& f,double x,double y){

        auto dfdx = f(detail::dvf{dual_number{x,1.0}},detail::dvf{dual_number{y,.0}})._dual;

        return dfdx._real_part;
    }

    static std::pair<double,double> gradient_at(const std::function<detail::dvf(detail::dvf, detail::dvf)>& f,double x,double y){

        auto dfdx = f(detail::dvf{dual_number{x,1.0}},detail::dvf{dual_number{y,.0}})._dual;
        auto dfdy = f(detail::dvf{dual_number{x,.0}},detail::dvf{dual_number{y,1.0}})._dual;

        assert(dfdx._real_part == dfdy._real_part);
        return {dfdx._dual_part,dfdy._dual_part};
    }

    static std::pair<double,std::pair<double,double>> value_and_gradient_at(const std::function<detail::dvf(detail::dvf, detail::dvf)>& f,double x,double y){

        auto dfdx = f(detail::dvf{dual_number{x,1.0}},detail::dvf{dual_number{y,.0}})._dual;
        auto dfdy = f(detail::dvf{dual_number{x,.0}},detail::dvf{dual_number{y,1.0}})._dual;

        assert(dfdx._real_part == dfdy._real_part);
        return {dfdx._real_part,{dfdx._dual_part,dfdy._dual_part}};
    }

    namespace detail{
        template<unsigned int n>
        using dvf_array = std::array<detail::dvf,n>;

        template<unsigned int n>
        using dvr_array = std::array<detail::dvr,n>;

        template<unsigned int n>
        using dvf_scalar_func_nd = std::function<detail::dvf(dvf_array<n>)>;

        template<unsigned int n>
        using dvr_scalar_func_nd = std::function<detail::dvr(dvr_array<n>)>;

        template<unsigned int n,unsigned int m>
        using dvf_md_func_nd = std::function<dvf_array<m>(dvf_array<n>)>;

        template<unsigned int n,unsigned int m>
        using dvr_md_func_nd = std::function<dvr_array<m>(dvr_array<n>)>;
    }
    //todo Really naive implementation...

    // value
    template<unsigned int n>
    static double value_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::array<double,n> args){
        return 0.0;
    }

    template<unsigned int n>
    static double value_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::initializer_list<double> args){

        std::array<double,n> arr{};
        int i = 0;
        for(double arg : args){
            arr[i++] = arg;
        }
        //return value_and_gradient_at<n>(f,arr);
    }

    //gradient for all
    template<unsigned int n>
    static std::array<double,n> gradient_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::array<double,n> args){
        std::array<double,n> grad{};

        std::array<detail::dvf,n> aargs{};

        std::transform(args.begin(),args.end(),aargs.begin(),[](const double arg){
            return detail::dvf{dual_number{arg,.0}};
        });

        double eval_result{};
        for(int i =0;i<args.size();++i){
            aargs[i]._dual._dual_part = 1;
            auto dual = f(aargs)._dual;
            aargs[i]._dual._dual_part = 0;
            assert(eval_result == 0 || eval_result == dual._real_part);
            eval_result = dual._real_part;
            grad[i] = dual._dual_part;
        }

        return {eval_result, grad};
    }

    template<unsigned int n>
    static std::array<double,n> gradient_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::initializer_list<double> args){

        std::array<double,n> arr{};
        int i = 0;
        for(double arg : args){
            arr[i++] = arg;
        }
        //return value_and_gradient_at<n>(f,arr);
    }

    //gradient for x_i
    template<unsigned int n>
    static double gradient_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::array<double,n> args,unsigned int grad_for){
        std::array<double,n> grad{};

        std::array<detail::dvf,n> aargs{};

        std::transform(args.begin(),args.end(),aargs.begin(),[](const double arg){
            return detail::dvf{dual_number{arg,.0}};
        });

        double eval_result{};
        for(int i =0;i<args.size();++i){
            aargs[i]._dual._dual_part = 1;
            auto dual = f(aargs)._dual;
            aargs[i]._dual._dual_part = 0;
            assert(eval_result == 0 || eval_result == dual._real_part);
            eval_result = dual._real_part;
            grad[i] = dual._dual_part;
        }

        return {eval_result, grad};
    }

    template<unsigned int n>
    static double gradient_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::initializer_list<double> args,unsigned int grad_for){
        std::array<double,n> grad{};

        std::array<detail::dvf,n> aargs{};

        std::transform(args.begin(),args.end(),aargs.begin(),[](const double arg){
            return detail::dvf{dual_number{arg,.0}};
        });

        double eval_result{};
        for(unsigned int i = 0U;i<args.size();++i){
            aargs[i]._dual._dual_part = 1;
            auto dual = f(aargs)._dual;
            aargs[i]._dual._dual_part = 0;
            assert(eval_result == 0 || eval_result == dual._real_part);
            eval_result = dual._real_part;
            grad[i] = dual._dual_part;
        }

        return {eval_result, grad};
    }

    //value and gradient for all

    template<unsigned int n>
    static std::pair<double,std::array<double,n>> value_and_gradient_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::array<double,n> args){
        std::array<double,n> grad{};

        std::array<detail::dvf,n> aargs{};

        std::transform(args.begin(),args.end(),aargs.begin(),[](const double arg){
            return detail::dvf{dual_number{arg,.0}};
        });

        double eval_result{};
        for(int i =0;i<args.size();++i){
            aargs[i]._dual._dual_part = 1;
            auto dual = f(aargs)._dual;
            aargs[i]._dual._dual_part = 0;
            assert(eval_result == 0 || eval_result == dual._real_part);
            eval_result = dual._real_part;
            grad[i] = dual._dual_part;
        }

        return {eval_result, grad};
    }

    template<unsigned int n>
    static std::pair<double,std::array<double,n>> value_and_gradient_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::initializer_list<double> args){

        std::array<double,n> arr{};
        int i = 0;
        for(double arg : args){
            arr[i++] = arg;
        }
        return value_and_gradient_at<n>(f,arr);
    }

    //value and gradient for x_i
    template<unsigned int n>
    static std::pair<double,double> value_and_gradient_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::array<double,n> args,unsigned int grad_for){

        std::array<detail::dvf,n> aargs{};

        std::transform(args.begin(),args.end(),aargs.begin(),[](const double arg){
            return detail::dvf{dual_number{arg,.0}};
        });

        aargs[grad_for]._dual._dual_part = 1;

        dual_number res = f(aargs)._dual;

        return {res._real_part, res._dual_part};
    }

    template<unsigned int n>
    static std::pair<double,double> value_and_gradient_at(
            const detail::dvf_scalar_func_nd<n>& f,
            std::initializer_list<double> args,unsigned int grad_for){

        std::array<detail::dvf,n> aargs{};

        std::transform(args.begin(),args.end(),aargs.begin(),[](const double arg){
            return detail::dvf{dual_number{arg,.0}};
        });

        aargs[grad_for]._dual._dual_part = 1;

        dual_number res = f(aargs)._dual;

        return {res._real_part, res._dual_part};
    }

    /*
     * REVERSE MODE OVERLOADING
     */
    static std::pair<double,double> value_and_gradient_at(const std::function<detail::dvr( detail::dvr)>& f,double x){
        detail::dvr _x {x,0U};
        auto z = f(_x);
        return {detail::eval(*z._node),detail::derivative_for(*z._node,0U)};
    }

    static std::pair<double,std::pair<double,double>> value_and_gradient_at(const std::function<detail::dvr( detail::dvr, detail::dvr)>& f,double x,double y){
        auto z = f({x,0U},{y,1U});
        auto out = detail::eval(*z._node);
        auto dfdx = detail::derivative_for(*z._node,0U);
        auto dfdy = detail::derivative_for(*z._node,1U);
        return {out,{dfdx,dfdy}};
    }

    //Really naive implementation...
    template<unsigned int n>
    static std::pair<double,std::array<double,n>> value_and_gradient_at(
            const std::function<detail::dvr(std::array<detail::dvr,n>)>& f,
            std::array<double,n> args){

        std::array<double,n> grad{};
        std::array<detail::dvr,n> aargs;

        unsigned int idx = 0U;

        for(const auto arg : args){
            aargs[idx] = {arg,idx};
            idx++;
        }

        auto z = f(aargs); //todo How to cache this?

        double eval_result = detail::eval(*z._node);

        for(int i =0;i<args.size();++i){
            grad[i] = detail::derivative_for(*z._node,i);
        }

        return {eval_result, grad};
    }

    template<unsigned int n>
    static std::pair<double,std::array<double,n>> value_and_gradient_at(
            const std::function<detail::dvr(std::array<detail::dvr,n>)>& f,
            std::initializer_list<double> args){

        std::array<double,n> arr{};
        int i = 0;
        for(double arg : args){
            arr[i++] = arg;
        }
        return gradient_at<n>(f,arr);
    }

}
#endif //AUTODIFF_AUTO_DIFF_H
