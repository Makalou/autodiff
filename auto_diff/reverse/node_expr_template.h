//
// Created by 王泽远 on 2023/1/12.
//

#ifndef AUTODIFF_NODE_EXPR_TEMPLATE_H
#define AUTODIFF_NODE_EXPR_TEMPLATE_H

#include <cmath>

namespace ad::expr_template{

        template <typename E>
        struct expr_node{};

        struct variable_node : expr_node<variable_node>{
            double _val{};
            unsigned int _idx;

            variable_node(double val,unsigned int idx) : _val(val), _idx(idx){};

            double eval() const{
                return _val;
            }

            double derivative_for(unsigned int idx) const{
                return idx == _idx?1:0;
            }
        };

        template<class E1, class E2>
        struct add_node : expr_node<add_node<E1,E2>>{
            const E1 & _l;
            const E2 & _r;

            add_node(E1 const& l, E2 const& r) : _l(l), _r(r) { }

            double eval() const{
                return _l.eval() + _r.eval();
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E2>
        struct add_node<double,E2> : expr_node<add_node<double,E2>>{
            double _l;
            const E2 & _r;

            add_node(double l, E2 const& r) : _l(l), _r(r) { }

            double eval() const{
                return _l + _r.eval();
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E1>
        struct add_node<E1,double> : expr_node<add_node<E1,double>>{
            const E1 & _l;
            double _r;

            add_node(E1 const& l, double r) : _l(l), _r(r) { }

            double eval() const{
                return _l.eval() + _r;
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E1, class E2>
        struct minus_node : expr_node<minus_node<E1,E2>>{
            const E1 & _l;
            const E2 & _r;

            minus_node(E1 const& l, E2 const& r) : _l(l), _r(r) { }

            double eval() const{
                return _l.eval() - _r.eval();
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E1>
        struct minus_node<E1,double> : expr_node<minus_node<E1,double>>{
            const E1 & _l;
            double _r;

            minus_node(E1 const& l, double r) : _l(l), _r(r) { }

            double eval() const{
                return _l.eval() - _r;
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E2>
        struct minus_node<double,E2> : expr_node<minus_node<double,E2>>{
            double _l;
            const E2 & _r;

            minus_node(double l, E2 const& r) : _l(l), _r(r) { }

            double eval() const{
                return _l - _r.eval();
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E1, class E2>
        struct multiply_node : expr_node<multiply_node<E1,E2>>{
            const E1 & _l;
            const E2 & _r;

            multiply_node(E1 const& l, E2 const& r) : _l(l), _r(r) { }

            double eval() const{
                auto l = _l.eval();
                auto r = _r.eval();
                auto res = l*r;
                return res;
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E1>
        struct multiply_node<E1,double> : expr_node<multiply_node<E1,double>>{
            const E1 & _l;
            double _r;

            multiply_node(E1 const& l, double r) : _l(l), _r(r) { }

            double eval() const{
                return _l.eval() * _r;
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E2>
        struct multiply_node<double,E2> : expr_node<multiply_node<double,E2>>{
            double _l;
            const E2 & _r;

            multiply_node(double l, E2 const& r) : _l(l), _r(r) { }

            double eval() const{
                return _l * _r.eval();
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E1, class E2>
        struct divide_node : expr_node<divide_node<E1,E2>>{
            const E1 & _l;
            const E2 & _r;

            divide_node(E1 const& l, E2 const& r) : _l(l), _r(r) { }

            double eval() const {
                return _l.eval() / _r.eval();
            }

            double derivative_for(unsigned int idx) const {

            }
        };

        template<class E1>
        struct divide_node<E1,double> : expr_node<divide_node<E1,double>>{
            const E1 & _l;
            double _r;

            divide_node(E1 const& l, double r) : _l(l), _r(r) { }

            double eval() {
                return _l.eval() / _r;
            }

            double derivative_for(unsigned int idx){

            }
        };

        template<class E2>
        struct divide_node<double,E2> : expr_node<divide_node<double,E2>>{
            double _l;
            const E2 & _r;

            divide_node(double l, E2 const& r) : _l(l), _r(r) { }

            double eval() const{
                return _l / _r.eval();
            }

            double derivative_for(unsigned int idx)const{

            }
        };

        template<class E>
        struct sin_node : expr_node<sin_node<E>>{
            const E & _l;

            explicit sin_node(E const& l) : _l(l){ }

            double eval() const{
                return std::sin(_l.eval());
            }

            double derivative_for(unsigned int idx)const{

            }
        };

        template<class E>
        struct cos_node : expr_node<cos_node<E>>{
            const E & _l;

            explicit cos_node(E const& l) : _l(l){ }

            double eval() const{
                return std::cos(_l.eval());
            }

            double derivative_for(unsigned int idx)const{

            }
        };

        template<class E>
        struct tan_node : expr_node<tan_node<E>>{
            const E & _l;

            explicit tan_node(E const& l) : _l(l){ }

            double eval() const{
                return std::tan(_l.eval());
            }

            double derivative_for(unsigned int idx)const{

            }
        };

        template<class E>
        struct pow_node : expr_node<pow_node<E>>{
            const E & _l;
            double _pow{};

            pow_node(E const& l, double pow) : _l(l),_pow(pow){}

            double eval()const {
                return std::pow(_l.eval(),_pow);
            }

            double derivative_for(unsigned int idx)const{

            }
        };

        template<class E>
        struct sqrt_node : expr_node<sqrt_node<E>>{
            const E & _l;

            explicit sqrt_node(E const& l) : _l(l){ }

            double eval()const {
                return std::sqrt(_l.eval());
            }

            double derivative_for(unsigned int idx)const{

            }
        };

        template<class E>
        struct exp_node : expr_node<exp_node<E>>{
            const E & _l;
            bool is_base_e{false};
            double _base{};

            explicit exp_node(E const& l) : _l(l),is_base_e(true){ }

            exp_node(double base, E const& l):_l(l),_base(base){}

            double eval() const{
                return is_base_e?std::exp(_l.eval()):std::pow(_base,_l.eval());
            }

            double derivative_for(unsigned int idx)const{

            }
        };

        template<class E>
        struct ln_node : expr_node<ln_node<E>>{
            const E & _l;

            explicit ln_node(E const& l) : _l(l){ }

            double eval() const{
                return std::log(_l.eval());
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E>
        struct asin_node : expr_node<asin_node<E>>{
            const E & _l;

            explicit asin_node(E const& l) : _l(l){ }

            double eval() const{
                return std::asin(_l.eval());
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E>
        struct acos_node : expr_node<acos_node<E>>{
            const E & _l;

            explicit acos_node(E const& l) : _l(l){ }

            double eval() const{
                return std::acos(_l.eval());
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template<class E>
        struct atan_node : expr_node<atan_node<E>>{
            const E & _l;

            explicit atan_node(E const& l) : _l(l){ }

            double eval() const{
                return std::atan(_l.eval());
            }

            double derivative_for(unsigned int idx) const{

            }
        };

        template <typename E1, typename E2>
        inline auto operator+(expr_node<E1> const& u, expr_node<E2> const& v) {
            return add_node<E1, E2>(static_cast<const E1&>(u), static_cast<const E2&>(v));
        }

        template <typename E1>
        inline auto operator+(expr_node<E1> const& u, double v) {
            return add_node<E1, double>(static_cast<const E1&>(u), v);
        }

        template <typename E2>
        inline auto operator+(double u, expr_node<E2> const& v) {
            return add_node<double, E2>(u, static_cast<const E2&>(v));
        }

        template <typename E1, typename E2>
        inline auto operator-(expr_node<E1> const& u, expr_node<E2> const& v) {
            return minus_node<E1, E2>(static_cast<const E1&>(u), static_cast<const E2&>(v));
        }

        template <typename E1>
        inline auto operator-(expr_node<E1> const& u, double v) {
            return minus_node<E1,double>(static_cast<const E1&>(u), v);
        }

        template <typename E2>
        inline auto operator-(double u, expr_node<E2> const& v) {
            return minus_node<double, E2>(u, static_cast<const E2&>(v));
        }

        template <typename E1, typename E2>
        inline auto operator*(expr_node<E1> const& u, expr_node<E2> const& v) {
            return multiply_node<E1, E2>(static_cast<const E1&>(u), static_cast<const E2&>(v));
        }

        template <typename E1>
        inline auto operator*(expr_node<E1> const& u, double v) {
            return multiply_node<E1, double>(static_cast<const E1&>(u),v);
        }

        template <typename E2>
        inline auto operator*(double u, expr_node<E2> const& v) {
            return multiply_node<double, E2>(u, static_cast<const E2&>(v));
        }

        template <typename E1, typename E2>
        inline auto operator/(expr_node<E1> const& u, expr_node<E2> const& v) {
            return divide_node<E1, E2>(static_cast<const E1&>(u), static_cast<const E2&>(v));
        }

        template <typename E1>
        inline auto operator/(expr_node<E1> const& u, double v) {
            return divide_node<E1, double>(static_cast<const E1&>(u), v);
        }

        template <typename E2>
        inline auto operator/(double u, expr_node<E2> const& v) {
            return divide_node<double, E2>(u, static_cast<const E2&>(v));
        }

        template<typename E>
        inline auto sin(expr_node<E> const& v){
            return sin_node<E>(static_cast<const E&>(v));
        }

        template<typename E>
        inline auto cos(expr_node<E> const& v){
            return cos_node<E>(static_cast<const E&>(v));
        }

        template<typename E>
        inline auto tan(expr_node<E> const& v){
            return tan_node<E>(static_cast<const E&>(v));
        }

        template<typename E>
        inline auto pow(expr_node<E> const& v,double n){
            return pow_node<E>(static_cast<const E&>(v),n);
        }

        template<typename E>
        inline auto pow(double  a,expr_node<E> const& v){
            return exp_node<E>(a,static_cast<const E&>(v));
        }

        template<typename E>
        inline auto sqrt(expr_node<E> const& v){
            return sqrt_node<E>(static_cast<const E&>(v));
        }

        template<typename E>
        inline auto exp(expr_node<E> const& v){
            return exp_node<E>(static_cast<const E&>(v));
        }

        template<typename E>
        inline auto ln(expr_node<E> const& v){
            return ln_node<E>(static_cast<const E&>(v));
        }

        template<typename E>
        inline auto asin(expr_node<E> const& v){
            return asin_node<E>(static_cast<const E&>(v));
        }

        template<typename E>
        inline auto acos(expr_node<E> const& v){
            return acos_node<E>(static_cast<const E&>(v));
        }

        template<typename E>
        inline auto atan(expr_node<E> const& v){
            return atan_node<E>(static_cast<const E&>(v));
        }
    }

#endif //AUTODIFF_NODE_EXPR_TEMPLATE_H
