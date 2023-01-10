//
// Created by 王泽远 on 2023/1/9.
//

#ifndef AUTODIFF_NODE_H
#define AUTODIFF_NODE_H

#include <utility>
#include <memory>

namespace ad{

    struct differential_node;

    namespace detail{
        using nsp = std::shared_ptr<differential_node>;
    }

    using uint = unsigned int;

    struct differential_node{
    public:
        detail::nsp _l;
        detail::nsp _r;

        double _eval_cache{};
        bool _eval_cache_init{false};

        double eval(){

            if(!_eval_cache_init){
                _eval_cache = do_eval();
                _eval_cache_init = true;
            }

            return _eval_cache;
        }

        virtual double do_eval()=0;

        virtual bool depend_on(uint idx)= 0;

        virtual double derivative_for(uint idx)=0;

        friend
        detail::nsp operator+(const detail::nsp& l,const detail::nsp& r);
        friend
        detail::nsp operator+(double l,const detail::nsp& r);
        friend
        detail::nsp operator+(const detail::nsp& l,double r);
    };

    struct constant_node : differential_node{
        double _val{};
        explicit constant_node(double c):_val(c){

        }

        double do_eval() override{
            return _val;
        }

        bool depend_on(uint idx) override{
            return false;
        }

        double derivative_for(uint idx) override{
            return 0;
        }
    };

    struct variable_node : differential_node{
        double _val{};
        uint _idx;
        explicit variable_node(double v,uint idx):_val(v),_idx(idx){

        }

        double do_eval() override{
            return _val;
        }

        bool depend_on(uint idx) override{
            return idx == _idx;
        }

        double derivative_for(uint idx) override{
            return idx == _idx? 1.0 : 0.0;
        }
    };

    struct add_node : differential_node{
        add_node(const detail::nsp& l,const detail::nsp& r){
            _l = l;
            _r = r;
        }

        double do_eval() override{
            return _l->do_eval() + _r->do_eval();
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx) || _r->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            double temp = 0;
            if(_l->depend_on(idx)){
                temp += _l->derivative_for(idx);
            }
            if(_r->depend_on(idx)){
                temp += _r->derivative_for(idx);
            }
            return temp;
        }
    };

    struct minus_node : differential_node{
        minus_node(const detail::nsp& l,const detail::nsp& r){
            _l = l;
            _r = r;
        }

        double do_eval() override{
            return _l->do_eval() - _r->do_eval();
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx) || _r->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            double temp = 0;
            if(_l->depend_on(idx)){
                temp += _l->derivative_for(idx);
            }
            if(_r->depend_on(idx)){
                temp -= _r->derivative_for(idx);
            }
            return temp;
        }
    };

    struct multiply_node : differential_node{
        multiply_node(const detail::nsp& l,const detail::nsp& r){
            _l = l;
            _r = r;
        }

        double do_eval() override{
            return _l->do_eval() * _r->do_eval();
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx) || _r->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            auto dfdx_mul_gx = _l->depend_on(idx)?_l->derivative_for(idx)*_r->eval():0;
            auto dgdx_mul_fx = _r->depend_on(idx)?_r->derivative_for(idx)*_l->eval():0;
            return dfdx_mul_gx + dgdx_mul_fx;
        }
    };

    struct divide_node : differential_node{
        divide_node(const detail::nsp& l,const detail::nsp& r){
            _l = l;
            _r = r;
        }

        double do_eval() override{
            return _l->do_eval() / _r->do_eval();
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx) || _r->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            auto gx = _r->eval();
            auto gx2 = gx*gx;
            auto dfdx_mul_gx = _l->depend_on(idx)?_l->derivative_for(idx)*_r->eval():0;
            auto dgdx_mul_fx = _r->depend_on(idx)?_r->derivative_for(idx)*_l->eval():0;
            return (dfdx_mul_gx-dgdx_mul_fx)/gx2;
        }
    };

    struct sin_node : differential_node{
        explicit sin_node(const detail::nsp& l){
            _l = l;
        }

        double do_eval() override{
            return std::sin(_l->do_eval());
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            auto cos = std::cos(_l->eval());
            auto dx =  _l->derivative_for(idx);
            return cos*dx;
        }
    };

    struct cos_node : differential_node{
        explicit cos_node(const detail::nsp& l){
            _l = l;
        }

        double do_eval() override{
            return std::cos(_l->do_eval());
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            return -std::sin(_l->eval())*_l->derivative_for(idx);
        }
    };

    struct tan_node : differential_node{
        explicit tan_node(const detail::nsp& l){
            _l = l;
        }

        double do_eval() override{
            return std::tan(_l->do_eval());
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            auto sec = 1/std::cos(_l->eval());
            auto sec2 = sec*sec;
            return sec2*_l->derivative_for(idx);
        }
    };

    struct pow_node : differential_node{
        double _pow;
        pow_node(const detail::nsp& l,double pow){
            _l = l;
            _pow = pow;
        }

        double do_eval() override{
            return std::pow(_l->do_eval(), _pow);
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            return _pow*std::pow(_l->eval(),_pow-1)*_l->derivative_for(idx);
        }
    };

    struct sqrt_node : differential_node{
        explicit sqrt_node(const detail::nsp& l){
            _l = l;
        }

        double do_eval() override{
            return std::sqrt(_l->do_eval());
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            return 0.5/std::sqrt(_l->eval())*_l->derivative_for(idx);
        }
    };

    struct exp_node : differential_node{
        double _base;
        bool is_base_e;

        explicit exp_node(const detail::nsp& l){
            _l = l;
            is_base_e = true;
        }

        exp_node(double base,const detail::nsp& l){
            _l = l;
            _base = base;
            is_base_e = false;
        }

        double do_eval() override{
            return is_base_e ? std::exp(_l->do_eval()) : std::pow(_base, _l->do_eval());
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            if(is_base_e){
                return std::exp(_l->eval())*_l->derivative_for(idx);
            }else{
                return std::log(_base)*std::pow(_base,_l->eval())*_l->derivative_for(idx);
            }
        }
    };

    struct ln_node : differential_node{
        explicit ln_node(const detail::nsp& l){
            _l = l;
        }

        double do_eval() override{
            return std::log(_l->do_eval());
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            return 1.0/_l->eval()*_l->derivative_for(idx);
        }
    };

    struct asin_node : differential_node{
        explicit asin_node(const detail::nsp& l){
            _l = l;
        }

        double do_eval() override{
            return std::asin(_l->do_eval());
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            auto x = _l->eval();
            auto dx = 1/(std::sqrt(1-x*x));
            return dx*_l->derivative_for(idx);
        }
    };

    struct acos_node : differential_node{
        explicit acos_node(const detail::nsp& l){
            _l = l;
        }

        double do_eval() override{
            return std::acos(_l->do_eval());
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            auto x = _l->eval();
            auto dx = -1/(std::sqrt(1-x*x));
            return dx*_l->derivative_for(idx);
        }
    };

    struct atan_node : differential_node{
        explicit atan_node(const detail::nsp& l){
            _l = l;
        }

        double do_eval() override{
            return std::atan(_l->do_eval());
        }

        bool depend_on(uint idx) override{
            return _l->depend_on(idx);
        }

        double derivative_for(uint idx) override{
            auto x = _l->eval();;
            auto dx = 1/(1+x*x);
            return dx*_l->derivative_for(idx);
        }
    };

    inline detail::nsp operator+(const  detail::nsp &l, const  detail::nsp &r) {
        return std::make_shared<add_node>(l,r);
    }

    inline detail::nsp operator-(const  detail::nsp &l, const  detail::nsp &r) {
        return std::make_shared<minus_node>(l,r);
    }

    inline detail::nsp operator*(const  detail::nsp &l, const  detail::nsp &r) {
        return std::make_shared<multiply_node>(l,r);
    }

    inline detail::nsp operator/(const  detail::nsp &l, const  detail::nsp &r) {
        return std::make_shared<divide_node>(l,r);
    }

    inline detail::nsp operator+(double l, const  detail::nsp &r) {
        return std::make_shared<add_node>(std::make_shared<constant_node>(l),r);
    }

    inline detail::nsp operator-(double l, const  detail::nsp&r) {
        return std::make_shared<minus_node>(std::make_shared<constant_node>(l),r);
    }

    inline detail::nsp operator*(double l, const  detail::nsp &r) {
        return std::make_shared<multiply_node>(std::make_shared<constant_node>(l),r);
    }

    inline detail::nsp operator/(double l, const  detail::nsp &r) {
        return std::make_shared<divide_node>(std::make_shared<constant_node>(l),r);
    }

    inline detail::nsp operator+(const  detail::nsp &l, double r){
        return std::make_shared<add_node>(l,std::make_shared<constant_node>(r));
    }

    inline detail::nsp operator-(const  detail::nsp &l, double r){
        return std::make_shared<minus_node>(l,std::make_shared<constant_node>(r));
    }

    inline detail::nsp operator*(const  detail::nsp &l, double r){
        return std::make_shared<multiply_node>(l,std::make_shared<constant_node>(r));
    }

    inline detail::nsp operator/(const  detail::nsp &l, double r){
        return std::make_shared<divide_node>(l,std::make_shared<constant_node>(r));
    }

    inline detail::nsp sin( const detail::nsp& d) {
        return std::make_shared<sin_node>(d);
    }

    inline detail::nsp cos( const detail::nsp& d) {
        return std::make_shared<cos_node>(d);
    }

    inline detail::nsp tan( const detail::nsp& d) {
        return std::make_shared<tan_node>(d);
    }

    inline detail::nsp pow( const detail::nsp& d, double n){
        return std::make_shared<pow_node>(d,n);
    }

    inline detail::nsp pow(double a,  const detail::nsp& d){
        return std::make_shared<exp_node>(a,d);
    }

    inline detail::nsp sqrt( const detail::nsp& d){
        return std::make_shared<sqrt_node>(d);
    }

    inline detail::nsp exp( const detail::nsp& d){
        return std::make_shared<exp_node>(d);
    }

    inline detail::nsp ln( const detail::nsp& d){
        return std::make_shared<ln_node>(d);
    }

    inline detail::nsp asin( const detail::nsp& d){
        return std::make_shared<asin_node>(d);
    }

    inline detail::nsp acos( const detail::nsp& d){
        return std::make_shared<acos_node>(d);
    }

    inline detail::nsp atan( const detail::nsp& d){
        return std::make_shared<atan_node>(d);
    }
}
#endif //AUTODIFF_NODE_H
