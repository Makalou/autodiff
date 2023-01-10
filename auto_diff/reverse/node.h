//
// Created by 王泽远 on 2023/1/9.
//

#ifndef AUTODIFF_NODE_H
#define AUTODIFF_NODE_H

#include <utility>
#include <memory>
#include <variant>

namespace ad{

    //struct differential_node;

    using uint = unsigned int;

    /*
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
     */

    struct constant_node;

    struct variable_node;

    struct add_node;

    struct minus_node;

    struct multiply_node;

    struct divide_node;

    struct sin_node;

    struct cos_node;

    struct tan_node;

    struct pow_node;

    struct sqrt_node;

    struct exp_node;

    struct ln_node;

    struct asin_node;

    struct acos_node;

    struct atan_node;

    using differential_node = std::variant<constant_node,variable_node,add_node,minus_node,multiply_node,divide_node,sin_node
                                            ,cos_node,tan_node,pow_node,sqrt_node,exp_node,ln_node,asin_node,acos_node,atan_node>;

    namespace detail{
        using nsp = std::shared_ptr<differential_node>;

        inline nsp make_differential_node(differential_node&& node){
            return std::make_shared<differential_node>(node);
        }
    }

    struct constant_node{
        double val{};
    };

    struct variable_node{
        double val{};
        unsigned int idx;
    };

    struct add_node{
        detail::nsp _l;
        detail::nsp _r;
    };

    struct minus_node{
        detail::nsp _l;
        detail::nsp _r;

    };

    struct multiply_node{
        detail::nsp _l;
        detail::nsp _r;
    };

    struct divide_node{
        detail::nsp _l;
        detail::nsp _r;
    };

    struct sin_node{
        detail::nsp _l;

    };

    struct cos_node{
        detail::nsp _l;
    };

    struct tan_node{
        detail::nsp _l;
    };

    struct pow_node{
        double _pow;
        detail::nsp _l;

        pow_node(const detail::nsp& l,double pow){
            _l = l;
            _pow = pow;
        }
    };

    struct sqrt_node{
        detail::nsp _l;
    };

    struct exp_node{
        double _base;
        bool is_base_e;

        detail::nsp _l;

        explicit exp_node(const detail::nsp& l){
            _l = l;
            is_base_e = true;
        }
        exp_node(double base,const detail::nsp& l){
            _l = l;
            _base = base;
            is_base_e = false;
        }
    };

    struct ln_node{
        detail::nsp _l;
    };

    struct asin_node{
        detail::nsp _l;
    };

    struct acos_node{
        detail::nsp _l;
    };

    struct atan_node{
        detail::nsp _l;
    };

    namespace detail{
        // helper type for the visitor #4
        template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
        // explicit deduction guide (not needed as of C++20)
        template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

        static double eval(const differential_node& node){
            return std::visit(overloaded{
                    [](const auto& node){
                        return 0.0;
                        },
                    [](const constant_node& node){
                        return node.val;
                    },
                    [](const variable_node& node){
                        return node.val;
                    },
                    [](const add_node& node){
                        return eval(*node._l) + eval(*node._r);
                        },
                    [](const minus_node& node){
                        return eval(*node._l) - eval(*node._r);
                        },
                    [](const multiply_node& node){
                        return eval(*node._l) * eval(*node._r);
                        },
                    [](const divide_node& node){
                        return eval(*node._l) / eval(*node._r);
                        },
                    [](const sin_node& node){
                        return std::sin(eval(*node._l));
                        },
                    [](const cos_node& node){
                        return std::cos(eval(*node._l));
                        },
                    [](const tan_node& node){
                        return std::tan(eval(*node._l));
                        },
                    [](const pow_node& node){
                        return std::pow(eval(*node._l),node._pow);
                        },
                    [](const sqrt_node& node){
                        return std::sqrt(eval(*node._l));
                        },
                    [](const exp_node& node) {
                        return node.is_base_e ? std::exp(eval(*node._l)) : std::pow(node._base, eval(*node._l));
                        },
                    [](const ln_node& node) {
                        return std::log(eval(*node._l));
                        },
                    [](const asin_node& node){
                        return std::asin(eval(*node._l));
                        },
                    [](const acos_node& node){
                        return std::acos(eval(*node._l));
                        },
                    [](const atan_node& node){
                        return std::atan(eval(*node._l));
                        }
            },node);
        }

        static bool depend_on(const differential_node& node,unsigned int idx){
            return std::visit(overloaded{
                    [](const auto& node){
                        return false;
                    },
                    [](const constant_node& node){
                        return false;
                    },
                    [idx](const variable_node& node){
                        return node.idx == idx;
                    },
                    [idx](const add_node& node){
                        return depend_on(*node._l,idx) || depend_on(*node._r,idx);
                    },
                    [idx](const minus_node& node){
                        return depend_on(*node._l,idx) || depend_on(*node._r,idx);
                    },
                    [idx](const multiply_node& node){
                        return depend_on(*node._l,idx) || depend_on(*node._r,idx);
                    },
                    [idx](const divide_node& node){
                        return depend_on(*node._l,idx) || depend_on(*node._r,idx);
                    },
                    [idx](const sin_node& node){
                        return depend_on(*node._l,idx);
                    },
                    [idx](const cos_node& node){
                        return depend_on(*node._l,idx);
                    },
                    [idx](const tan_node& node){
                        return depend_on(*node._l,idx);
                    },
                    [idx](const pow_node& node){
                        return depend_on(*node._l,idx);
                    },
                    [idx](const sqrt_node& node){
                        return depend_on(*node._l,idx);
                    },
                    [idx](const exp_node& node) {
                        return depend_on(*node._l,idx);
                    },
                    [idx](const ln_node& node) {
                        return depend_on(*node._l,idx);
                    },
                    [idx](const asin_node& node){
                        return depend_on(*node._l,idx);
                    },
                    [idx](const acos_node& node){
                        return depend_on(*node._l,idx);
                    },
                    [idx](const atan_node& node){
                        return depend_on(*node._l,idx);
                    }
            },node);
        }

        static double derivative_for(const differential_node& node, unsigned int idx){
            return std::visit(overloaded{
                    [idx](const auto& node){
                        return 0.0;
                    },
                    [idx](const constant_node& node){
                        return 0.0;
                    },
                    [idx](const variable_node& node){
                        return idx == node.idx? 1.0 : 0.0;
                    },
                    [idx](const add_node& node){
                        double temp = 0;
                        if(depend_on(*node._l,idx)){
                            temp += derivative_for(*node._l,idx);
                        }
                        if(depend_on(*node._r,idx)){
                            temp += derivative_for(*node._r,idx);
                        }
                        return temp;
                    },
                    [idx](const minus_node& node){
                        double temp = 0;
                        if(depend_on(*node._l,idx)){
                            temp += derivative_for(*node._l,idx);
                        }
                        if(depend_on(*node._r,idx)){
                            temp -= derivative_for(*node._r,idx);
                        }
                        return temp;
                    },
                    [idx](const multiply_node& node){
                        auto dfdx_mul_gx = depend_on(*node._l,idx)?derivative_for(*node._l,idx)* eval(*node._r):0;
                        auto dgdx_mul_fx = depend_on(*node._r,idx)?derivative_for(*node._r,idx)* eval(*node._l):0;
                        return dfdx_mul_gx + dgdx_mul_fx;
                    },
                    [idx](const divide_node& node){
                        auto gx = eval(*node._r);
                        auto gx2 = gx*gx;
                        auto dfdx_mul_gx = depend_on(*node._l,idx)?derivative_for(*node._l,idx)* eval(*node._r):0;
                        auto dgdx_mul_fx = depend_on(*node._r,idx)?derivative_for(*node._r,idx)* eval(*node._l):0;
                        return (dfdx_mul_gx-dgdx_mul_fx)/gx2;
                    },
                    [idx](const sin_node& node){
                        auto cos = std::cos(eval(*node._l));
                        auto dx =  derivative_for(*node._l,idx);
                        return cos*dx;
                    },
                    [idx](const cos_node& node){
                        return -std::sin(eval(*node._l))*derivative_for(*node._l,idx);
                    },
                    [idx](const tan_node& node){
                        auto sec = 1/std::cos(eval(*node._l));
                        auto sec2 = sec*sec;
                        return sec2*derivative_for(*node._l,idx);
                    },
                    [idx](const pow_node& node){
                        return node._pow*std::pow(eval(*node._l),node._pow-1)*derivative_for(*node._l,idx);
                    },
                    [idx](const sqrt_node& node){
                        return 0.5/std::sqrt(eval(*node._l))*derivative_for(*node._l,idx);
                    },
                    [idx](const exp_node& node) {
                        return node.is_base_e? std::exp(eval(*node._l))*derivative_for(*node._l,idx):std::log(node._base)*std::pow(node._base,eval(*node._l))*derivative_for(*node._l,idx);;
                    },
                    [idx](const ln_node& node) {
                        return 1.0/eval(*node._l)*derivative_for(*node._l,idx);
                    },
                    [idx](const asin_node& node){
                        auto x = eval(*node._l);
                        auto dx = 1/(std::sqrt(1-x*x));
                        return dx*derivative_for(*node._l,idx);
                    },
                    [idx](const acos_node& node){
                        auto x = eval(*node._l);
                        auto dx = -1/(std::sqrt(1-x*x));
                        return dx*derivative_for(*node._l,idx);
                    },
                    [idx](const atan_node& node){
                        auto x = eval(*node._l);;
                        auto dx = 1/(1+x*x);
                        return dx*derivative_for(*node._l,idx);
                    }
            },node);
        }
    }

    inline detail::nsp operator+(const  detail::nsp &l, const  detail::nsp &r) {
        return detail::make_differential_node(add_node{l,r});
    }

    inline detail::nsp operator-(const  detail::nsp &l, const  detail::nsp &r) {
        return detail::make_differential_node(minus_node{l,r});
    }

    inline detail::nsp operator*(const  detail::nsp &l, const  detail::nsp &r) {
        return detail::make_differential_node(multiply_node{l,r});
    }

    inline detail::nsp operator/(const  detail::nsp &l, const  detail::nsp &r) {
        return detail::make_differential_node(divide_node{l,r});
    }

    inline detail::nsp operator+(double l, const  detail::nsp &r) {
        return detail::make_differential_node(add_node{detail::make_differential_node({constant_node{l}}),r});
    }

    inline detail::nsp operator-(double l, const  detail::nsp&r) {
        return detail::make_differential_node(minus_node{detail::make_differential_node(constant_node{l}),r});
    }

    inline detail::nsp operator*(double l, const  detail::nsp &r) {
        return detail::make_differential_node(multiply_node{detail::make_differential_node(constant_node{l}),r});
    }

    inline detail::nsp operator/(double l, const  detail::nsp &r) {
        return detail::make_differential_node(divide_node{detail::make_differential_node(constant_node{l}),r});
    }

    inline detail::nsp operator+(const  detail::nsp &l, double r){
        return detail::make_differential_node(add_node{l, detail::make_differential_node(constant_node{r})});
    }

    inline detail::nsp operator-(const  detail::nsp &l, double r){
        return detail::make_differential_node(minus_node{l, detail::make_differential_node(constant_node{r})});
    }

    inline detail::nsp operator*(const  detail::nsp &l, double r){
        return detail::make_differential_node(multiply_node{l, detail::make_differential_node(constant_node{r})});
    }

    inline detail::nsp operator/(const  detail::nsp &l, double r){
        return detail::make_differential_node(divide_node{l, detail::make_differential_node(constant_node{r})});
    }

    inline detail::nsp sin( const detail::nsp& d) {
        return detail::make_differential_node(sin_node{d});
    }

    inline detail::nsp cos( const detail::nsp& d) {
        return detail::make_differential_node(cos_node{d});
    }

    inline detail::nsp tan( const detail::nsp& d) {
        return detail::make_differential_node(tan_node{d});
    }

    inline detail::nsp pow( const detail::nsp& d, double n){
        return detail::make_differential_node(pow_node{d,n});
    }

    inline detail::nsp pow(double a,  const detail::nsp& d){
        return detail::make_differential_node(exp_node{a,d});
    }

    inline detail::nsp sqrt( const detail::nsp& d){
        return detail::make_differential_node(sqrt_node{d});
    }

    inline detail::nsp exp( const detail::nsp& d){
        return detail::make_differential_node(exp_node{d});
    }

    inline detail::nsp ln( const detail::nsp& d){
        return detail::make_differential_node(ln_node{d});
    }

    inline detail::nsp asin( const detail::nsp& d){
        return detail::make_differential_node(asin_node{d});
    }

    inline detail::nsp acos( const detail::nsp& d){
        return detail::make_differential_node(acos_node{d});
    }

    inline detail::nsp atan( const detail::nsp& d){
        return detail::make_differential_node(atan_node{d});
    }
}
#endif //AUTODIFF_NODE_H
