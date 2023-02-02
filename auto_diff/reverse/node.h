//
// Created by 王泽远 on 2023/1/9.
//

#ifndef AUTODIFF_NODE_H
#define AUTODIFF_NODE_H

#include <utility>
#include <memory>
#include <variant>

namespace ad{

    using uint = unsigned int;

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
        using nup = std::unique_ptr<differential_node>;

        inline nsp make_differential_node_shared(differential_node && node){
            return std::make_shared<differential_node>(node);
        }

        inline nup make_differential_node_unique(differential_node&& node){
            return std::make_unique<differential_node>(node);
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
        double constant{};
        add_node(detail::nsp  l, detail::nsp  r) : _l(std::move(l)),_r(std::move(r)){};
        add_node(double l, detail::nsp  r) : constant(l),_r(std::move(r)){};
        add_node(detail::nsp l, double  r) : _l(std::move(l)),constant(r){};
    };

    struct minus_node{
        detail::nsp _l;
        detail::nsp _r;
        double constant{};
        minus_node(detail::nsp  l, detail::nsp  r) : _l(std::move(l)),_r(std::move(r)){};
        minus_node(double  l, detail::nsp  r) : constant(l),_r(std::move(r)){};
        minus_node(detail::nsp l, double  r) : _l(std::move(l)),constant(r){};
    };

    struct multiply_node{
        detail::nsp _l;
        detail::nsp _r;
        double constant{};
        multiply_node(detail::nsp l, detail::nsp r):_l(std::move(l)),_r(std::move(r)){};
        multiply_node(double l, detail::nsp r):constant(l),_r(std::move(r)){};
        multiply_node(detail::nsp l, double r):_l(std::move(l)),constant(r){};
    };

    struct divide_node{
        detail::nsp _l;
        detail::nsp _r;
        double constant{};
        divide_node(detail::nsp l, detail::nsp r):_l(std::move(l)),_r(std::move(r)){};
        divide_node(double l, detail::nsp r):constant(l),_r(std::move(r)){};
        divide_node(detail::nsp l, double r):_l(std::move(l)),constant(r){};
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
        double _base{};
        bool is_base_e{false};

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
                    [](const constant_node& node){
                        return node.val;
                    },
                    [](const variable_node& node){
                        return node.val;
                    },
                    [](const add_node& node) {
                        if (node._l && node._r) {
                            return eval(*node._l) + eval(*node._r);
                        } else if (node._l) {
                            return eval(*node._l) + node.constant;
                        } else if (node._r) {
                            return node.constant + eval(*node._l);
                        }else{
                            return 0.0;
                        }
                    },
                    [](const minus_node& node){
                        if (node._l && node._r) {
                            return eval(*node._l) - eval(*node._r);
                        } else if (node._l) {
                            return eval(*node._l) - node.constant;
                        } else if (node._r) {
                            return node.constant - eval(*node._l);
                        }else {
                            return 0.0;
                        }
                        },
                    [](const multiply_node& node){
                        if (node._l && node._r) {
                            return eval(*node._l) * eval(*node._r);
                        } else if (node._l) {
                            return eval(*node._l) * node.constant;
                        } else if (node._r) {
                            return node.constant * eval(*node._l);
                        } else{
                            return 0.0;
                        }
                        },
                    [](const divide_node& node){
                        if (node._l && node._r) {
                            return eval(*node._l) / eval(*node._r);
                        } else if (node._l) {
                            return eval(*node._l) / node.constant;
                        } else if (node._r) {
                            return node.constant / eval(*node._l);
                        }else{
                            return 0.0;
                        }
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
                    [](const constant_node& node){
                        return false;
                    },
                    [idx](const variable_node& node){
                        return node.idx == idx;
                    },
                    [idx](const add_node& node){
                        if(node._l&&node._r){
                            return depend_on(*node._l,idx) || depend_on(*node._r,idx);
                        } else if(node._l){
                            return depend_on(*node._l,idx);
                        }else if(node._r){
                            return depend_on(*node._r,idx);
                        }else{
                            return false;
                        }
                    },
                    [idx](const minus_node& node){
                        if(node._l&&node._r){
                            return depend_on(*node._l,idx) || depend_on(*node._r,idx);
                        } else if(node._l){
                            return depend_on(*node._l,idx);
                        }else if(node._r){
                            return depend_on(*node._r,idx);
                        }else{
                            return false;
                        }
                    },
                    [idx](const multiply_node& node){
                        if(node._l&&node._r){
                            return depend_on(*node._l,idx) || depend_on(*node._r,idx);
                        } else if(node._l){
                            return depend_on(*node._l,idx);
                        }else if(node._r){
                            return depend_on(*node._r,idx);
                        }else{
                            return false;
                        }
                    },
                    [idx](const divide_node& node){
                        if(node._l&&node._r){
                            return depend_on(*node._l,idx) || depend_on(*node._r,idx);
                        } else if(node._l){
                            return depend_on(*node._l,idx);
                        }else if(node._r){
                            return depend_on(*node._r,idx);
                        }else{
                            return false;
                        }
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
                    [](const constant_node& node){
                        return 0.0;
                    },
                    [idx](const variable_node& node){
                        return idx == node.idx? 1.0 : 0.0;
                    },
                    [idx](const add_node& node){
                        double temp = .0;
                        if(node._l&&depend_on(*node._l,idx)){
                            temp += derivative_for(*node._l,idx);
                        }
                        if(node._r&&depend_on(*node._r,idx)){
                            temp += derivative_for(*node._r,idx);
                        }
                        return temp;
                    },
                    [idx](const minus_node& node){
                        double temp = 0.;
                        if(node._l&&depend_on(*node._l,idx)){
                            temp += derivative_for(*node._l,idx);
                        }
                        if(node._r&&depend_on(*node._r,idx)){
                            temp -= derivative_for(*node._r,idx);
                        }
                        return temp;
                    },
                    [idx](const multiply_node& node){
                        auto dfdx_mul_gx =
                                (node._l&&depend_on(*node._l,idx))?
                                derivative_for(*node._l,idx) * (node._r?eval(*node._r):node.constant):
                                0.0;
                        auto dgdx_mul_fx =
                                (node._r&&depend_on(*node._r,idx))?
                                derivative_for(*node._r,idx)* (node._l?eval(*node._l):node.constant):
                                0.0;
                        return dfdx_mul_gx + dgdx_mul_fx;
                    },
                    [idx](const divide_node& node){
                        auto gx = node._r?eval(*node._r):node.constant;
                        auto gx2 = gx*gx;
                        auto dfdx_mul_gx =
                                (node._l&&depend_on(*node._l,idx))?
                                derivative_for(*node._l,idx) * (node._r?eval(*node._r):node.constant):
                                0.0;
                        auto dgdx_mul_fx =
                                (node._r&&depend_on(*node._r,idx))?
                                derivative_for(*node._r,idx)* (node._l?eval(*node._l):node.constant):
                                0.0;
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
        return detail::make_differential_node_shared(add_node{l, r});
    }

    inline detail::nsp operator-(const  detail::nsp &l, const  detail::nsp &r) {
        return detail::make_differential_node_shared(minus_node{l, r});
    }

    inline detail::nsp operator*(const  detail::nsp &l, const  detail::nsp &r) {
        return detail::make_differential_node_shared(multiply_node{l, r});
    }

    inline detail::nsp operator/(const  detail::nsp &l, const  detail::nsp &r) {
        return detail::make_differential_node_shared(divide_node{l, r});
    }

    inline detail::nsp operator+(double l, const  detail::nsp &r) {
        return detail::make_differential_node_shared(add_node{r, l});
    }

    inline detail::nsp operator-(double l, const  detail::nsp&r) {
        return detail::make_differential_node_shared(minus_node{l, r});
    }

    inline detail::nsp operator*(double l, const  detail::nsp &r) {
        return detail::make_differential_node_shared(multiply_node{r, l});
    }

    inline detail::nsp operator/(double l, const  detail::nsp &r) {
        return detail::make_differential_node_shared(divide_node{l, r});
    }

    inline detail::nsp operator+(const  detail::nsp &l, double r){
        return detail::make_differential_node_shared(add_node{l, r});
    }

    inline detail::nsp operator-(const  detail::nsp &l, double r){
        return detail::make_differential_node_shared(minus_node{l, r});
    }

    inline detail::nsp operator*(const  detail::nsp &l, double r){
        return detail::make_differential_node_shared(multiply_node{l, r});
    }

    inline detail::nsp operator/(const  detail::nsp &l, double r){
        return detail::make_differential_node_shared(divide_node{l, r});
    }

    inline detail::nsp sin( const detail::nsp& d) {
        return detail::make_differential_node_shared(sin_node{d});
    }

    inline detail::nsp cos( const detail::nsp& d) {
        return detail::make_differential_node_shared(cos_node{d});
    }

    inline detail::nsp tan( const detail::nsp& d) {
        return detail::make_differential_node_shared(tan_node{d});
    }

    inline detail::nsp pow( const detail::nsp& d, double n){
        return detail::make_differential_node_shared(pow_node{d, n});
    }

    inline detail::nsp pow(double a,  const detail::nsp& d){
        return detail::make_differential_node_shared(exp_node{a, d});
    }

    inline detail::nsp sqrt( const detail::nsp& d){
        return detail::make_differential_node_shared(sqrt_node{d});
    }

    inline detail::nsp exp( const detail::nsp& d){
        return detail::make_differential_node_shared(exp_node{d});
    }

    inline detail::nsp ln( const detail::nsp& d){
        return detail::make_differential_node_shared(ln_node{d});
    }

    inline detail::nsp asin( const detail::nsp& d){
        return detail::make_differential_node_shared(asin_node{d});
    }

    inline detail::nsp acos( const detail::nsp& d){
        return detail::make_differential_node_shared(acos_node{d});
    }

    inline detail::nsp atan( const detail::nsp& d){
        return detail::make_differential_node_shared(atan_node{d});
    }
}
#endif //AUTODIFF_NODE_H
