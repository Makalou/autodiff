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

    struct differential_node{
    public:
        std::pair<detail::nsp,detail::nsp> dependencies{};

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
    };

    struct variable_node : differential_node{
        
    };

    struct add_node : differential_node{
        add_node(const detail::nsp& l,const detail::nsp& r){}
    };

    struct minus_node : differential_node{

    };

    struct multiply_node : differential_node{

    };

    struct divide_node : differential_node{

    };

    detail::nsp operator+(const detail::nsp& l,const detail::nsp& r){
        return std::make_shared<add_node>(l,r);
    }

    detail::nsp operator+(double l,const detail::nsp& r){
        return std::make_shared<constant_node>(l) + r;
    }

    detail::nsp operator+(const detail::nsp& l,double r){
        return l + std::make_shared<constant_node>(r);
    }
}
#endif //AUTODIFF_NODE_H
