//
// Created by 王泽远 on 2023/1/8.
//

#ifndef AUTODIFF_DUAL_H
#define AUTODIFF_DUAL_H

#include "numeric"
#include <cmath>
#include "../auto_diff.h"

namespace ad {

    struct dual_number {
    public:
        double _dual_part;
        double _real_part;
    public:

        dual_number(double real_part,double dual_part): _real_part(real_part),_dual_part(dual_part){};

        //todo Caution : FP comparison
        bool operator == (const dual_number &other) const {
            return _real_part == other._real_part;
        }

        bool operator != (const dual_number &other) const{
            return _real_part != other._real_part;
        }

        bool operator < (const dual_number &other) const{
            return _real_part < other._real_part;
        }

        bool operator <= (const dual_number &other) const{
            return _real_part <= other._real_part;
        }

        bool operator > (const dual_number &other) const{
            return _real_part > other._real_part;
        }

        bool operator >= (const dual_number &other) const{
            return _real_part >= other._real_part;
        }

        // Basic algebra operation
        dual_number operator-() const{
            return dual_number{-_real_part,_dual_part};
        }

        void operator+=(const dual_number &other) {
            _real_part += other._real_part;
            _dual_part += other._dual_part;
        }

        void operator-=(const dual_number &other) {
            _real_part -= other._real_part;
            _dual_part -= other._dual_part;
        }

        void operator*=(const dual_number &other) {
            _real_part *= other._real_part;
            _dual_part = other._real_part * _dual_part + _real_part * other._dual_part;
        }

        void operator/=(const dual_number &other) {
            _real_part /= other._real_part;
            _dual_part = other._real_part / _dual_part + _real_part / other._dual_part;
        }

    };

    inline bool operator == (double x,const dual_number &dual){
        return x == dual._real_part;
    }

    inline bool operator == (const dual_number &dual,double x){
        return x == dual._real_part;
    }

    inline bool operator != (double x,const dual_number &dual){
        return x != dual._real_part;
    }

    inline bool operator != (const dual_number &dual,double x){
        return x != dual._real_part;
    }

    inline bool operator < (double x,const dual_number &dual){
        return x < dual._real_part;
    }

    inline bool operator < (const dual_number &dual,double x){
        return x < dual._real_part;
    }

    inline bool operator <= (double x,const dual_number &dual){
        return x <= dual._real_part;
    }

    inline bool operator <= (const dual_number &dual,double x){
        return x <= dual._real_part;
    }

    inline bool operator > (double x,const dual_number &dual){
        return x > dual._real_part;
    }

    inline bool operator > (const dual_number &dual,double x){
        return x > dual._real_part;
    }

    inline bool operator >= (double x,const dual_number &dual){
        return x >= dual._real_part;
    }

    inline bool operator >= (const dual_number &dual,double x){
        return x >= dual._real_part;
    }

    inline dual_number operator+(const  dual_number &l, const  dual_number &r) {
        return dual_number{l._real_part + r._real_part, l._dual_part + r._dual_part};
    }

    inline dual_number operator-(const  dual_number &l, const  dual_number &r) {
        return dual_number{l._real_part - r._real_part, l._dual_part - r._dual_part};
    }

    inline dual_number operator*(const  dual_number &l, const  dual_number &r) {
        return dual_number{l._real_part * r._real_part, r._real_part * l._dual_part + l._real_part * r._dual_part};
    }

    inline dual_number operator/(const  dual_number &l, const  dual_number &r) {
        return dual_number{l._real_part / r._real_part, r._real_part / l._dual_part + l._real_part / r._dual_part};
    }

    inline dual_number operator+(double l, const  dual_number &r) {
        return dual_number{l + r._real_part, r._dual_part};
    }

    inline dual_number operator-(double l, const  dual_number&r) {
        return dual_number{l - r._real_part, - r._dual_part};
    }

    inline dual_number operator*(double l, const  dual_number &r) {
        return dual_number{l * r._real_part, l * r._dual_part};
    }

    inline dual_number operator/(double l, const  dual_number &r) {
        return dual_number{l / r._real_part,  l / r._dual_part};
    }

    inline dual_number operator+(const  dual_number &l, double r){
        return dual_number{l._real_part + r, l._dual_part};
    }

    inline dual_number operator-(const  dual_number &l, double r){
        return dual_number{l._real_part - r, l._dual_part};
    }

    inline dual_number operator*(const  dual_number &l, double r){
        return dual_number{l._real_part * r, l._dual_part * r};
    }

    inline dual_number operator/(const  dual_number &l, double r){
        return dual_number{l._real_part / r,  l._dual_part / r};
    }

    inline dual_number sin( dual_number d) {
        return dual_number{std::sin(d._real_part),std::cos(d._real_part)*d._dual_part};
    }

    inline dual_number cos( dual_number d) {
        return dual_number{std::cos(d._real_part),-std::sin(d._real_part)*d._dual_part};
    }

    inline dual_number tan( dual_number d) {
        auto sec = 1/std::cos(d._real_part);
        auto sec2 = sec*sec;
        return dual_number{std::tan(d._real_part),sec2*d._dual_part};
    }

    inline dual_number pow( dual_number d, double n){
        return dual_number{std::pow(d._real_part,n),n*std::pow(d._real_part,n-1)*d._dual_part};
    }

    inline dual_number pow(double a,  dual_number d){
        return dual_number{std::pow(a,d._real_part),std::log(a)*std::pow(a,d._real_part)*d._dual_part};
    }

    inline dual_number pow(dual_number d1,dual_number d2){
        return dual_number{std::pow(d1._real_part,d2._real_part),(pow(d1._real_part,d2)+ pow(d1,d2._real_part))._dual_part};
    }

    inline dual_number sqrt( dual_number d){
        return dual_number{std::sqrt(d._real_part),0.5/std::sqrt(d._real_part)*d._dual_part};
    }

    inline dual_number exp( dual_number d){
        return dual_number{std::exp(d._real_part),std::exp(d._real_part)*d._dual_part};
    }

    inline dual_number ln( dual_number d){
        return dual_number{std::log(d._real_part),1/d._real_part*d._dual_part};
    }

    inline dual_number asin( dual_number d){
        auto dx = 1/(std::sqrt(1-d._real_part*d._real_part));
        return dual_number{std::asin(d._real_part),dx*d._dual_part};
    }

    inline dual_number acos( dual_number d){
        auto dx = -1/(std::sqrt(1-d._real_part*d._real_part));
        return dual_number{std::acos(d._real_part),dx*d._dual_part};
    }

    inline dual_number atan( dual_number d){
        auto dx = 1/(1+d._real_part*d._real_part);
        return dual_number{std::atan(d._real_part),dx*d._dual_part};
    }

}//end namespace ad

#endif //AUTODIFF_DUAL_H
