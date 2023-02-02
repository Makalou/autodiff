//
// Created by 王泽远 on 2023/2/2.
//

#ifndef AUTODIFF_INTRUSIVE_PTR_H
#define AUTODIFF_INTRUSIVE_PTR_H

template<typename T>
struct intrusive_ptr{
    T* ptr;

    intrusive_ptr(T* _ptr) : ptr(_ptr){

    }

    intrusive_ptr(intrusive_ptr const& other){
        ptr = other.ptr;
        ptr->inc_count();
    }

    intrusive_ptr& operator=(intrusive_ptr const& other){
        ptr = other.ptr;
        ptr->inc_count();
        return *this;
    }

    ~intrusive_ptr(){
        if(ptr&&ptr->dec_count() <= 0){
            delete ptr;
        }
    }
};

#endif //AUTODIFF_INTRUSIVE_PTR_H
