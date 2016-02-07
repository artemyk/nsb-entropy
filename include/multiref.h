/*
 *  File - multiref.h
 *
 *  Declarations of the class used to do NSB numerics
 *
 *  Copyright (c) 2001, 2002 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef MULTIREF_H
#define MULTIREF_H

class Ptr_base { 
public: 
  class Base { 
  private: 
    int m_refs; 
    friend class Ptr_base; 
  protected: 
    Base() : m_refs(0) {} 
    virtual ~Base(); 
  }; 

private: 
  Base * m_ptr; 
  
  void addref() 
  { if (m_ptr != 0) m_ptr -> m_refs++; } 
  void subref() { 
    if (m_ptr != 0 && --(m_ptr -> m_refs) == 0) {
      delete m_ptr; 
      m_ptr = 0; 
    }
  } 

protected: 
  Ptr_base(Base * p) : m_ptr(p) { addref(); } 
  Ptr_base(const Ptr_base& p) : m_ptr(p.m_ptr) { addref(); } 
  ~Ptr_base() { subref(); } 

  Ptr_base& operator = (Base * p) { 
    if (p != m_ptr) { 
      subref(); 
      m_ptr = p; 
      addref(); 
    } 
    return *this;
  }; 
    
  Base * value() { return m_ptr; } 
  const Base * value() const { return m_ptr; } 
}; 

template <class Derived> 
class Ptr : public Ptr_base { 
public: 
  Ptr() : Ptr_base(0) {} 
  Ptr(Derived* d) : Ptr_base(d) {} 
  Ptr(const Ptr& p) : Ptr_base(p) {} 
  // no need for explicit destructor here. 

  operator Derived * () { return value(); } 
  operator const Derived * () const{ return value(); }	

  Derived& operator *(){ return *value(); } 
  const Derived& operator *()const{ return *value(); } 
  Derived* operator->(){ return value(); } 
  //  const Derived* operator->(){ return value(); } 


}; 

#endif
