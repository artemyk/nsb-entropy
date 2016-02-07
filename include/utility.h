/*
 *  File - utility.h
 *  Project - Entropy_parsing
 * 
 *  Numbers to strings conversions utilities; checks for appropriate compiler
 *  and includes appropriate headers.
 *
 *  Created by nemenman on Fri Nov 02 2001.
 *  Copyright (c) 2001 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef UTILITY_H
#define UTILITY_H

#include <stdio.h>
#include <stdexcept>
#include <string>
#include <iostream>

// see which compiler we are using
#include "compiler.h"
#ifdef GCC_2
#include "utility_2.h"
#else
#include "utility_3.h"
#endif



//concatenating a string and a long
inline std::string operator+(const std::string& s, const long i){
    return s+num2str(i);
}


//concatenating a string and an unsigned long
inline std::string operator+(const std::string& s, const unsigned long i){
    return s+num2str(i);
}


//concatenating a long and a string
inline std::string operator+(const long i, const std::string& s){
    return num2str(i)+s;
}

//concatenating a string and a int 
inline std::string operator+(const std::string& s, const int i){
    return s+num2str(i);
}

//concatenating an int and a string
inline std::string operator+(const int i, const std::string& s){
    return num2str(i)+s;
}

//concatenating a string and a double
inline std::string operator+(const std::string& s, const double d){
    return s+num2str(d);
}

//concatenating an double and a string
inline std::string operator+(const double d, const std::string& s){
    return num2str(d)+s;
}


#endif
