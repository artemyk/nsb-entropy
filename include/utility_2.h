/*
 *  File - utility.h
 *  Project - Entropy_parsing
 * 
 *  Numbers to strings conversions utilities; for GCC 2.95
 *
 *  Created by nemenman on Fri Nov 02 2001.
 *  Copyright (c) 2001 Ilya Nemenman. All rights reserved.
 *
 */

#include <strstream>

const int NUM2STR_LEN = 50;

//conversion of a long to a string
inline string num2str(const long i){
    char  str[NUM2STR_LEN];
    int pr=snprintf(str,NUM2STR_LEN,"%D",i);
    if(pr>=NUM2STR_LEN) throw out_of_range("Out of range in 'string num2str(const long)'");
    return string(str);
}

//conversion of an unsigned long to a string
inline string num2str(const unsigned long i){
    char  str[NUM2STR_LEN];
    int pr=snprintf(str,NUM2STR_LEN,"%D",i);
    if(pr>=NUM2STR_LEN) throw out_of_range("Out of range in 'string num2str(const long)'");
    return string(str);
}

//conversion of an int to a string
inline string num2str(const int i){
    char  str[NUM2STR_LEN];
    int pr=snprintf(str,NUM2STR_LEN,"%i",i);
    if(pr>=NUM2STR_LEN) throw out_of_range("Out of range in 'string num2str(const int)'");
    return string(str);
}

//conversion of a double to a string
inline string num2str(const double d){
    char  str[NUM2STR_LEN];
    int pr=snprintf(str,NUM2STR_LEN,"%f",d);
    if(pr>=NUM2STR_LEN) throw out_of_range("Out of range in 'string num2str(const double)'");
    return string(str);
}

