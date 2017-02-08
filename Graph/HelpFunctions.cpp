/* This is a header file for the graph library */

#include "./stdc++.h" // include standard c++ library, this file is included in the project folder

/* help functions: (1) generate random number between (0, kS) inside of file*/ 
int getRandomNumber(int kS) {
    int o = rand() % kS;
    //cout<<"randome number is:" << o<<endl;
    return o; // kS;
    }
    
int randomDistance() {
    const double kStart = 1.0;
    const double kEnd = 10.0;
    int k = (int)(rand() % 100 + 1)/10;  // set the random number edges' value is from 1 to 10
    //cout<<"randome distance is:" << k<<endl;
    return k;
    }