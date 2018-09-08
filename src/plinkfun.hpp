//
//  plinkfun.hpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef plinkfun_hpp
#define plinkfun_hpp

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
using namespace std;


void readPlink(string stringname,int N, int P, char* X);
int getLineNum(string filename);

#endif /* plinkfun_hpp */
