//
//  plinkfun.cpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#include "plinkfun.hpp"
#include <stdio.h>
#include <math.h>
#include <bitset>

int getLineNum(string filename){
    FILE *pf = fopen(filename.c_str(), "r"); // 打开文件
    char buf[10000];
    int lineCnt = 0;
    if (!pf) // 判断是否打开成功
        return -1;
    while (fgets(buf, 10000, pf)) // fgets循环读取，直到文件最后，才会返回NULL
        lineCnt++; // 累计行数
    fclose(pf);
    return lineCnt;
}

void getFourGentype(char* geno, std::bitset<8> bits){
  int idx = 0;
  for (int j=0; j < 8; j = j + 2) {
    if(bits[j] && bits[j+1]){
      geno[idx] = 0;
    }else if(!bits[j] && !bits[j+1]){
      geno[idx] = 2;
    }else if(!bits[j] && bits[j+1]){
      geno[idx] = 1;
    }else if(bits[j] && !bits[j+1]){
      geno[idx] = 0;
    }
    idx++;
  }
}


void readPlink(string plinkfile,int N, int P, char* X){
  FILE *fp;
  unsigned char buff[3];
  string bedfile = plinkfile +".bed";
  fp = fopen(bedfile.c_str(), "rb");
  if (!fp) return;
  fread(buff, sizeof(char), 3, fp);

  std::bitset<8> magic1(buff[0]);
  std::bitset<8> magic2(buff[1]);
  std::bitset<8> mode0(buff[2]);

  if(magic1.to_ulong() != 108 || magic2.to_ulong() != 27){
    //   cout <<"Error Identifier of plink binary file" << endl;
  }
  unsigned long mode =  mode0.to_ulong();
  if(mode == 0){
    printf ("individual-Major Order:improper type of plink file");
    exit (EXIT_FAILURE);
  }
  long n = 0;
  long long charNum = ceil(N*1.0/4)*10000;
  long long leftGenoNum = ceil(N*1.0/4)*P;
  long nblock = ceil(N*1.0/4);
  long nSNP = 0;
  while (!feof(fp)) {
    if(leftGenoNum <= 0)
      break;
    if(leftGenoNum <= charNum){
      charNum  = leftGenoNum;
    }
    char* genotype = new char[charNum];
    fread(genotype, sizeof(char), charNum, fp);
    char* geno = new char[4];
    long nSNPc = long(charNum / nblock); //number of SNPs of this iteration
    long long idx = 0;
    for (long i=0; i < nSNPc; i++) {
      for(long j=0; j < nblock - 1; j++){
        long long indx = (long long)(nSNP) * (long long)(N) + (long long)(j*4);
        std::bitset<8> bits(genotype[idx]);
        getFourGentype(geno,bits);
        memcpy(X + indx, geno, 4*sizeof(char));
        idx++;
        leftGenoNum -= 1;
      }
      long left = N - (nblock - 1)*4;
      std::bitset<8> bits(genotype[idx]);
      getFourGentype(geno,bits);

      long long indx2 = (long long)(nSNP) * (long long)(N) + (long long)(nblock - 1)*4;
      long long indx3 = left*sizeof(char);
      memcpy(X + indx2, geno, indx3);
      idx++;
      leftGenoNum -= 1;
      nSNP ++;
    }

    delete[] geno;
    delete[] genotype;
    n++;
  }
}
