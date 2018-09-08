//
//  readPlink.cpp
//  ReadPlinkGit
//
//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef readPlink_hpp
#define readPlink_hpp

#include<iostream>
#include <sstream>
#include <stdio.h>
#include <boost/range.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/algorithm/string.hpp>

#include "plinkfun.hpp"
#include "aux.hpp"
#include <map>

using namespace std;
using namespace arma;
using namespace boost;

//GenoInfo ReadDataFromFile(string stringname);
class Chroms;
class SNP;

class SNP{
public:
    SNP(){

    }
    SNP(string name, int idx, int type){
        this -> name = name;
        this -> idx = idx;
        this -> type = type;
    }
    SNP(const SNP& snp){
        this -> name = snp.name;
        this -> idx = snp.idx;
        this -> type = snp.type;
    }
    string name;
    int idx;
    int type;
    bool operator<(const SNP&) const;
    bool operator>(const SNP&) const;
    bool operator != (const SNP&) const;
    bool operator == (const SNP&) const;

    //   friend ostream& operator<<(ostream& os, const SNP& obj);
};


const int const_chrom_no = 22;
Chroms read_snpnames(string filename, int P);
class Chroms{
public:
    //      chromsomes[const_chrom_no];
    Col<int> chromsome;
    Col<int> A1Array;
    Col<int> A2Array;
    vector<SNP> snps;
  //  vector<string> snpidvector;
    int chrom_no;
    int P;
    Chroms(){

    }
    arma::Col<int> index;
    //to allocate the chromsomes for the chroms, chr_idx = -1 for one chromsome, chr_idx != -1 for const_chrom_no 22
    Chroms(int P){
        chromsome.resize(P);
        A1Array.resize(P);
        A2Array.resize(P);
    //    snps = new SNP[P];
        this -> P = P;
    }

    Chroms(string bimfile, int P){
        *this = read_snpnames(bimfile, P);
    }

    Chroms(const Chroms& chrom){
      this -> chromsome = chrom.chromsome;
      this -> A1Array = chrom.A1Array;
      this -> A2Array = chrom.A2Array;
      this -> snps = chrom.snps;
      this -> index = chrom.index;
      this -> chrom_no = chrom.chrom_no;
      this -> P = chrom.P;
    }

    Chroms(vector<string>& snpnames){
      this -> P = snpnames.size();
      for(int i = 0; i < this -> P; i++){
        SNP snp(snpnames[i], i, 1);
        snps.push_back(snp);
      }
    }

    void clear(){
      chromsome.reset();
      A1Array.reset();
      A2Array.reset();
      index.reset();
    }
};

#define pvalue_v 0
#define zvalue_v 1

#define beta_ss 0
#define pvalue_ss 1

#define from_ss 0
#define from_x 1

/*get the positions of the identifiers in the column names*/
Col<int> getPositions(vector <string> fields, vector<string>& identifiers);


vec read_phenotypes(string filename, int N, int phenopos = 5);

float normalCFD(float value);

int snps_overlap(vector<SNP>& chrom_x_i, vector<SNP>& chrom_s_i, Col<uword>& xindex, Col<uword>& sindex);

template <typename T>
void keep_indices(arma::Mat<T>& mat_obj, int type, uvec index);
Chroms getChromsY(string stringname, vec& y, int& N, int& P);

#endif
