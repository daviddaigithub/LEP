
#include "readPlink.hpp"



float normalCFD(float value)
{
    return 0.5 * erfc(-value * M_SQRT1_2);
}


/*get the positions of the identifiers in the column names*/
Col<int> getPositions(vector <string> fields, vector<string>& identifiers){

    uword size = identifiers.size();
    Col<int> pos(size);
    for(uword i = 0; i < size; i++){
        pos[i] = -1;
        for(uword j = 0; j < fields.size(); j++){
            if(fields[j].compare(identifiers[i]) == 0){
                pos[i] = (int)j ;
                break;
            }
        }

    }
    return pos;

}

Chroms getChromsY(string stringname, vec& y, int& N, int& P){
  string famfile = stringname;
  famfile += ".fam";
  string bimfile = stringname;
  bimfile += ".bim";
  N =  getLineNum(famfile);
  P =  getLineNum(bimfile);
  Chroms chroms(bimfile, P);
  int phenotype_pos = 5;
  y = read_phenotypes(famfile, N, phenotype_pos);
  return chroms;
}




int snps_overlap(vector<SNP>& chrom_x_i, vector<SNP>& chrom_s_i, Col<uword>& xindex, Col<uword>& sindex){

    vector<SNP> common_snp_in_x;
    vector<SNP> common_snp_in_ss;
    sort(chrom_s_i.begin(), chrom_s_i.end());
    sort(chrom_x_i.begin(), chrom_x_i.end());

    set_intersection(chrom_x_i.begin(),chrom_x_i.end(),chrom_s_i.begin(),chrom_s_i.end(), back_inserter(common_snp_in_x));

    set_intersection(chrom_s_i.begin(),chrom_s_i.end(), chrom_x_i.begin(),chrom_x_i.end(), back_inserter(common_snp_in_ss));

    Mat<uword> indexPair(common_snp_in_x.size(), 2);
    for(int i = 0; i < common_snp_in_x.size(); i++){
        indexPair(i,0) = (uword)common_snp_in_x[i].idx;
        indexPair(i,1) = (uword)common_snp_in_ss[i].idx;
    }

    uvec sorted_index = sort_index(indexPair.col(0));
    indexPair = indexPair.rows(sorted_index);
    xindex = indexPair.col(0);
    sindex = indexPair.col(1);

    return (int)common_snp_in_x.size();

}

vec read_phenotypes(string filename, int N, int phenopos){
    std::ifstream ifs(filename.c_str());

    std::string line;
    vec y(N);
    string snpname;
    float pheno = 0;
    string tmp_str;
    int idx = 0;
    while(std::getline(ifs, line)) // read one line from ifs
    {
        std::istringstream iss(line); // access line as a stream
        vector<string> fields;
        boost::split( fields, line, boost::is_any_of(" \t"));
        pheno = (double)atof(fields[phenopos].c_str());
        y[idx] = pheno;
        idx++;
    }
    ifs.close();
    return y;
}

Chroms read_snpnames(string filename, int P){
  std::ifstream ifs(filename.c_str());

  std::string line;
  Chroms chroms(P);//snps(P);
  int chromsome;
  string snpname;
  vector <string> fields;
  int i = 0;
  while(std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    boost::split( fields, line, boost::is_any_of(" \t *"));
    chromsome = (int)atoi(fields[0].c_str());
    snpname = fields[1];
    int a1 = *fields[4].c_str();
    int a2 = *fields[5].c_str();
    chroms.chromsome[i] = chromsome;
    SNP snp(snpname, (int)i, from_x);
    chroms.snps.push_back(snp);
    chroms.A1Array[i] = a1;
    chroms.A2Array[i] = a2;
    i++;

  }
  ifs.close();
  return chroms;
}






bool SNP::operator<(const SNP& obj) const{
    return this -> name < obj.name;
}
bool SNP::operator>(const SNP& obj) const{
    return this -> name > obj.name;
}

bool SNP::operator != (const SNP& obj) const{
    return this -> name.compare(obj.name) > 0;
}


bool SNP::operator == (const SNP& obj) const{
    return this -> name.compare(obj.name) == 0;
}
