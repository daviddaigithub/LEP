//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai,eeyang. All rights reserved.
//

#include "method.hpp"

/**shuffle the index for cross validation*/
arma::Col<uword> cross_valind(arma::uword N, arma::uword nfold){
    arma::Col<uword> indices(N);
    arma_rng::set_seed(123);
    arma::Col<uword> vec_n = arma::shuffle(arma::linspace < arma::Col <uword> >(1, N, N));

    indices.fill(nfold);
    for(uword n = 1; n <= nfold-1; n++){
        arma::Col<uword> in = vec_n.rows((n-1)*N/nfold,n*N/nfold - 1);
        indices.elem(in - 1 ).fill(n);
    }
    return indices;
}
