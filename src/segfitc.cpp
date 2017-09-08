/**
 * @file   main.cpp
 * @author yuhonglin <yuhonglin1986@gmail.com>
 * @date   Sun Dec 22 10:39:41 2013
 *
 * @brief  SegFit: an algorithm combining segmentation and fitting
 *
 *
 */

#include <vector>
#include <string>
#include <iostream>

#include <cstdlib>

#include "option.hpp"
#include "logger.hpp"
#include "SegAlgFactory.hpp"
#include "SegAlg.hpp"

#include <R.h>
#include <Rinternals.h>


using namespace std;


extern "C" SEXP yasegfitc(SEXP s, SEXP alg, SEXP fit, SEXP smp, SEXP nseg) {

    LOGGER->enable_exception = false;
    
    // compute
    string algtype(CHAR(STRING_ELT(alg,0)));
    string fittype(CHAR(STRING_ELT(fit,0)));

    SegAlgFactory factory;
    unique_ptr<SegAlg> segalg = factory.make(algtype);

    if (algtype=="topdown") {
	segalg->set_parameter("numseg", *REAL(nseg));
    }

    if (algtype=="dp") {
	segalg->set_parameter("eta", *REAL(smp));
    }

    segalg->set_fitalg(fittype);

    shared_ptr<vector<double>> seqdbptr(new vector<double>(length(s)));
    memcpy(seqdbptr->data(), REAL(s), sizeof(double)*seqdbptr->size());

    segalg->set_string(seqdbptr);

    segalg->run();

    auto result = segalg->get_result();

    // allocate return value
    SEXP ret;
    PROTECT(ret = Rf_allocMatrix(REALSXP, 7, result->size()));

    int i = 0;
    for ( vector<Segment>::const_iterator iter = result->begin();
	  iter != result->end(); iter++ )
    {
	REAL(ret)[i]   = static_cast<double>(iter->headIndex) + 1;
	REAL(ret)[++i] = static_cast<double>(iter->tailIndex) + 1;
	REAL(ret)[++i] = iter->a;
	REAL(ret)[++i] = iter->b;
	REAL(ret)[++i] = iter->c;
	REAL(ret)[++i] = static_cast<double>(iter->order);
	REAL(ret)[++i] = static_cast<double>(iter->loss);
	
	++i;
    }

    UNPROTECT(1);
    
    return ret;
}
