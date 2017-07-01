// Copyright 2017 The lcCMSA-ES Authors.  All Rights Reserved.
//
// This file is part of lcCMSA-ES.
//
// lcCMSA-ES is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lcCMSA-ES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lcCMSA-ES.  If not, see <http://www.gnu.org/licenses/>.

#include "mex.h"

#include <vector>

#include <lpsolve/lp_lib.h>

static void solveNonNegativeLp(const mxArray *c,
                               const mxArray *A,
                               const mxArray *b,
                               mxArray *x) {
    const size_t m = mxGetM(A);
    const size_t n = mxGetN(A);

    const double *cData = mxGetPr(c);
    const double *AData = mxGetPr(A);
    const double *bData = mxGetPr(b);
    double *xData = mxGetPr(x);

    // + 1 because first value of every column is the objective
    std::vector<double> column(m + 1, 0);

    // by default, we have a nonnegative LP
    // columns added iteratively below
    lprec *lp = make_lp(m, 0);
    if (!lp) {
        mexErrMsgTxt("Unable to create new LP model");
    }

    set_verbose(lp, NEUTRAL);

    for (size_t col = 0; col < n; ++col) {
        column.at(0) = cData[col]; // objective value
        for (size_t row = 0; row < m; ++row) {
            size_t index = col * m + row; // column major
            column.at(row + 1) = AData[index];
        }
        if (!add_column(lp, &column[0])) {
            mexPrintf("Unable to add column\n");
            goto done;
        }
    }

    for (size_t row = 0; row < m; ++row) {
        size_t index = row;
        column.at(row + 1) = bData[index];
    }
    set_rh_vec(lp, &column[0]);

    for (size_t row = 0; row < m; ++row) {
        if (!set_constr_type(lp, row + 1, EQ)) {
            mexPrintf("Unable to set constraint type\n");
            goto done;
        }
    }

    if (0 != solve(lp)) {
        mexPrintf("Unable to find optimal solution to LP\n");
        goto done;
    }

    for (size_t i = 0; i < n; ++i) {
        xData[i] = get_var_primalresult(lp, m + 1 + i);
    }

 done:
    delete_lp(lp);
}

// x = solveNonNegativeLpWithCgal(c, A, b)
// Solves min c'*x s.t. Ax = b, x >= 0.
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (3 != nrhs) {
        mexErrMsgTxt("Three inputs required.");
    }
    const mxArray *c = prhs[0];
    const mxArray *A = prhs[1];
    const mxArray *b = prhs[2];

    if (!mxIsDouble(c) ||
          mxIsComplex(c)) {
        mexErrMsgTxt("Input vector c must be type double.");
    }

    if (!mxIsDouble(A) ||
          mxIsComplex(A)) {
        mexErrMsgTxt("Input matrix A must be type double.");
    }

    if (!mxIsDouble(b) ||
           mxIsComplex(b)) {
        mexErrMsgTxt("Input vector b must be type double.");
    }

    if (!(1 == mxGetN(b) && mxGetM(A) == mxGetM(b))) {
        mexErrMsgTxt("Input vector b must be a column vector with the length "
                     "being the same as the number of rows of A.");
    }

    if (!(1 == mxGetN(c) && mxGetN(A) == mxGetM(c))) {
        mexErrMsgTxt("Input vector c must be a column vector with the length "
                     "being the same as the number of columns of A.");
    }

    if (1 != nlhs) {
        mexErrMsgTxt("One output required.");
    }

    plhs[0] = mxCreateDoubleMatrix(mxGetM(c), 1, mxREAL);
    mxArray *x = plhs[0];
    solveNonNegativeLp(c, A, b, x);
}

