% Copyright 2017 The lcCMSA-ES Authors.  All Rights Reserved.
%
% This file is part of lcCMSA-ES.
%
% lcCMSA-ES is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% lcCMSA-ES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with lcCMSA-ES.  If not, see <http://www.gnu.org/licenses/>.

function xProjected = projectToPositiveOrthantL1(x, A, b)
  % xProjected = projectToPositiveOrthantL1(x, A, b)
  %    Projects the given x into the positive orthant, i.e.,
  %    solves the optimization problem
  %    min_{x_opt}{||x_opt - x||_1}
  %    such that x >= 0 and Ax = b.
  %    Informally it computes an x_opt in the positive orthant
  %    with minimal L1 distance to the given
  %    x and fulfilling the constraints.
  %
  %    This optimization can be written as a standard
  %    linear programming problem.
  %
  %    min_{x_opt}{||x_opt - x||_1}
  %    = min_{x_opt}{sum_i{|(x_opt_i - x_i)|}}
  %
  %    min 1' * z
  %    such that z >= 0
  %              z >= x_opt - x
  %              z >= -x_opt + x
  %              A * x_opt = b
  %              x_opt >= 0
  %
  %    Rewriting the constraints for the variables to optimize to be
  %    on the lhs yields
  %
  %    min 1' * z
  %    such that z >= 0
  %              z - x_opt >= -x
  %              z + x_opt >= x
  %              A * x_opt = b
  %              x_opt >= 0.
  %
  %    Introducing slack variable vectors s_1, s_2, s_3 we arrive at the
  %    final lp
  %
  %    min 1' * z
  %    such that z >= 0
  %              A *  x_opt                       =  b
  %                  -x_opt + z - s_1             = -x
  %                   x_opt + z       - s_2       =  x
  %                   x_opt                 - s_3 = 0.

  dimension = length(x);
  I = eye(dimension);
  Z = zeros(dimension);
  ZA = zeros(size(A));
  M = [ A, ZA,  ZA,  ZA,  ZA;
        -I,  I,  -I,   Z,   Z;
        I,  I,   Z,  -I,   Z;
        I,  Z,   Z,   Z,  -I];
  y = [ b;
        -x;
        x;
        zeros(dimension, 1)];
  c = [ zeros(dimension, 1);
        ones(dimension, 1);
        zeros(dimension, 1);
        zeros(dimension, 1);
        zeros(dimension, 1)];
  lowerBound = zeros(length(c), 1);
  upperBound = [];
  %[opt, ~, ~, ~] = glpk(c, M, y, lowerBound, upperBound);
  opt = solveNonNegativeLpWithLpSolve(c, M, y);
  xProjected = opt(1:dimension, 1);
end
