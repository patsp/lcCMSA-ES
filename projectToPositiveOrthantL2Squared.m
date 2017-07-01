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

function xProjected = projectToPositiveOrthantL2Squared(x, A, b, ...
                                                        xInit, upperBoundScalar)
  % xProjected = projectToPositiveOrthantL2Squared(x, A, b, ...
  %                                                xInit, upperBoundScalar)
  %    Projects the given x into the positive orthant, i.e.,
  %    solves the optimization problem
  %    min_{x_opt}{||x_opt - x||^2}
  %    such that x >= 0 and Ax = b.
  %    Informally it computes an x_opt in the positive orthant
  %    with minimal squared euclidean distance to the given
  %    x and fulfilling the constraints.
  %
  %    This optimization can be written as a standard
  %    quadratic programming problem.
  %
  %    min_{x_opt}{||x_opt - x||^2}
  %    = min_{x_opt}{sum_i{(x_opt_i - x_i)^2}}
  %    = min_{x_opt}{sum_i{x_opt_i^2 - 2*x_i*x_opt_i + x_i^2}}
  %    = min_{x_opt}{x_opt' * I * x_opt - 2*x'*x_opt + x'*x}
  %    = min_{x_opt}{0.5 * x_opt' * I * x_opt - x_opt' * x + 0.5 * x'*x}
  %    = min_{x_opt}{0.5 * x_opt' * I * x_opt - x_opt' * x}
  %      (constant term does not matter in optimization)

  dimension = length(xInit);
  x0 = xInit;
  H = eye(dimension);
  q = -x;
  lowerBound = zeros(dimension, 1);
  upperBound = upperBoundScalar * ones(dimension, 1);
  [xProjected, ~, ~, ~] = qp(x0, H, q, A, b, lowerBound, upperBound);
end

