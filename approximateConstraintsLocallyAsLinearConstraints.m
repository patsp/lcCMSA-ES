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

function [Aineq, bineq, Aeq, beq] = ...
         approximateConstraintsLocallyAsLinearConstraints(x, sigma, g, h)
  % [Aineq, bineq, Aeq, beq] = ...
  %       approximateConstraintsLocallyAsLinearConstraints(x, sigma, g, h)
  %
  %    This function approximates the inequality constraint function g
  %    and equality constraint function h. It is done by sampling vectors
  %    from a multivariate normal distribution with mean x and identity
  %    covariance matrix scaled with sigma.
  %
  %    x is expected to be a vector of some dimension D
  %
  %    sigma is expected to be a scalar
  %
  %    g is the inequality constraint function. It is expected to be
  %        a function taking a vector of dimension D as input
  %        and returning a vector of K values where K represents the
  %        number of inequality constraints. The K values returned
  %        by a call to g(x) represent the values of the K inequalities
  %        at the position x in parameter space.
  %
  %    h is the equality constraint function. It is expected to be
  %        a function taking a vector of dimension D as input
  %        and returning a vector of M values where M represents the
  %        number of equality constraints. The M values returned
  %        by a call to h(x) represent the values of the M equalities
  %        at the position x in parameter space.
  %
  %    The function returns the linear approximations of g and h in matrix
  %    form: Aineq x <= bineq for a D-dimensional x
  %    and Aeq x = beq for a D-dimensional x.

  D = size(x, 1);
  K = size(g(x), 1);
  M = size(h(x), 1);
  L = 10 * (D + 1);
  Y = zeros(L, D + 1);
  G = zeros(L, K);
  H = zeros(L, M);
  for l = 1:L
    ylhat = x + sigma * randn(D, 1);
    yl = [ylhat; 1];
    Y(l, :) = yl;
    G(l, :) = g(ylhat)';
    H(l, :) = h(ylhat)';
  end
  Ypinv = pinv(Y);
  Wineq = zeros(D + 1, K);
  for k = 1:K
    Wineq(:, k) = Ypinv * G(:, k);
  end
  AineqPrime = Wineq';
  Aineq = AineqPrime(:, 1:D);
  bineq = -AineqPrime(:, D + 1);
  Weq = zeros(D + 1, M);
  for m = 1:M
    Weq(:, m) = Ypinv * H(:, m);
  end
  AeqPrime = Weq';
  Aeq = AeqPrime(:, 1:D);
  beq = -AeqPrime(:, D + 1);
end
