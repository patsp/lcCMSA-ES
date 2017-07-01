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

function [A, b, c] = createKleeMintyCubeConstraintSystem(dimension)
  % [A, b] = createKleeMintyCubeConstraintSystem(dimension)
  %    Creates a linear system of equations
  %        Ax = b
  %        x >= 0
  %    that represents the feasible region.
  %
  %    The constraints represent the Klee-Minty cube.
  %    Klee, V., Minty, G.J.:
  %    How good is the simplex algorithm? In: Shisha, O. (ed.) Inequalities III,
  %    pp. 159–175. Academic, New York (1972)
  %
  %    The actual constraints are taken from
  %    https://en.wikipedia.org/wiki/Klee%E2%80%93Minty_cube.
  %
  %     x_1              <= 5
  %    4x_1 +  x_2       <= 25
  %    8x_1 + 4x_2 + x_3 <= 125
  %    .
  %    .
  %    .
  %    2^Dx_1 + 2^(D-1)x_2 + ... + 4x_{D-1} + x_D <= 5^D
  %    x_1 >= 0, ..., x_D >= 0
  %
  %    Introducing a slack variable for every inequality we can
  %    transform this into standard form.
  %
  %    The objective is to maximize
  %    2^(D-1)x_1 + 2^(D-2)x_2 + ... + 2x_{D-1} + x_D.
  %    These coefficients are returned in the 'c' vector.

  A = [];
  b = [];
  c = [];

  if (dimension <= 0)
    error('dimension must be greater than zero');
  end

  numSlackVariables = dimension;
  A = zeros(dimension, dimension + numSlackVariables);
  b = zeros(dimension, 1);
  c = zeros(dimension + numSlackVariables, 1);
  for i = 1:dimension
    A(i, i) = 1;
    for j = 1:(i - 1)
      A(i, j) = 2^(i - j + 1);
    end
    % slack variable for inequality i
    A(i, dimension + i) = 1;
    b(i, 1) = 5^i;
    c(i, 1) = -2^(dimension - i);
  end

end

