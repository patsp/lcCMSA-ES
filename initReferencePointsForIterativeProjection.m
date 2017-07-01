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

function pOpts = initReferencePointsForIterativeProjection(x, nPOpts, A, b)
  % pOpts = initReferencePointsForIterativeProjection(x, nOpts, A, b)
  %    This function initializes the set of nOpts reference points for
  %    the Iterative Projection. The are sampled from a uniform distribution
  %    bounded by the norm of x and projected onto the feasible region
  %    represented by Ax = b, x >= 0.

  dimension = size(A, 2);
  pOpts = zeros(dimension, nPOpts);
  lo = -100 * norm(x);
  hi = 100 * norm(x);
  for k = 1:nPOpts
    p = lo + (hi - lo) * rand(dimension, 1);
    pOpt = projectToPositiveOrthantL1(p, A, b);
    pOpts(:, k) = pOpt;
  end
end
