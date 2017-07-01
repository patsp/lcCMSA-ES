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

function xProjected = projectToPositiveOrthantIter(x, A, b, p)
  % xProjected = projectToPositiveOrthantIteratively(x, A, b, pOpt)
  %    This function projects the given x into the positive orthant.
  %    The projection is done by moving towards p
  %    in the null space of A.
  %    This function assumes that A * p = b and p >= 0.

  dimension = length(x);
  xProjected = x;
  d = p - x;
  alpha = 0;
  for k = 1:dimension
    tmp = xProjected(k);
    if tmp < -1e-20
      if abs(d(k)) > 1e-20
        alpha = max(alpha, (-1e-20 - tmp) / d(k));
      end
    end
  end
  xProjected = xProjected + alpha * d;

end
