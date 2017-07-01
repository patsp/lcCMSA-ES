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

function sqrtCnormalized = computeSqrtCNormalized(C, targetCConditionNumber)
  % Cnormalized = computeSqrtCNormalized(C, targetCConditionNumber)
  %    This function computes the normalized covariance matrix.

  dimensionOfNullSpace = size(C, 1);

  C = 0.5 * (C + C');
  [U, D] = eig(C);
  diagD = diag(D);
  diagD(diagD < 0) = 1e-15;
  r = 0;
  if cond(C) > targetCConditionNumber
    eigSorted = sort(diagD);
    r = -sqrt(eigSorted(1)) + ...
        (sqrt(eigSorted(end)) / targetCConditionNumber) + ...
        sqrt((sqrt(eigSorted(1)) - ...
              sqrt(eigSorted(end)) / targetCConditionNumber) ^ 2 - ...
             eigSorted(1) + eigSorted(end) / targetCConditionNumber);
  end
  sqrtC = U * sqrt(diag(diagD)) + r * eye(dimensionOfNullSpace);
  Ctilde = sqrtC * sqrtC';
  sqrtCnormalized = (det(Ctilde) ^ (-1 / (2 * size(C, 1)))) * sqrtC;
end

