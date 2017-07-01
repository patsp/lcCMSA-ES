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

function y = backtransformStandardFormVector(x, D, K)
  % y = backtransformStandardFormVector(x, D, K)
  %    This function transforms a vector in standard form back
  %    to a vector in the original system.

  I = eye(D, D);
  Z = zeros(D, D);
  Z2 = zeros(D, K);
  F = [I -I Z2 Z Z];
  y = F * x;
  y = y(1:D);
end
