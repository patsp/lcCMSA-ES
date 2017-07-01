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

function [A, b] = transformToStandardFormEq(W, c, lbnds, ubnds, K)
  % [A, b] = transformToStandardFormEq(W, c, lbnds, ubnds, K)
  %    This function transforms a system
  %        Wy = c, lbnds <= y <= ubnds
  %    into
  %        Ax = b, x >= 0.

  [M, D] = size(W);
  A = [W         -W         zeros(M, K) zeros(M, D) zeros(M, D);
       eye(D, D) -eye(D, D) zeros(D, K) -eye(D, D)  zeros(D, D);
       eye(D, D) -eye(D, D) zeros(D, K) zeros(D, D) eye(D, D)];
  b = [c;
       lbnds;
       ubnds];
end
