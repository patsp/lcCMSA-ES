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

[A, b, c] = createKleeMintyCubeConstraintSystem(5);

[xopt, fmin, ~, ~] = glpk(c, A, b, zeros(length(c), 1));
disp('glpk xopt'), disp(xopt);
disp('glpk fmin'), disp(fmin);

config.dimension = size(A, 2);
config.lambda = 4 * config.dimension;
config.mu = config.lambda / 4;
config.gStop = 10000;
config.sigmaInit = 1 / sqrt(config.dimension);
config.sigmaStop = 1e-10;
config.sigmaMax = Inf;
config.centroidDiffEpsilonAbs = 1e-9;
config.centroidDiffEpsilonRel = 1e-9;
config.centroidDiffGenerations = 10;
config.upperBound = Inf;
config.distance = 'it';

[x, info] = lcCMSAES(@(x) c' * x, A, b, config);

disp('x'), disp(x);
disp('info'), disp(info);
fBest = fmin;
if abs(fBest) < 1e-6
  disp('abs(fBest - (f of best point))');
  disp(sprintf('%.20e', abs(fBest - info.aBest.f)));
else
  disp('abs(fBest - (f of best point)) / fBest');
  disp(sprintf('%.20e', abs(fBest - info.aBest.f) / fBest));
end

