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

function [fBest, aBest] = runlcCMSAESOnRotatedKleeMintyProblem(problem, ...
                                                               budget, ...
                                                               lbnds, ...
                                                               ubnds, ...
                                                               input)

  config.dimension = input.dim;
  assert(config.dimension == length(lbnds));
  assert(config.dimension == length(ubnds));

  config.lambda = 4 * config.dimension;
  config.mu = floor(config.lambda / 4);
  % numFunctionEvaluations = numGenerations * config.lambda + numGenerations + 1
  %                        = (numGenerations + 1) * config.lambda + 1
  % numGenerations = floor((numFunctionEvaluations - 1) / config.lambda) - 1
  config.gStop = floor((budget - 1) / config.lambda) - 1;
  config.sigmaInit = 1 / sqrt(config.dimension);
  config.sigmaStop = 1e-10;
  config.sigmaMax = 1e4;
  config.centroidDiffEpsilonAbs = 1e-9;
  config.centroidDiffEpsilonRel = 1e-9;
  config.centroidDiffGenerations = 10;
  config.upperBound = Inf;
  config.distance = 'l2';

  [A, b, f] = preprocessCocoProblem(@(x) fHelper(x, problem), ...
                                    @(x) gHelper(x, problem), ...
                                    lbnds(:), ubnds(:));
  config.dimension = size(A, 2);
  [x, info] = ...
  lcCMSAES(f, A, b, config);
  %disp('info.terminationCriterionStr'), disp(info.terminationCriterionStr);
  %disp('info.aBest.f'), disp(info.aBest.f);
  fBest = info.aBest.f;
  aBest = struct();
  DPrime = size(lbnds(:), 1);
  KPrime = size(gHelper(zeros(size(lbnds(:), 1), 1), problem), 1);
  backtransform = @(x) x(1:DPrime) + lbnds(:);
  xBestBacktransformed = backtransform(info.aBest.x);
  aBest.y = xBestBacktransformed;
  aBest.val = info.aBest.f;
  gv = gHelper(xBestBacktransformed, ...
               problem);
  aBest.conv = sum(gv .* (gv > 0));
end

function f = fHelper(x, problem)
  [f, ~] = evaluateRotatedKleeMintyProblem(x, problem);
end

function g = gHelper(x, problem)
  [~, g] = evaluateRotatedKleeMintyProblem(x, problem);
end

