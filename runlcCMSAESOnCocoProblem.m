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

function runlcCMSAESOnCocoProblem(...
                                   problem, lbnds, ubnds, budget)
  % runlcCMSAESOnCocoProblem(...
  %                          problem, lbnds, ubnds, budget)
  %    This functions runs the lcCMSAES with the given budget on
  %    the given COCO problem with lower bounds lbnds and upper bounds
  %    ubnds.

  config.dimension = length(lbnds);
  config.lambda = 4 * config.dimension;
  config.mu = floor(config.lambda / 4);
  % Set gStop according to budget. Note that this is not exact
  % because the pre-processing below also takes up some budget.
  % This is okay for us because in the order of 10^5 a few hundred
  % more or less is not an issue.
  % numFunctionEvaluations = numGenerations * config.lambda + numGenerations + 1
  %                        = (numGenerations + 1) * config.lambda + 1
  % numGenerations = floor((numFunctionEvaluations - 1) / config.lambda) - 1
  config.gStop = floor((budget - 1) / config.lambda) - 1;
  config.sigmaInit = 1 / sqrt(config.dimension);
  config.sigmaStop = 1e-6;
  config.sigmaMax = 1e4;
  config.centroidDiffEpsilonAbs = 1e-9;
  config.centroidDiffEpsilonRel = 1e-9;
  config.centroidDiffGenerations = 10;
  config.upperBound = Inf;
  config.distance = 'it';

  [A, b, f] = preprocessCocoProblem(@(x) fWrapper(problem, x), ...
                                    @(x) gWrapper(problem, x), ...
                                    lbnds(:), ubnds(:));
  config.dimension = size(A, 2);
  [x, info] = lcCMSAES(f, A, b, config);
  disp('info.terminationCriterionStr'), disp(info.terminationCriterionStr);
  disp('info.aBest.f'), disp(info.aBest.f);
end

function y = fWrapper(problem, x)
  y = cocoEvaluateFunction(problem, x(:));
end

function y = gWrapper(problem, x)
  y = [];
  if cocoProblemGetNumberOfConstraints(problem) > 0
    y = cocoEvaluateConstraint(problem, x(:));
  end
  y = y(:);
end

