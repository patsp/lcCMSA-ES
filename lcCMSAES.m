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

function [x, info] = lcCMSAES(f, A, b, config)
  % [x, info] = lcCMSAES(f, A, b, config)
  %    Attempts to minimize f(x)
  %    such that Ax = b and x >= 0
  %
  %    A is expected to be a (K x D)-matrix.
  %    b is expected to be a column vector of size K.
  %    f is expected to be a handle to a fitness function.
  %    config can contain configuration parameters.
  %
  %    The idea is to use an evolution strategy ((\mu/\mu_I, \lambda)-CMSA-ES)
  %    for optimization.
  %    First a solution to the linear system is computed.
  %    If there is no solution or one solution we are done.
  %    If there are multiple solutions the idea is to move
  %    inside the feasible region. For this, mutation is performed
  %    in the null space because Ay = 0 for y \in null(A).
  %    Projection is performed to make infeasible offspring feasible.

  dimension = size(A, 2);
  x = zeros(dimension, 1);

  info = struct;
  info.hasSolution = true;
  info.isSolutionUnique = false;
  info.numFitnessEvaluations = 0;
  info.numGenerations = 0;
  info.fDynamics = [];
  info.sigmaDynamics = [];

  if (~(config.lambda >= config.mu))
    error('config.lambda must be greater or equal to config.mu');
  end

  gStop = config.gStop;
  sigma = config.sigmaInit;
  sigmaStop = config.sigmaStop;
  sigmaMax = config.sigmaMax;
  mu = config.mu;
  lambda = config.lambda;
  targetCConditionNumber = 1e12;
  centroidDiffEpsilonAbs = config.centroidDiffEpsilonAbs;
  centroidDiffEpsilonRel = config.centroidDiffEpsilonRel;
  G = config.centroidDiffGenerations;
  upperBound = config.upperBound;

  x = A \ b;
  if ~(norm(abs(A * x - b)) < 1e-6)
    warning('~(norm(abs(A * x - b)) < 1e-6');
    warning(sprintf('~(%f < %f)', norm(abs(A * x - b)), 1e-6));
    %info.hasSolution = false;
    %return;
  end

  baseOfNullSpace = null(A);
  dimensionOfNullSpace = size(baseOfNullSpace, 2);
  if (dimensionOfNullSpace == 0)
    info.isSolutionUnique = true;
    return;
  end

  tau = 1 / sqrt(2 * dimensionOfNullSpace);
  tauC = 1 + (dimensionOfNullSpace * (dimensionOfNullSpace - 1)) / (2 * mu);
  covUpdate = min(floor(tauC), floor(dimensionOfNullSpace / 2));
  gLag = 50 * dimensionOfNullSpace;

  B = baseOfNullSpace;
  BInv = B';

  project = @(x, A, b, xInit, upperBound) ...
             projectToPositiveOrthantL1(x, A, b);
  if isfield(config, 'distance')
    if strcmp(config.distance, 'l1')
      project = @(x, A, b, xInit, upperBound) ...
                 projectToPositiveOrthantL1(x, A, b);
    elseif strcmp(config.distance, 'l2')
      project = @(x, A, b, xInit, upperBound) ...
                 projectToPositiveOrthantL2Squared(x, A, b, xInit, upperBound);
    elseif strcmp(config.distance, 'it')
      nPOpts = dimensionOfNullSpace * 10;
      pOpts = initReferencePointsForIterativeProjection(x, nPOpts, A, b);
      project = @(x, A, b, xInit, upperBound) ...
                 projectToPositiveOrthantIter(...
                             x, A, b, pOpts(:, randi([1, nPOpts], 1, 1)));
    else
      error(sprintf('Unknown projection method %s', config.distance));
    end
  end

  if isfield(config, 'xInit')
    x = config.xInit;
  else
    x = x + (B * norm(x) * randn(dimensionOfNullSpace, 1));
  end

  if (min(x) < -1e-20)
    x = project(x, A, b, x, upperBound);
  end

  a = struct;
  a.f = f(x);
  info.numFitnessEvaluations = info.numFitnessEvaluations + 1;
  a.x = x;
  a.z = zeros(dimension, 1);
  a.s = zeros(dimensionOfNullSpace, 1);
  a.sigma = sigma;
  aBest = a;
  aBestG = 0;

  xMemory = cell(G, 1);
  xMemory{1} = x;
  C = eye(dimensionOfNullSpace);
  sqrtC = eye(dimensionOfNullSpace);
  g = 0;
  do
    info.fDynamics(end + 1) = a.f;
    info.sigmaDynamics(end + 1) = a.sigma;

    offsprings = cell(lambda, 1);
    fitnessesAndSigmas = zeros(lambda, 2);
    for k = 1:lambda
      offspring.sigma = sigma * e^(tau * randn(1, 1));
      offspring.s = sqrtC * randn(dimensionOfNullSpace, 1);
      offspring.z = offspring.sigma * B * offspring.s;
      offspring.x = x + offspring.z;
      if (min(offspring.x) < -1e-20)
        prevOffspringZ = offspring.z;
        prevOffspringS = offspring.s;
        offspring.x = project(offspring.x, A, b, x, upperBound);
        offspring.z = offspring.x - x;
        offspring.s = (BInv * offspring.z) / offspring.sigma;
      end
      offspring.f = f(offspring.x);
      info.numFitnessEvaluations = info.numFitnessEvaluations + 1;

      fitnessesAndSigmas(k, 1) = offspring.f;
      fitnessesAndSigmas(k, 2) = offspring.sigma;

      if (offspring.f < aBest.f)
        aBest = offspring;
        aBestG = g + 1;
      end

      offsprings{k} = offspring;
    end
    [~, sortedIndices] = sortrows(fitnessesAndSigmas, 1);%[1, 2]);

    weight = 1 / mu;
    zCentroid = zeros(dimension, 1);
    sigmaCentroid = 0;
    ssCentroid = zeros(dimensionOfNullSpace, dimensionOfNullSpace);
    for k = 1:mu
      zCentroid = zCentroid + weight * offsprings{sortedIndices(k)}.z;
      sigmaCentroid = sigmaCentroid + ...
                      weight * offsprings{sortedIndices(k)}.sigma;
      s = offsprings{sortedIndices(k)}.s;
      ssCentroid = ssCentroid + weight * (s * s');
    end

    x = x + zCentroid;
    if (sigmaCentroid <= sigmaMax)
      sigma = sigmaCentroid;
    else
      sigma = sigmaMax;
    end
    C = (1 - (1 / tauC)) * C + (1 / tauC) * ssCentroid;

    if mod(g, covUpdate) == 0 || g == 0
      sqrtC = computeSqrtCNormalized(C, targetCConditionNumber);
      C = sqrtC * sqrtC';
      if iscomplex(sqrtC) || any(any(isnan(sqrtC) | isinf(sqrtC))) || ...
         iscomplex(C) || any(any(isnan(C) | isinf(C)))
        warning('stopped because C and sqrtC became unstable');
        break;
      end
    end

    % +1 because 1-based indexing
    xMemory{mod(g, G + 1) + 1} = x;

    centroidDiff = Inf;
    relativeCentroidDiff = Inf;
    % +1 because 1-based indexing
    if (g - G + 1 >= 1)
      centroidDiff = norm(xMemory{mod(g, G + 1) + 1} - ...
                          xMemory{mod(g - G, G + 1) + 1});
      relativeCentroidDiff = abs(norm(xMemory{mod(g, G + 1) + 1}) / ...
                                 norm(xMemory{mod(g - G, G + 1) + 1}) - 1);
      centroidDiff = -Inf;
      relativeCentroidDiff = -Inf;
      for k = 2:G
        centroidDiff = max(centroidDiff, norm(xMemory{k} - xMemory{k - 1}));
        relativeCentroidDiff = ...
          max(centroidDiff, abs(norm(xMemory{k}) / norm(xMemory{k - 1}) - 1));
      end
    end

    a.f = f(x);
    info.numFitnessEvaluations = info.numFitnessEvaluations + 1;
    a.x = x;
    a.z = zCentroid;
    a.sigma = sigma;

    g = g + 1;
  until (g > gStop || sigma < sigmaStop || ...
         centroidDiff < centroidDiffEpsilonAbs || ...
         relativeCentroidDiff < centroidDiffEpsilonRel || ...
         g - aBestG >= gLag)

  info.terminationCriterionStr = '';
  info.terminationCriterion = 0;
  if iscomplex(sqrtC)
    info.terminationCriterionStr = ['Covariance matrix became complex'];
    info.terminationCriterion = 0;
  elseif g > gStop
    info.terminationCriterionStr = ...
      ['Maximum number of generations reached, g > gStop ', ...
       sprintf('(%d > %d)', g, gStop)];
    info.terminationCriterion = 1;
  elseif sigma < sigmaStop
    info.terminationCriterionStr = ...
      ['Mutation step size decreased below threshold,', ...
       ' sigma < sigmaStop ', ...
       sprintf('(%f < %f)', sigma, sigmaStop)];
    info.terminationCriterion = 2;
  elseif centroidDiff < centroidDiffEpsilonAbs
    info.terminationCriterionStr = ...
      ['Absolute difference between centroids', ...
       ' decreased below threshold,', ...
       ' centroidDiff < centroidDiffEpsilonAbs ', ...
       sprintf('(%f < %f)', centroidDiff, centroidDiffEpsilonAbs)];
    info.terminationCriterion = 3;
  elseif relativeCentroidDiff < centroidDiffEpsilonRel
    info.terminationCriterionStr = ...
      ['Relative difference between centroids', ...
       ' decreased below threshold,', ...
       ' relativeCentroidDiff < centroidDiffEpsilonRel ', ...
       sprintf('(%f < %f)', relativeCentroidDiff, centroidDiffEpsilonRel)];
    info.terminationCriterion = 4;
  elseif g - aBestG >= gLag
    info.terminationCriterionStr = ...
      ['Best-so-far has not been updated for the last gLag generations ', ...
       sprintf('(gLag = %d)', gLag)];
    info.terminationCriterion = 5;
  end

  info.aBest = aBest;
  info.numGenerations = g;

end
