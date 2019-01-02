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

function [A, b, f] = preprocessCocoProblem(fPrime, ...
                                           constraintFun, ...
                                           lbnds, ubnds)
  % [A, b, f] = preprocessCocoProblem(fPrime, constraintFun, ...
  %                                   lbnds, ubnds)
  %    This function pre-processes a BBOB COCO bbob-constrained suite
  %    optimization problem of the form
  %        fPrime(x) -> min!
  %        s.t. constraintFun(x) <= 0
  %             lbnds <= x <= ubnds
  %    into a problem of the form
  %        f(x) -> min!
  %        s.t. Ax = b
  %             x >= 0.

  DPrime = size(lbnds, 1);
  KPrime = size(constraintFun(zeros(size(lbnds, 1), 1)), 1);
  hdummy = @(x) 0;
  [Aineq, bineq, Aeq, beq] = ...
      approximateConstraintsLocallyAsLinearConstraints(zeros(DPrime, 1), ...
                                                  1, constraintFun, hdummy);
  %[A, b] = transformToStandardFormIneq(Aineq, bineq, lbnds, ubnds);
  %f = @(x) fPrime(backtransformStandardFormVector(x, DPrime, KPrime));

  % This is an alternative to transformFormIneq and
  % backtransformStandardFormVector to have less variables
  % in the transformed system.
  % This makes use of the fact that there are lower and upper bounds
  % for all the variables. Then, every variable x_i in the original
  % system is replaced by x'_i + lbnds_i, then x'_i >= 0 implies
  % x_i - lbnds_i >= 0 which implies further that x_i >= lnbds_i.
  % Further a constraint is introduced for handling the upper bound:
  % for every x_i a constraint x_i - lbnds_i <= ubnds_i - lbnds_i is introduced
  % that yields x_i <= ubnds_i.
  % The following code puts this into matrix form.
  W = Aineq;
  c = bineq;
  % Special case for the case that we already have standard bounds
  % (note: this interprets numbers larger than 1e20 as infinity).
  if all(lbnds == 0) && all(ubnds >= 1e20)
    A = [W eye(KPrime, KPrime)];
    b = [c - W * lbnds];
  else
    A = [W                    eye(KPrime, KPrime),   zeros(KPrime, DPrime);
         eye(DPrime, DPrime), zeros(DPrime, KPrime), eye(DPrime, DPrime)];
    b = [c - W * lbnds;
         ubnds - lbnds];
  end
  % Backtransformation is done by adding lbnds_i to every x'_i because
  % x_i = x'_i + lbnds_i.
  f = @(x) fPrime(x(1:DPrime) + lbnds);
end

