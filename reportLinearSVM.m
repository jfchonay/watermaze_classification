function S = reportLinearSVM(M, featureNames, doCrossVal, topK)
% REPORTLINEARSVM  Summarize a linear fitcsvm() model (version-robust).
% Usage:
%   S = reportLinearSVM(M)
%   S = reportLinearSVM(M, featureNames)
%   S = reportLinearSVM(M, featureNames, true)
%   S = reportLinearSVM(M, featureNames, true, 20)
%
% Inputs
%   M             : ClassificationSVM from fitcsvm (Linear kernel)
%   featureNames  : cellstr of feature names (optional)
%   doCrossVal    : logical, run crossval(M) to get kfoldLoss (default: false)
%   topK          : integer, how many top-|w| features to show (default: 15)

if nargin < 2 || isempty(featureNames), featureNames = []; end
if nargin < 3 || isempty(doCrossVal),   doCrossVal   = false; end
if nargin < 4 || isempty(topK),         topK         = 15; end

%% Kernel sanity check
if isprop(M,'KernelParameters') && isfield(M.KernelParameters,'Function')
    kern = M.KernelParameters.Function;
else
    kern = 'linear'; % best guess
end
if ~strcmpi(kern,'linear')
    warning('Model kernel is %s, but this reporter assumes linear.', kern);
end

%% Box constraint (C) and solver details
Cvec = M.BoxConstraints;
if isscalar(Cvec)
    C = Cvec;
else
    C = mode(Cvec); % representative if per-observation C used
end

Solver = '';
Iter   = NaN;
if isprop(M,'ConvergenceInfo') && ~isempty(M.ConvergenceInfo)
    CI = M.ConvergenceInfo;
    if isfield(CI,'Solver'),     Solver = CI.Solver;        end
    if isfield(CI,'Iterations'), Iter   = CI.Iterations;    end
end

%% Support vectors / indices (robust across versions)
SV = M.SupportVectors;                 % [nSV x d], always available
nSV_total = size(SV,1);

if isprop(M,'IsSupportVector') && ~isempty(M.IsSupportVector)
    svIdx = find(M.IsSupportVector);
elseif isprop(M,'SupportVectorIndices') && ~isempty(M.SupportVectorIndices)
    svIdx = M.SupportVectorIndices;
else
    svIdx = []; % not exposed in this release
end

%% Labels of support vectors
% Many releases expose SupportVectorLabels; if not, we try to recover via Alpha sign.
y_sv = [];
if isprop(M,'SupportVectorLabels') && ~isempty(M.SupportVectorLabels)
    y_sv = double(M.SupportVectorLabels(:));
end

alpha = M.Alpha(:);                    % [nSV x 1]
b     = M.Bias;

%% Compute weight vector w in the training space
% If Standardize=true, the optimization ran on Z=(X-Mu)./Sigma. Use the same space.
useStandardized = isprop(M,'Mu') && isprop(M,'Sigma') && ~isempty(M.Mu) && ~isempty(M.Sigma);

if useStandardized
    ZSV = bsxfun(@rdivide, bsxfun(@minus, SV, M.Mu), M.Sigma);
    if ~isempty(y_sv)
        w = ZSV' * (alpha .* y_sv);
    else
        % Fall back: sign of alpha is consistent with labels for binary SVMs
        w = ZSV' * alpha;
        warning('SupportVectorLabels not found; w computed using Alpha only.');
    end
else
    if ~isempty(y_sv)
        w = SV' * (alpha .* y_sv);
    else
        w = SV' * alpha;
        warning('SupportVectorLabels not found; w computed using Alpha only.');
    end
end

normw  = norm(w);
margin = 2 / normw;

%% Per-class SV counts (if binary and labels available)
nSV_A = NaN; nSV_B = NaN; classA = ""; classB = "";
if ~isempty(y_sv)
    classes = M.ClassNames;
    if numel(classes)==2
        classA = string(classes(1)); classB = string(classes(2));
        nSV_A = sum(y_sv == double(classes(1)));
        nSV_B = sum(y_sv == double(classes(2)));
    end
else
    classes = M.ClassNames;
end

%% Hyperparameter optimization snapshot (if present)
BestHyperparams = struct();
if isprop(M,'HyperparameterOptimizationResults') && ~isempty(M.HyperparameterOptimizationResults)
    try
        bo = M.HyperparameterOptimizationResults;
        BestHyperparams = struct(bo.XAtMinObjective);
    catch
        % tolerate structural differences across releases
    end
end

%% Optional cross-validation
CV = struct(); cvLoss = NaN;
if doCrossVal
    CVSVM = crossval(M);
    cvLoss = kfoldLoss(CVSVM);
    CV.K = CVSVM.KFold;
    CV.Loss = cvLoss;
end

%% Top features by |w|
absw = abs(w);
[abswSorted, idx] = sort(absw, 'descend');
k = min(topK, numel(w));
topIdx = idx(1:k);
topW   = w(topIdx);
if ~isempty(featureNames) && numel(featureNames)==numel(w)
    topNames = featureNames(topIdx);
else
    topNames = arrayfun(@(j) sprintf('x_%d', j), topIdx, 'UniformOutput', false);
end

%% Print compact report
fprintf('\n===== Linear SVM Report (fitcsvm) =====\n');
fprintf('Kernel                : %s\n', kern);
fprintf('Box Constraint (C)    : %g\n', C);
if ~isempty(Solver), fprintf('Solver                : %s\n', Solver); end
if ~isnan(Iter),     fprintf('Iterations            : %d\n', Iter);   end
fprintf('Standardized training : %s\n', string(useStandardized));
fprintf('Bias (b)              : %.6g\n', b);
fprintf('||w||                 : %.6g\n', normw);
fprintf('Geometric margin      : %.6g\n', margin);
fprintf('# Support Vectors     : %d\n', nSV_total);
if ~isempty(svIdx)
    fprintf('  (Indices available; first 10): %s\n', mat2str(svIdx(1:min(10,end))));
else
    fprintf('  (SV indices not exposed in this MATLAB release)\n');
end
if numel(classes)==2 && ~isnan(nSV_A)
    fprintf('  By class [%s]: %d,  [%s]: %d\n', classA, nSV_A, classB, nSV_B);
end
if ~isempty(fieldnames(BestHyperparams))
    fprintf('Bayes Opt best params :\n');
    fns = fieldnames(BestHyperparams);
    for i=1:numel(fns)
        val = BestHyperparams.(fns{i});
        if isnumeric(val), valStr = mat2str(val);
        elseif isstring(val) || ischar(val), valStr = char(val);
        elseif islogical(val), valStr = string(val);
        else, valStr = '<obj>';
        end
        fprintf('  - %s: %s\n', fns{i}, valStr);
    end
end
if doCrossVal
    fprintf('K-fold CV loss        : %.6g\n', cvLoss);
end
fprintf('\nTop %d features by |weight|:\n', k);
for i = 1:k
    fprintf('%3d) %-20s  w=%.6g  |w|=%.6g\n', i, topNames{i}, topW(i), abswSorted(i));
end
fprintf('=======================================\n\n');

%% Package outputs
S = struct();
S.C                    = C;
S.Kernel               = kern;
S.Solver               = Solver;
S.Iterations           = Iter;
S.Standardized         = useStandardized;
S.w                    = w(:);
S.b                    = b;
S.norm_w               = normw;
S.margin               = margin;
S.SupportVectors       = SV;
S.SupportVectorLabels  = y_sv;            % may be empty if not exposed
S.SupportVectorIndices = svIdx;           % robust fallback
S.nSV_total            = nSV_total;
S.nSV_by_class         = struct('class1', classA, 'n', nSV_A, ...
                                'class2', classB, 'm', nSV_B);
S.TopFeatureIdx        = topIdx(:);
S.TopFeatureNames      = topNames(:);
S.TopFeatureWeights    = topW(:);
S.AbsWeightsSorted     = abswSorted(:);
S.BestHyperparameters  = BestHyperparams;
S.CV                   = CV;
end
