
% =========
% ndPH software, 
% v0.1 
% 
% Copyright (C) 2016, 2017 Ben Cassidy
% released under GPL license
% =========


function [tPermStr] = NDPH_compare_AUtCtest( in_1, in_2 , pairedORunpaired, nperm, q, DO_PLOT)

nd_1 = in_1.funcFit.yhatfd_D1_all_alt;
nd_2 = in_2.funcFit.yhatfd_D1_all_alt;

args = in_1.funcFit.evalarg_all;
assert(all(args == in_2.funcFit.evalarg_all))

NTHRESH = 21;
STDTHRESH_IND = 11;
% we will use 0.001 threshold, but examine around for assessing stability
meanTHRESHrng = logspace(-4,-2,NTHRESH); 

switch lower(pairedORunpaired)
    case 'paired'
        tPermStr = paired_tperm_fd(nd_1, nd_2, nperm, q, args, meanTHRESHrng, STDTHRESH_IND);
    case 'unpaired'
        tPermStr = unpaired_tperm_fd(nd_1, nd_2, nperm, q, args, meanTHRESHrng, STDTHRESH_IND);
    otherwise
        error('specify paired or unpaired t test')     
end

if DO_PLOT
    AUC = tPermStr.AUCstr.AUC_full(STDTHRESH_IND);
    validInds = tPermStr.AUCstr.validInds_full(:,STDTHRESH_IND);
    plot_ndph_test(tPermStr, AUC, validInds )
    AUC_null = tPermStr.AUCstr.AUC_null(:,STDTHRESH_IND);
    nullthresh = tPermStr.AUCstr.AUC_nullThresh;
    plot_AUC(AUC_null, nullthresh, AUC)
end


end

function tpermStr = paired_tperm_fd(x1fd, x2fd, nperm, q, argvals, meanTHRESHrng, STDTHRESH_IND)
% This function partially based on tperm_fd.m by Prof. J. Ramsay
% see http://www.psych.mcgill.ca/misc/fda/

%  check first two arguments

if nargin < 2
    error('Arguments X1FD and X2FD not supplied.');
end

if ~isa_fd(x1fd) || ~isa_fd(x2fd)
    error('x1fd and x2fd must both be functional data objects');
end

if length(size(getcoef(x1fd))) > 2 || length(size(getcoef(x2fd))) > 2
    error('Both of X1FD and X2FD are not univariate.');
end

%  Set default arguments

if nargin < 5,  argvals = [];  end
if nargin < 4,  q = 0.05;      end
if nargin < 3,  nperm = 200;   end

range1 = getbasisrange(getbasis(x1fd));
range2 = getbasisrange(getbasis(x2fd));


if ~all(range1 == range2)
    error('x1fd and x2fd do not have the same range.');
end

if isempty(argvals)
    narg = 101;
    argvals = linspace(range1(1),range1(2),narg)';
else
    narg = length(argvals);
end

q = 1-q;

x1mat = eval_fd(argvals,x1fd);
x2mat = eval_fd(argvals,x2fd);

n1 = size(x1mat,2);
n2 = size(x2mat,2);
assert(n1==n2);

Tnull = zeros(nperm,1);

Tnullvals = zeros(length(argvals),nperm);

NTHRESH = length(meanTHRESHrng);


% smooth slightly to prevent numerical precision errors when assessing
% non-normality
% this is only a heuristic for thresholding, empirically robust to choice.
windowSizeData = ceil(narg./25); 
bD = (1/windowSizeData)*ones(1,windowSizeData); 
smoothD1_1 = filtfilt(bD,1, x1mat);
smoothD1_2 = filtfilt(bD,1, x2mat);

[h_1, ADStat_1, CV_1] = myADtestNormal(smoothD1_1', 0.01);
[h_2, ADStat_2, CV_2] = myADtestNormal(smoothD1_2', 0.01);
NONNORMAL_INDICATOR = ~(~h_1 | ~h_2);

AUC_null = zeros(nperm,NTHRESH);
for i=1:nperm
    doPerm = (randn(n1,1) > 0);
    
    x1matP = zeros(size(x1mat));
    x2matP = zeros(size(x2mat));
    x1matP(:,~doPerm) = x1mat(:,~doPerm);
    x1matP(:,doPerm) = x2mat(:,doPerm);
    x2matP(:,~doPerm) = x2mat(:,~doPerm);
    x2matP(:,doPerm) = x1mat(:,doPerm);
    
    XmatPDiff = x2matP - x1matP;
    
    tPmean = mean(XmatPDiff,2);
    tPvar = var(XmatPDiff,0,2) ./ n1;
    Tnullvals(:,i) = abs(tPmean./sqrt(tPvar));
    Tnull(i) = max(Tnullvals(:,i));
    
    [AUCstrNull] = fanova_tperm_AUC(Tnullvals(:,i), tPmean, meanTHRESHrng, NONNORMAL_INDICATOR);
    AUC_null(i,:) = AUCstrNull.AUC_full;
end


XmatDiff = x2mat - x1mat;
tmean = mean(XmatDiff,2);
tvar = var(XmatDiff,0,2) ./ n1;

Tvals = abs(tmean)./sqrt(tvar);

% [AUCstr] = fanova_paired_tperm_AUC(Tvals, tmean, meanTHRESHrng, smoothD1_1, smoothD1_2);
[AUCstr] = fanova_tperm_AUC(Tvals, tmean, meanTHRESHrng, NONNORMAL_INDICATOR);
AUC_nullMean = mean(AUC_null(:,STDTHRESH_IND));
AUC_nullStd = std(AUC_null(:,STDTHRESH_IND));
AUCstr.AUC_nullThresh = norminv(q, AUC_nullMean, AUC_nullStd);
AUCstr.AUC_null = AUC_null;
AUCstr.AUC_nullMean = AUC_nullMean;
AUCstr.AUC_nullStd = AUC_nullStd;

Tobs = max(Tvals);

pval = mean( Tobs < Tnull );
qval = quantile(Tnull, q(1));

pvals_pts = zeros(narg,1);
qvals_pts = zeros(narg,1);
for i=1:narg
    pvals_pts(i) = mean(Tvals(i) < Tnullvals(i,:));
    qvals_pts(i) = quantile(Tnullvals(i,:), q(1));
end

tpermStr.tmean = tmean;
tpermStr.tvar = tvar;

tpermStr.q            = q;
tpermStr.pval         = pval;
tpermStr.qval         = qval;
tpermStr.Tobs         = Tobs;
tpermStr.Tnull        = Tnull;
tpermStr.Tvals        = Tvals;
tpermStr.Tnullvals    = Tnullvals;
tpermStr.pvals_pts    = pvals_pts;
tpermStr.qvals_pts    = qvals_pts;
tpermStr.argvals      = argvals;

tpermStr.XmatDiff = XmatDiff;
tpermStr.x1mat = x1mat;
tpermStr.x2mat = x2mat;

tpermStr.AUCstr = AUCstr;


end

function tpermStr = unpaired_tperm_fd(x1fd, x2fd, nperm, q, argvals, meanTHRESHrng, STDTHRESH_IND)
% This function partially based on tperm_fd.m by Prof. J. Ramsay
% see http://www.psych.mcgill.ca/misc/fda/

if nargin < 2
    error('Arguments X1FD and X2FD not supplied.');
end

if ~isa_fd(x1fd) || ~isa_fd(x2fd)
    error('x1fd and x2fd must both be functional data objects');
end

if length(size(getcoef(x1fd))) > 2 || length(size(getcoef(x2fd))) > 2
    error('Both of X1FD and X2FD are not univariate.');
end

%  Set default arguments

% if nargin < 6,  plotres = 1;   end
if nargin < 5,  argvals = [];  end
if nargin < 4,  q = 0.05;      end
if nargin < 3,  nperm = 200;   end


range1 = getbasisrange(getbasis(x1fd));
range2 = getbasisrange(getbasis(x2fd));


if ~all(range1 == range2)
    error('x1fd and x2fd do not have the same range.');
end

if isempty(argvals)
    narg = 101;
    argvals = linspace(range1(1),range1(2),narg)';
else
    narg = length(argvals);
end

q = 1-q;

x1mat = eval_fd(argvals,x1fd);
x2mat = eval_fd(argvals,x2fd);

Xmat = [x1mat,x2mat];

n1 = size(x1mat,2);
n2 = size(x2mat,2);

Tnull = zeros(nperm,1);

Tnullvals = zeros(length(argvals),nperm);


NTHRESH = length(meanTHRESHrng);
% smooth slightly to prevent numerical precision errors
windowSizeData = ceil(narg./25); % this is only a heuristic for thresholding, fairly invariant to choice.
bD = (1/windowSizeData)*ones(1,windowSizeData); 
smoothD1_1 = filtfilt(bD,1, x1mat);
smoothD1_2 = filtfilt(bD,1, x2mat);

[h_1, ADStat_1, CV_1] = myADtestNormal(smoothD1_1', 0.01);
[h_2, ADStat_2, CV_2] = myADtestNormal(smoothD1_2', 0.01);
NONNORMAL_INDICATOR = ~(~h_1 | ~h_2);

% permInit = [true(n1,1); false(n2,1)];

AUC_null = zeros(nperm,NTHRESH);
for i = 1:nperm
    
    permOrder = randperm(n1+n2);
   
    x1matP = Xmat(:,permOrder(1:n1));
    x2matP = Xmat(:,permOrder(n1+1:end));

    tPmean1 = mean(x1matP,2);
    tPmean2 = mean(x2matP,2);
    
    tPvar1 = var(x1matP, 0, 2) ./ n1;
    tPvar2 = var(x2matP, 0, 2) ./ n2;
    
    Tnullvals(:,i) = abs(tPmean1-tPmean2)./sqrt(tPvar1+tPvar2);
    Tnull(i) = max(Tnullvals(:,i));
    
%     [AUCstrNull] = fanova_paired_tperm_AUC(Tnullvals(:,i), [tPmean1, tPmean2], meanTHRESHrng, x1matP_sm, x2matP_sm);
    [AUCstrNull] = fanova_tperm_AUC(Tnullvals(:,i), [tPmean1, tPmean2], meanTHRESHrng, NONNORMAL_INDICATOR);
    AUC_null(i,:) = AUCstrNull.AUC_full;

end

tmean1 = mean(x1mat,2);
tmean2 = mean(x2mat,2);

tvar1 = var(x1mat, 0, 2) ./ n1;
tvar2 = var(x2mat, 0, 2) ./ n2;

Tvals = abs(tmean1-tmean2)./sqrt(tvar1+tvar2);

% [AUCstr] = fanova_tperm_AUC(Tvals, [tmean1, tmean2], meanTHRESHrng, NONNORMAL_INDICATOR);
% IMPORTANT : :: :: :  STABILITY CORRECT ON ABS(DIFF(MEANS)) RATHER THAN
% JUST INTERSECTION OF CONDITIONS ON EACH MEAN. both answers are fine, but
% this is a better treatment of the model validity condition. valid when
% abs(difference of means) is above a threshold, rather than just when both
% means are above. slightly small AUC as a result so this is more
% conservative for hypothesis testing. 
[AUCstr] = fanova_tperm_AUC(Tvals, abs(tmean1-tmean2), meanTHRESHrng, NONNORMAL_INDICATOR); 
AUC_nullMean = mean(AUC_null(:,STDTHRESH_IND));
AUC_nullStd = std(AUC_null(:,STDTHRESH_IND));
AUCstr.AUC_nullThresh = norminv(q, AUC_nullMean, AUC_nullStd);
AUCstr.AUC_null = AUC_null;
AUCstr.AUC_nullMean = AUC_nullMean;
AUCstr.AUC_nullStd = AUC_nullStd;

Tobs  = max(Tvals);

pval = mean( Tobs < Tnull );
qval = quantile(Tnull, q(1));

pvals_pts = zeros(narg,1);
qvals_pts = pvals_pts;
for i=1:narg
    pvals_pts(i) = mean(Tvals(i) < Tnullvals(i,:));
    qvals_pts(i) = quantile(Tnullvals(i,:), q(1));
end


tpermStr.tmean1 = tmean1;
tpermStr.tmean2 = tmean2;
tpermStr.tvar1 = tvar1;
tpermStr.tvar2 = tvar2;

tpermStr.q            = q;
tpermStr.pval         = pval;
tpermStr.qval         = qval;
tpermStr.Tobs         = Tobs;
tpermStr.Tnull        = Tnull;
tpermStr.Tvals        = Tvals;
tpermStr.Tnullvals    = Tnullvals;
tpermStr.pvals_pts    = pvals_pts;
tpermStr.qvals_pts    = qvals_pts;
tpermStr.argvals      = argvals;

% tpermStr.XmatDiff = XmatDiff;
tpermStr.x1mat = x1mat;
tpermStr.x2mat = x2mat;

tpermStr.AUCstr = AUCstr;

end

function [AUCstr] = fanova_tperm_AUC(Tvals, tmean, meanTHRESHrng, NONNORMAL_INDICATOR)

NTHRESH = length(meanTHRESHrng);
AUC_full = zeros(NTHRESH,1);

% non-normality condition - this is only for assessing validity of the 
% test statistic through violation of the Gaussian assumption and does not 
% affect the threshold selection. 

windowSizeMean = ceil(numel(Tvals)./250); % this is only a heuristic for thresholding, 
bM = (1/windowSizeMean)*ones(1,windowSizeMean); % this is only a heuristic for thresholding, 
smoothAbsMean = ( filtfilt(bM,1, abs(tmean) ) ); % this is only a heuristic for thresholding, 

validInds_full = NaN*zeros(numel(Tvals), NTHRESH);
for lp2=1:NTHRESH
    nullDomainMean =  (smoothAbsMean < repmat( meanTHRESHrng(lp2) * max(smoothAbsMean), size(smoothAbsMean,1),1 ));    
    validInds_full(:,lp2) = ( ~(all(nullDomainMean,2)) & ~NONNORMAL_INDICATOR);
    AUC_full(lp2) = trapz(Tvals(validInds_full(:,lp2)>0)) ./ numel(Tvals);
end

AUCstr.AUC_full                 = AUC_full;
AUCstr.validInds_full           = validInds_full;
AUCstr.NONNORMAL_INDICATOR      = NONNORMAL_INDICATOR;

end


function [H, ADStat, CV] = myADtestNormal(x, alpha)

% BC - Trimmed and vectorized version of matlab adtest for normal distribution

n = size(x,1);

z = normcdf(x, repmat(mean(x),n,1), repmat(std(x),n,1));
ADStat = ComputeADStat(z,n); % 234

[alphas, CVs] = computeCriticalValues_norm(n);

% 1-D interpolation into the tabulated results. In the upper tail,
% CV vs log(alpha) is approximately linear
pp = pchip(log(alphas), CVs);
CV = ppval(pp,log(alpha));

H = (ADStat > CV)';
%------------------------------------------
function ADStat = ComputeADStat(z,n)
% Compute Anderson-Darling Statistic
% Sort the data and compute the statistic
z = sort(z);
w = 2*(1:n) - 1;
ADStat = -w*(log(z)+ log(1-z(end:-1:1,:)))./n - n;

end
%------------------------------------------
function [alphas, CVs] = computeCriticalValues_norm(n)

alphas = [0.0005 0.001 0.0015 0.002 0.005 0.01 0.025 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5...
    0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99];

% An improved version of the Petitt method for the composite normal case.
% The model used is A_n = A_inf (1+b_0/n+b_1/n^2), where different
% estimates b_0 and b_1 are used for each significance level. This allows
% us to model the entire range of the distribution.
CVs = ...
    [ 1.5649  1.4407  1.3699  1.3187  1.1556    1.0339    0.8733    0.7519    0.6308    0.5598    0.5092    0.4694    0.4366    0.4084...
    0.3835    0.3611    0.3405    0.3212    0.3029    0.2852    0.2679    0.2506    0.2330    0.2144...
    0.1935    0.1673    0.1296] +...
    [-0.9362 -0.9029  -0.8906  -0.8865  -0.8375   -0.7835   -0.6746   -0.5835   -0.4775   -0.4094   -0.3679   -0.3327   -0.3099   -0.2969...
    -0.2795   -0.2623   -0.2464   -0.2325   -0.2164   -0.1994   -0.1784   -0.1569   -0.1377   -0.1201...
    -0.0989   -0.0800   -0.0598]./n +...
    [-8.3249  -6.6022 -5.6461  -4.9685  -3.2208   -2.1647   -1.2460   -0.7803   -0.4627   -0.3672   -0.2833   -0.2349   -0.1442   -0.0229...
    0.0377    0.0817    0.1150    0.1583    0.1801    0.1887    0.1695    0.1513    0.1533    0.1724...
    0.2027    0.3158    0.6431]./n^2;
end
end

function plot_ndph_test(tpermStr, AUC, validInds )

Tvals_valid = tpermStr.Tvals; Tvals_valid(~validInds) = NaN;
pltsty = {'b-', 'b--'};
figure();
ylim = [ min([tpermStr.Tvals;tpermStr.qvals_pts]),max([tpermStr.Tobs;tpermStr.qval])];
plot(NaN,NaN, pltsty{1}, 'linewidth', 2)  % get overlapping correct order
hold on
plot(NaN,NaN, pltsty{2}) 
plot(tpermStr.argvals, tpermStr.Tvals, pltsty{2});
plot(tpermStr.argvals, Tvals_valid, pltsty{1}, 'linewidth', 2 );
ahArea = area(tpermStr.argvals, Tvals_valid);
ahArea.FaceColor = [1 0 0];
ahArea.FaceAlpha = 0.2;
ah = gca();
ah.XLabel.String = 'Filtration Distance';
ah.YLabel.String = 'Functional Difference';
ah.FontSize = 22;
ah.Title.String = ['ndPH test ' '\color{red}AUC= ' num2str(AUC)];

lh = legend('Observed functional difference', 'Invalid');
lh.Location = 'NorthWest';

end

function plot_AUC(AUC_null, nullthresh, AUC)

[nHistCounts, histEdges] = histcounts(AUC_null, 'Normalization', 'probability');
histEdges = sort([histEdges, nullthresh]);
figure()
histogram(AUC_null, histEdges,  'Normalization', 'probability');
hold on

ah = gca();

line([nullthresh,nullthresh], [0 ah.YLim(2)],'color', 'r', 'linewidth', 2);
scatter(AUC, 0, 300, '^r', 'filled')

ah = gca();
ah.FontSize = 22;
ah.XLabel.String = 'ndPH area under curve';
ah.YLabel.String = 'Density';

fh = gcf();
fh.Color = 'w';
end



