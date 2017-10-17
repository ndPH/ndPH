
% =========
% ndPH software, 
% v0.1 
% 
% Copyright (C) 2016, 2017 Ben Cassidy
% released under GPL license
% =========

% distMat is a #Nodes x #Nodes x #Networks array, 
% lower triangular within each network
% distances in range [0, 1]

% [out ] = NDPH_prepare(distMats, [MST_ONLY], [CI_MAX_IGNORE])
% MST_ONLY -> do not do curve fitting
% CI_MAX_IGNORE -> this value corresponds to a link strength to ignore wrt
% finding structure in MST. typically should be ==1 for statistical networks.


function [out ] = NDPH_fit(distMats, varargin)

assert(nargin < 4)
if nargin < 2
    MST_ONLY = false;
else
    MST_ONLY = varargin{1};
    assert(islogical(MST_ONLY));
end
if nargin < 3
    CI_MAX_IGNORE = NaN;
else
    CI_MAX_IGNORE = varargin{2};
    assert(isnumeric(CI_MAX_IGNORE));
    assert(numel(CI_MAX_IGNORE)==1);
end

assert(isnumeric(distMats));
assert(ndims(distMats) ==2 || ndims(distMats) == 3);

[nNodes, nNodes2, nNets] = size(distMats);
assert(nNodes == nNodes2)

MST = zeros(nNodes-1, 3, nNets);
filtration = NaN*zeros(nNodes+1,nNets);

for lp = 1:nNets
    fprintf('.')
    NW = minspantree(graph(distMats(:,:,lp), 'lower', 'OmitSelfLoops'));
    MST(:,:,lp) = table2array( sortrows(NW.Edges, 2) );
%     $ SLM via   or similar
%     [SLM(:,:,lp), MST2(:,:,lp)] = single_linkage_matrix(distMats(:,:,lp));
    filtration(:,lp) = [0; MST(:,3,lp); 1];    
end

if ~(MST_ONLY) 
    
betti0 = repmat( (linspace(0,1,nNodes+1))',1,nNets); % normalize for numerical stability
betti0(1) = betti0(2); % cannot be zero connected components, limit to scaled of 1 connected component
betti0 = flip(betti0);

[yhatcoef_D0_all,yhatfd_D0_all, yhatfd_D1_all, yhatfd_D1_all_alt, ...
    evalarg_all, opts, opts_all] ...
    = PH_Betti0_functionFit_sparseDenseNetwork(filtration, betti0, 1);
D0_buf_all = eval_fd(evalarg_all, yhatfd_D0_all);

out.funcFit.yhatcoef_D0_all = yhatcoef_D0_all;
out.funcFit.yhatfd_D0_all = yhatfd_D0_all;
out.funcFit.yhatfd_D1_all = yhatfd_D1_all;
out.funcFit.yhatfd_D1_all_alt = yhatfd_D1_all_alt;
out.funcFit.evalarg_all = evalarg_all;
out.funcFit.opt = opts;
out.funcFit.opts_all = opts_all;

end

out.MST = MST;

end


function [yhatcoef_D0_all,yhatfd_D0_all, ...
    yhatfd_D1_all, yhatfd_D1_all_alt, ...
    evalarg_all, opts, opts_all] ...
    = PH_Betti0_functionFit_sparseDenseNetwork(filtrationMatrix, betti0Matrix, vb)

NTRIALS = size(filtrationMatrix,2);

yhatmat_D0 = cell(1,NTRIALS);
yhatmat_D1 = cell(1,NTRIALS);

warning off

% tuning parameters
% for numerical stability overcompensate with pseudo-interpolation just to 
% surely handle the tricky cases that have both
% extremely small and large gradients within the same curve. 
% obviously difficult to remove ringing in cases with a near step-change.
opts.NINT_h=200;                        
opts.NINT_v=50;                         
opts.domainLimExt = [-0.1, 0.1];        
opts.rangeLimExt = [-0.0001, 0.0001];
opts.NknotSkip = 8;
opts.norderF = 4;
opts.nlambda = 20;
opts.domainLim = [0, 1];
opts.NEVALPTS_indivFn = 15001;
opts_all.Nbasis_allFn = 4001;
opts_all.domainLimFinal = opts.domainLim;

for lp = 1:NTRIALS
    if vb,fprintf('\nFit barcode : %d \n', lp), end
    filtBuf = filtrationMatrix(:,lp);
    betti0Buf = betti0Matrix(:,lp);
    
    [yhatmat_D0{lp}, yhatmat_D1{lp}, narg_DOMAINEXTMASK, evalarg_DOMAINEXTMASK] = ...
        barcodeFunctionFit(filtBuf, betti0Buf, opts);
end

[yhatcoef_D0_all,evalarg_all,...
    yhatfd_D0_all, yhatfd_D1_all, yhatfd_D1_all_alt] = ...
    barcodeFunctionProject_all(yhatmat_D0, yhatmat_D1, narg_DOMAINEXTMASK, evalarg_DOMAINEXTMASK, opts_all);

warning on
end

function [yhatmat_D0_buf, yhatmat_D1_buf, narg_DOMAINEXTMASK, evalarg_DOMAINEXTMASK,...
    filtBuf_trunc, betti0Buf_reversed_trunc] = ...
    barcodeFunctionFit(filtBuf, betti0Buf, opts)

NINT_h = opts.NINT_h;
NINT_v = opts.NINT_v;
domainLimExt = opts.domainLimExt;
rangeLimExt = opts.rangeLimExt;
NknotSkip = opts.NknotSkip;
norderF = opts.norderF;
nlambda = opts.nlambda;
NEVALPTS_indivFn = opts.NEVALPTS_indivFn;

betti0Buf_reversed = max(betti0Buf) - betti0Buf;

% % assumes min(filtBuf) == 0, filtBuf(end)==1 for sparse or dense until here
endInd = find(filtBuf==1,1, 'first');
filtBuf_trunc = filtBuf(1:endInd); % truncate, or leave untouched if dense network
% must be wanting monotonic INCREASING fn for this: (not decreasing)
betti0Buf_reversed_trunc = [betti0Buf_reversed(1:(endInd-1)); betti0Buf_reversed(endInd-1)]; % truncate, or leave untouched if dense network

% special case for a bug, where filtBuf has duplicates (after sparse
% cleaning above)
filtBuf_trunc = filterDupes(filtBuf_trunc); % add \eps to duplicates
betti0Buf_reversed_trunc = filterDupes(betti0Buf_reversed_trunc); % add \eps to duplicates

filtBufRng = range(filtBuf_trunc);
betti0BufRng = range(betti0Buf_reversed_trunc);

filtRngExtLow = filtBuf_trunc(1)-abs(domainLimExt(1))*filtBufRng;
filtRngExtHigh = filtBuf_trunc(end)+abs(domainLimExt(end))*filtBufRng;
betti0RngExtLow = betti0Buf_reversed_trunc(1)-abs(rangeLimExt(1))*betti0BufRng;
betti0RngExtHigh = betti0Buf_reversed_trunc(end)+abs(rangeLimExt(end))*betti0BufRng;

fraw = [filtRngExtLow; filtBuf_trunc; filtRngExtHigh];
braw = [betti0RngExtLow; betti0Buf_reversed_trunc; betti0RngExtHigh];

FI_h = griddedInterpolant(fraw, braw, 'linear');
hfine_f = linspace(fraw(1),fraw(end),NINT_h)';
hfine_f = unique(sort(hfine_f));
vfine_f = FI_h(hfine_f);

FI_v = griddedInterpolant(braw, fraw, 'linear');
vfine_b = linspace(braw(1), braw(end), NINT_v)';
vfine_b = unique(sort(vfine_b));
hfine_b = FI_v(vfine_b);

[hfinesort, finesortidx] = sort([hfine_f; hfine_b ; fraw ]);
vfine_tosort = [vfine_f; vfine_b; braw];
vfinesort = vfine_tosort(finesortidx);
[hFineSortUnique, hfineuniqueidx] = unique(hfinesort);
vFineSortUnique = vfinesort(hfineuniqueidx);

filtToFit = hFineSortUnique;
betti0ToFit = vFineSortUnique;

knotsF = unique(filtToFit(1:NknotSkip:end));
if knotsF(end)~=filtRngExtHigh, knotsF(end+1) = filtRngExtHigh; end

nbasisF = length(knotsF)+norderF-2;

domainLimFineInit = [knotsF(1) knotsF(end)];
wbasisF = create_bspline_basis(domainLimFineInit, nbasisF, norderF, knotsF);

lambda= logspace(-10,-1,nlambda);
SSE = zeros(nlambda,1);
WfdF_buf = cell(nlambda,1);
betaF_buf = cell(nlambda,1);
yhatfdF_buf = cell(nlambda,1);
Fstr_buf = cell(nlambda,1);
argvals_buf = cell(nlambda,1);
y_buf = cell(nlambda,1);
y2cMap_buf = cell(nlambda,1);
for lp2= 1:nlambda
    WfdParF = fdPar(fd(ones(nbasisF,1),wbasisF), 0, lambda(lp2));
    [WfdF_buf{lp2}, betaF_buf{lp2}, yhatfdF_buf{lp2},...
        Fstr_buf{lp2}, argvals_buf{lp2}, y_buf{lp2}, y2cMap_buf{lp2}]...
        = smooth_monotone(filtToFit, betti0ToFit, WfdParF,...
        [], ones(numel(filtToFit),1), 0.0001, 50, [], 0);
    SSE(lp2) = Fstr_buf{lp2}.fn;
end
[~,minInd] = min(SSE);

WfdF = WfdF_buf{minInd};
betaF= betaF_buf{minInd};

rangeval = [filtToFit(1) filtToFit(end)];
fooprime = primes(NEVALPTS_indivFn);
narg_primeReduce     = fooprime(end);

evalarg = create_evalarg_withSpecificInteriorPts(rangeval(1), rangeval(2), narg_primeReduce, sort(opts.domainLim) );

hmat_D0     = eval_mon(evalarg, WfdF);
hmat_D1     = eval_mon(evalarg, WfdF,1);
DOMAINEXTMASK = evalarg>=min(opts.domainLim) & evalarg<= max(opts.domainLim);
evalarg_DOMAINEXTMASK = evalarg(DOMAINEXTMASK);
narg_DOMAINEXTMASK = numel(evalarg_DOMAINEXTMASK);
yhatmat_D0_buf = betaF(1) + betaF(2)*hmat_D0(DOMAINEXTMASK);
yhatmat_D1_buf = betaF(2)*hmat_D1(DOMAINEXTMASK);


end

function [yhatcoef_D0_all,evalarg_all,...
    yhatfd_D0_all, yhatfd_D1_all, yhatfd_D1_all_alt] = ...
        barcodeFunctionProject_all(yhatmat_D0, yhatmat_D1, narg_DOMAINEXTMASK, evalarg_DOMAINEXTMASK, opts_all)

Nbasis_allFn = opts_all.Nbasis_allFn;
domainLimFinal = opts_all.domainLimFinal;
evalarg_all = evalarg_DOMAINEXTMASK;
fooprime_all = primes(Nbasis_allFn);
nbasisF_all = fooprime_all(end);
norderF_all = 4;
wbasisF_all = create_bspline_basis(domainLimFinal, nbasisF_all, norderF_all);
yhatcoef_D0_all = project_basis(cell2mat(yhatmat_D0), evalarg_all, wbasisF_all, 1);
yhatcoef_D1_all = project_basis(cell2mat(yhatmat_D1), evalarg_all, wbasisF_all, 1);
yhatfd_D0_all = fd(yhatcoef_D0_all, wbasisF_all);
yhatfd_D1_all = fd(yhatcoef_D1_all, wbasisF_all);
% different approach to deal with ringing artefacts
yhatfd_D1_all_alt = deriv_fd(yhatfd_D0_all, 1);  

end

function evalarg = create_evalarg_withSpecificInteriorPts(lowerLim, upperLim, narg, interiorPts)
evalarg  = linspace(lowerLim, upperLim, narg)';
evalarg = sort([evalarg; interiorPts(:)]);
end

function data = filterDupes(data)

indexToDupes = @(dataIn, idxDataIn) find(not(ismember(1:numel(dataIn),idxDataIn)));
[~, iUniqueData, ~] = unique(data, 'first');
iDupes = indexToDupes(data, iUniqueData);
while ~isempty(iDupes) % case for multiple duplicates
    data(iDupes) = data(iDupes)+eps(data(iDupes));
    [~, iUniqueData, ~] = unique(data, 'first');
    iDupes = indexToDupes(data, iUniqueData);
end

end
