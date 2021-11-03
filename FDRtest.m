function [sigbins_all, sig_regions] = FDRtest(group1, group2, alpha_level, iterations, ploton)
%GLOBALPERMTEST Computes global acceptance bounds and regions of significance
%for a given statistic for two sets of multidimensional observations
%
%   Usage:
%   globalpermtest() RUNS DEMO
%   [sig_regions, acceptance_bounds] = globalpermtest(group1, group2, x_bins, alpha_level, statfcn, iterations, group1_name, group2_name, ploton)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required
%   alpha_level:  global acceptance alpha (Default: 0.05)x a
%   statfcn: handle statistic to compute across trials (Default: @(x)mean(x,2) MUST COMPUTE ACROSS DIM 2)
%   ploton: 0 no plot, 1 plot traces (default), 2 plot mean and 2*standard errors     % Previously boolean for plot output (Default: true)
%
%   Output:
%   sig_regions: A cell array of contiguous significant dimensions
%   acceptance_bounds: the 2 x <dimensions> global acceptance bounds at a level defined by alpha_level
%
%   Example:
%      globalpermtest(); %RUNS DEMO
%
%   Copyright 2021 Michael J. Prerau, Ph.D.
%
%   Last modified 11/01/2021
%********************************************************************

%Call the examples for no input
if nargin==0
    demo;
    return;
end


if nargin<3 || isempty(alpha_level)
    alpha_level = 0.1;
end

if nargin<4 || isempty(iterations)
    iterations = false;
end

if nargin <5 || isempty(ploton)
    ploton = true;
end

%Remove nan dimensions
p1 = group1(:,any(~isnan(group1)));
p2 = group2(:,any(~isnan(group2)));

if iterations
    pvalues = mattest(p1,p2,'permute',iterations);
else
    [~,pvalues]=ttest2(p1',p2');
end

padj = mafdr(pvalues,'BHFDR',true);

%Find the significant bins
sigbins_all=padj<alpha_level;
[cons_all,sig_regions]=consecutive(sigbins_all);

%Plot the results
if ploton
    figure('units','normalized','position',[0 0 1 1],'color','w');
    hold all;
    xvals = 1:length(padj);
    plot(xvals, padj);
    hline(alpha_level);
    
    %Plot significant regions
    yl=ylim;
    for i=1:length(cons_all)
        inds=sig_regions{i};
        h=fill(xvals([inds(1) inds(1) inds(end) inds(end)]),[yl(1) yl(2) yl(2) yl(1)],'g','edgecolor','none');
        uistack(h,'bottom');
    end
    axis tight;
end
end


function demo

%Define dataset
N1 = 200;
N2 = 303;
N = N1+N2;

%Set time resolution
T = 100;

%Initialize data and null matrices
g1 = zeros(T,N1);
g2 = zeros(T,N2);

%Create functions
f1 =@(x)normpdf(x)*(rand*5+3);
f2 =@(x)normpdf(x)*(rand*5+2);

x=linspace(-5,5,T)';

%Generate data
for ii = 1:N
    if ii<=N1
        g1(:,ii) = smooth(f1(x)+randn(size(x))*.125);
    else
        g2(:,ii-N1) = smooth(f2(x)+randn(size(x))*.125);
    end
end

FDRtest(g1, g2);
end