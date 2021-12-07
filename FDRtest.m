function [sigbins_all, p_adj, p_values] = FDRtest(group1, group2, alpha_level, iterations, ploton)
%FDRTEST Computes FDR regions of significance
%
%   [sigbins, p_adj, p_values] =  FDRtest(group1, group2, alpha_level, iterations, ploton)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required
%   alpha_level:  level for FDR acceptance (Default: 0.1)
%   iterations: if numeric, perform a permulation t-test with the given
%   number of iterations (default: false = single t-test)
%   ploton: (default: true)
%
%   Output:
%   sigbins: A matrix of signifiance bins
%   p_adj: Adjusted p-values
%   p_values: raw p-values
%
%   Example:
%      FDRtest(); %RUNS DEMO
%
%   Copyright 2021 Michael J. Prerau, Ph.D.
%
%   Last modified 11/03/2021
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
    p_values = mattest(p1,p2,'permute',iterations);
else
    [~,p_values]=ttest2(p1',p2');
end

p_adj = mafdr(p_values,'BHFDR',true);

%Find the significant bins
sigbins_all = p_adj<alpha_level;


%Plot the results
if ploton
    [cons_all,sig_regions]=consecutive(sigbins_all);
    
    figure('units','normalized','position',[0 0 1 1],'color','w');
    hold all;
    xvals = 1:length(p_adj);
    
    %Plot significant regions
    yl=ylim;
    h_r = [];
    
    for ii = 1:length(cons_all)
        inds = sig_regions{ii};
        if ~isempty(inds)
            h_r(ii) = fill(xvals([inds(1) inds(1) inds(end) inds(end)]),[yl(1) yl(2) yl(2) yl(1)],'g','edgecolor','none');
            uistack(h_r,'bottom');
        end
    end
    
    h_pa = plot(xvals, p_adj,'linewidth',2,'color','b');
    h_t = hline(alpha_level,'color','k','linestyle','--');
    
    if isempty(h_r)
        legend([h_pa, h_t],{'Adjusted p-values','Threshold'});
    else
        legend([h_pa, h_t, h_r(1)],{'Adjusted p-values','Threshold','Significant Regions'});
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