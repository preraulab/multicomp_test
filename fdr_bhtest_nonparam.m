function [sigbins_all, p_adj, p_values] = fdr_bhtest_nonparam(varargin)
%fdr_bhtest_nonparam Computes FDR regions of significance without independence assumption
%
%   [sigbins, p_adj, p_values] =  fdr_bhtest_nonparam(group1, group2, FDR)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required
%   FDR:  level for FDR acceptance (Default: 0.1)
%   paired: logical flag for doing paired sample ranksign or ranksum as main hypothesis
%       testing statistics (default: false)
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
%   Last modified 12/08/2021
%********************************************************************

%Call the examples for no input
if nargin==0
    demo;
    return;
end

p = inputParser;
addRequired(p,'group1',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addRequired(p,'group2',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addOptional(p,'FDR',0.1,@(x)validateattributes(x,{'numeric','1d'},{'nonempty','positive','<=',1}));
addOptional(p,'paired',false,@islogical);
addOptional(p,'ploton',true,@islogical);

parse(p,varargin{:});

input_arguments = struct2cell(p.Results);
input_flags = fieldnames(p.Results);
eval(['[', sprintf('%s ', input_flags{:}), '] = deal(input_arguments{:});']);

%Change infs to nans
group1(isinf(group1)) = nan;
group2(isinf(group2)) = nan;

if paired
    assert(size(group1,2) == size(group2,2),'There must be the same number of observations for a paired test.')
    testfun = @signrank;
else
    testfun = @ranksum;
end

p_values = zeros(1,size(group1,1));
for ii = 1:length(p_values)
    g1 = group1(ii,:);
    g2 = group2(ii,:);

    if all(isnan(g1)) && all(isnan(g2))
        p_values(ii) = nan;
    else
        if any(isnan(g1)) || all(isnan(g2))
            g1(isnan(g1)) = 0;
            g2(isnan(g2)) = 0;
        end

        p_values(ii)=testfun(g1,g2);
    end
end


[~,~,~,p_adj] = fdr_bh(p_values, FDR, 'dep');

%Find the significant bins
sigbins_all = p_adj<FDR;

%Plot the results
if ploton
    figure
    [cons_all,sig_regions]=consecutive_runs(sigbins_all);
    [cons_nan,nan_regions]=consecutive_runs(isnan(p_values));

    %     figure('units','normalized','color','w');
    hold all;
    xvals = 1:length(p_adj);

    %Plot significant regions
    yl=ylim;

    %Plot regions of significance
    h_sigregions = [];

    for ii = 1:length(cons_all)
        inds = sig_regions{ii};

        if ~isempty(inds)
            h_sigregions(ii) = fill([xvals(inds(1))-.5 xvals(inds(1))-.5 xvals(inds(end))+.5 xvals(inds(end))+.5],[yl(1) yl(2) yl(2) yl(1)],'g','edgecolor','none');
        end
    end
    uistack(h_sigregions,'bottom');

    for ii = 1:length(cons_nan)
        inds = nan_regions{ii};

        if ~isempty(inds)
            h_nanregions(ii) = fill([xvals(inds(1))-.5 xvals(inds(1))-.5 xvals(inds(end))+.5 xvals(inds(end))+.5],[yl(1) yl(2) yl(2) yl(1)],'r','edgecolor','none');

        end
    end
    uistack(h_nanregions,'bottom');

    %Plot adjusted pvalues and threshold line
    h_pasj = plot(xvals, p_adj,'linewidth',2,'color','b');
    axis tight;
    h_threshold = hline(FDR,'color','k','linestyle','--');

    if isempty(h_sigregions)
        legend([h_pasj, h_threshold],{'Adjusted p-values','Threshold'});
    else
        legend([h_pasj, h_threshold, h_sigregions(1)],{'Adjusted p-values','Threshold','Significant Regions','No Data Regions'});
    end

    xlabel('Bin Number');
    ylabel('Adjusted p-value');
end
end


function demo
close all;
%Define dataset
N1 = 200;
N2 = 200;
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
        data1 = smooth(f1(x)+randn(size(x))*.125)+1;
        data1(end-10:end) = nan;
        g1(:,ii) = data1;
    else
        data2 = smooth(f2(x)+randn(size(x))*.125)+1;
        data2(end-3:end) = nan;
        g2(:,ii-N1) = data2;
    end
end

fdr_bhtest_nonparam(g1, g2, 'paired', false);
end