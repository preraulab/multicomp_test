function [sigbins_all, p_adj, p_values] = FDRtest(varargin)
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

p = inputParser;
addRequired(p,'group1',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addRequired(p,'group2',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addOptional(p,'alpha_level',0.1,@(x)validateattributes(x,{'numeric','1d'},{'nonempty','positive','<=',1}));
addOptional(p,'iterations',false,@islogical);
addOptional(p,'ploton',true,@islogical);

parse(p,varargin{:});

input_arguments = struct2cell(p.Results);
input_flags = fieldnames(p.Results);
eval(['[', sprintf('%s ', input_flags{:}), '] = deal(input_arguments{:});']);


%Change infs to nans
group1(isinf(group1)) = nan;
group2(isinf(group2)) = nan;

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
        
    %Plot regions of significance
    h_sigregions = [];
    
    for ii = 1:length(cons_all)
        inds = sig_regions{ii};
        
        if ~isempty(inds)
            h_sigregions(ii) = fill(xvals([inds(1) inds(1) inds(end) inds(end)]),[yl(1) yl(2) yl(2) yl(1)],'g','edgecolor','none');
            uistack(h_sigregions,'bottom');
        end
    end
    
    %Plot adjusted pvalues and threshold line
    h_pasj = plot(xvals, p_adj,'linewidth',2,'color','b');
    h_threshold = hline(alpha_level,'color','k','linestyle','--');
    
    if isempty(h_sigregions)
        legend([h_pasj, h_threshold],{'Adjusted p-values','Threshold'});
    else
        legend([h_pasj, h_threshold, h_sigregions(1)],{'Adjusted p-values','Threshold','Significant Regions'});
    end
    
    xlabel('Bin Number');
    ylabel('Adjusted p-value');
    
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