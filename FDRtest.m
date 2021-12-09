function [sigbins_all, p_adj, p_values] = FDRtest(varargin)
%FDRTEST Computes FDR regions of significance
%
%   Usage:
%   FDRtest2() RUNS DEMO
%   [sigbins, p_adj, p_values] =  FDRtest(group1, group2, FDR)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required
%   FDR:  level for FDR acceptance (Default: 0.1)
%   use_mattest: logical flag for whether to use mattest (default: false)
%   mattest_options: Name-Value argument cell array for mattest function
%   calls (default: {'permute',false}) 
%   mafdr_options: Name-Value argument cell array for mafdr function calls
%   (default: {'BHFDR',true})
%   paired: logical flag for doing paired sample ttest or two sample ttest
%   as main hypothesis testing statistics (default: true)
%   ploton: (default: true)
%
%   Output:
%   sigbins: A matrix of signifiance bins
%   p_adj: Adjusted p-values
%   p_values: raw p-values
%
%   Example:
%       FDRtest2(); %RUNS DEMO
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
addOptional(p,'use_mattest',false,@islogical);
addOptional(p,'mattest_options',{'permute',false},@iscell);
addOptional(p,'mafdr_options',{'BHFDR',true},@iscell);
addOptional(p,'paired',true,@islogical);
addOptional(p,'ploton',true,@islogical);

parse(p,varargin{:});

input_arguments = struct2cell(p.Results);
input_flags = fieldnames(p.Results);
eval(['[', sprintf('%s ', input_flags{:}), '] = deal(input_arguments{:});']);

%Change infs to nans
group1(isinf(group1)) = nan;
group2(isinf(group2)) = nan;

if use_mattest
    assert(~paired, 'Cannot run paired ttest using mattest.')
    p_values = mattest(p1,p2,mattest_options{:});
else
    if paired 
        [~,p_values]=ttest(group1',group2');
    else
        [~,p_values]=ttest2(group1',group2');
    end
end

p_adj = mafdr(p_values,mafdr_options{:});

%Find the significant bins
sigbins_all = p_adj<FDR;

%Plot the results
if ploton
    [cons_all,sig_regions]=consecutive(sigbins_all); % this currently does not plot significant single time points (i.e., length(sig_regions{ii})<2)
    
    figure('units','normalized','color','w');
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
    h_threshold = hline(FDR,'color','k','linestyle','--');
    
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
        g1(:,ii) = smooth(f1(x)+randn(size(x))*.125);
    else
        g2(:,ii-N1) = smooth(f2(x)+randn(size(x))*.125);
    end
end

FDRtest(g1, g2, 'paired', true);
end