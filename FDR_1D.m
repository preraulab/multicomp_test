function [sigbins_all, p_adj, p_values] = FDR_1D(varargin)
% FDR_1D Perform statistical comparisons between 1D vectors using a False
% Discover Rate (FDR) approach
%
% Usage:
%   FDR_1D(group1, group2, <options>)
%
% Input:
%   - group1: Numeric 1D vector representing the first group of data.
%   - group2: Numeric 1D vector representing the second group of data.
%   - FDR: (Optional) False Discovery Rate (FDR) threshold for multiple testing correction. Default is 0.1.
%   - method: (Optional) 'dependent' will use the Benjamini & Yekutieli (2001) procedure, 
%     and 'independent' will use the Benjamini & Hochberg (1995) procedure that assumes data are independent or positively
%     dependent. Default is 'dependent'
%   - paired: (Optional) Boolean indicating whether the data in group1 and group2 are paired. Default is false.
%   - nonparam: (Optional) Boolean indicating whether to use nonparametric test. 
%               For nonparametric tests, a ranksum test is use for unpaired
%               an a signrank test is used for paired. Paired an unpaired
%               t-tests are used otherwise
%               Default is true.
%   - ploton: (Optional) Boolean indicating whether to plot the results. Default is true.
%
% Output:
%   - sigbins: Vector indicating significant regions between group1 and group2.
%   - p_adj: Vector of adjusted p-values after multiple testing correction.
%   - p_values: Vector of raw p-values.
%
% Example:
%   %Run FDR_1D() for demo data
%
%     %Define dataset
%     N1 = 200;
%     N2 = 200;
%     N = N1+N2;
%     
%     %Set time resolution
%     T = 100;
%     
%     %Initialize data and null matrices
%     g1 = zeros(T,N1);
%     g2 = zeros(T,N2);
%     
%     %Create functions
%     f1 =@(x)normpdf(x)*(rand*5+3);
%     f2 =@(x)normpdf(x)*(rand*5+2);
%     
%     x=linspace(-5,5,T)';
%     
%     %Generate data
%     for ii = 1:N
%         if ii<=N1
%             data1 = smooth(f1(x)+randn(size(x))*.25)+1;
%             data1(end-10:end) = nan;
%             g1(:,ii) = data1;
%         else
%             data2 = smooth(f2(x)+randn(size(x))*.25)+1;
%             data2(end-3:end) = nan;
%             g2(:,ii-N1) = data2;
%         end
%     end
% 
%     %Run nonparametric test
%     FDR_1D(group1,group2,'FDR', .1,'paired', true, 'method', 'dependent','nonparam',false);
%
%     %Run parametric test
%     FDR_1D(group1,group2,'FDR', .1,'paired', false,'method', 'independent','nonparam',true);
%
% See also:
%   FDR_2D, fdr_bh
%
%   Copyright 2023 Michael J. Prerau Laboratory. - http://www.sleepEEG.org

% DEMO
if nargin == 0
    FDR = .1;
    paired = false;
    demo_func(FDR,true,true,'dependent');
    demo_func(FDR,false,false,'independent');
    return;
end

p = inputParser;
addRequired(p,'group1',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addRequired(p,'group2',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addOptional(p,'FDR',0.1,@(x)validateattributes(x,{'numeric','1d'},{'nonempty','positive','<=',1}));
addOptional(p,'method', 'dependent', @(x)ismember(x,{'dependent','independent'}));
addOptional(p,'paired',false,@islogical);
addOptional(p,'nonparam',true,@islogical);
addOptional(p,'ploton',true,@islogical);

parse(p,varargin{:});

input_arguments = struct2cell(p.Results);
input_flags = fieldnames(p.Results);
eval(['[', sprintf('%s ', input_flags{:}), '] = deal(input_arguments{:});']);

%Change infs to nans
group1(isinf(group1)) = nan;
group2(isinf(group2)) = nan;
if nonparam
    if paired
        assert(size(group1,2) == size(group2,2),'There must be the same number of observations for a paired test.')
        testfun = @signrank;
    else
        testfun = @ranksum;
    end
else   
    if paired 
        [~,p_values]=ttest(group1',group2');
    else
        [~,p_values]=ttest2(group1',group2');
    end
    
end

if nonparam
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
   
end

switch method
    case 'dependent'
        method_string = 'dep';
    case 'independent'
        method_string = 'pdep';
    otherwise
        error('Invalid method string')
end

[~,~,~,p_adj] = fdr_bh(p_values, FDR, method_string);

%Find the significant bins
sigbins_all = p_adj<FDR;

%Plot the results
if ploton
    figure
    ax = figdesign(2,1,'type','usletter','orient','landscape');
    axes(ax(1))
    hold all;
    plot(group1,'color',[1 0 0 .1])
    plot(group2,'color',[0 0 1 .1])
    legend('Group 1','Group 2')

    axes(ax(2))
    [cons_all,sig_regions]=consecutive_runs(sigbins_all);
    [cons_nan,nan_regions]=consecutive_runs(isnan(p_values));
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

    h_nanregions = [];
    for ii = 1:length(cons_nan)
        inds = nan_regions{ii};

        if ~isempty(inds)
            h_nanregions(ii) = fill([xvals(inds(1))-.5 xvals(inds(1))-.5 xvals(inds(end))+.5 xvals(inds(end))+.5],[yl(1) yl(2) yl(2) yl(1)],'r','edgecolor','none');
        end
    end

    if ~isempty(h_nanregions)
        uistack(h_nanregions,'bottom');
    end

    %Plot adjusted pvalues and threshold line
    h_pasj = plot(xvals, p_adj,'linewidth',2,'color','b');
    axis tight;
    h_threshold = hline(FDR,'color','k','linestyle','--');

    if isempty(h_sigregions)
        legend([h_pasj, h_threshold],{'Adjusted p-values','Threshold'});
    else
        if ~isempty(h_nanregions)
            legend([h_pasj, h_threshold, h_sigregions(1), h_nanregions(1)],{'Adjusted p-values','Threshold','Significant Regions','No Data Regions'});
        else
            legend([h_pasj, h_threshold, h_sigregions(1)],{'Adjusted p-values','Threshold','Significant Regions'});
        end
    end

    xlabel('Bin Number');
    ylabel('Adjusted p-value');

    if paired
        pstring = 'Paired';
    else
        pstring = 'Unpaired';
    end

    if nonparam
        npstring = 'Nonparametric';
    else
        npstring = 'Parametric';
    end

    mstring = [upper(method(1)) method(2:end)];

    suptitle([mstring ' ' pstring ' ' npstring ' Test with FDR of ' num2str(FDR)])
end
end


function demo_func(FDR,paired, nonparam, method)
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
        data1 = smooth(f1(x)+randn(size(x))*.25)+1;
        data1(end-10:end) = nan;
        g1(:,ii) = data1;
    else
        data2 = smooth(f2(x)+randn(size(x))*.25)+1;
        data2(end-3:end) = nan;
        g2(:,ii-N1) = data2;
    end
end

FDR_1D(g1,g2, 'FDR',FDR,'method',method,'paired', paired,'nonparam',nonparam);
end