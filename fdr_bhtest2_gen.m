function [sigbins, p_adj, p_values] =  fdr_bhtest2_gen(varargin)
%FDRtest2 Computes FDR regions of significance for a 2d matrix without independence assumption
%
%   Usage:
%   FDRtest2() RUNS DEMO
%   [sigbins, p_adj, p_values] =  FDRtest2(group1, group2, 'FDR', FDR)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required for non
%   demo
%   FDR:  level for FDR acceptance (Default: 0.1)
%   paired: logical flag for doing paired sample ttest or two sample ttest
%   as main hypothesis testing statistics (default: true)
%   nonparam: logical flag for doing nonparametric test (default:true)
%   ploton: (default: true)
%   demo_on: logical flag for running demo (default:false)
%
%   Output:
%   sigbins: A matrix of signifiance bins
%   p_adj: Adjusted p-values
%   p_values: raw p-values
%
%   Example:
%      FDRtest2(); %RUNS DEMO
%
%   Copyright 2021 Michael J. Prerau, Ph.D.
%
%   Last modified 12/08/2021
%********************************************************************

%%
p = inputParser;
addOptional(p,'group1',[],@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addOptional(p,'group2',[],@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addOptional(p,'FDR',0.1,@(x)validateattributes(x,{'numeric','1d'},{'nonempty','positive','<=',1}));
addOptional(p,'paired',false,@islogical);
addOptional(p,'nonparam',true,@islogical);
addOptional(p,'ploton',true,@islogical);
addOptional(p,'demo_on',false,@islogical);

parse(p,varargin{:});

input_arguments = struct2cell(p.Results);
input_flags = fieldnames(p.Results);
eval(['[', sprintf('%s ', input_flags{:}), '] = deal(input_arguments{:});']);

% DEMO
if demo_on
    demo_func(FDR,paired,nonparam);
    return;
else
    assert(~isempty(group1) & ~isempty(group2),'group1 and group 2 required for non demo')
end

%Get group sizes
[R,C,N1] = size(group1);
[R2,C2,N2] = size(group2);

%Check that they are the same size and are valid
assert(R == R2 && C == C2,'Group dims must be the same');
assert(any(isfinite(group1),'all') && any(isfinite(group2),'all'), 'Groups must have valid numeric data')

%Reshape to linear
g1_redim = reshape(group1,size(group1,1)*size(group1,2),size(group1,3));
g2_redim = reshape(group2,size(group2,1)*size(group2,2),size(group2,3));

%Set inf values to nan
g1_redim(isinf(g1_redim)) = nan;
g2_redim(isinf(g2_redim)) = nan;

%Compute test
[linear_sigbins, linear_p_adj, linear_p_values] = fdr_bhtest_gen('group1',g1_redim, 'group2',g2_redim,  ...
    'FDR',FDR,'paired',paired,'nonparam',nonparam,'ploton',ploton);

%Reshape output
sigbins = reshape(linear_sigbins, R,C);
p_adj = reshape(linear_p_adj, R,C);
p_values = reshape(linear_p_values, R,C);

%Plot results
if ploton
    g1_mean = mean(reshape(g1_redim, R, C, N1),3,'omitnan');
    g2_mean = mean(reshape(g2_redim, R, C, N2),3,'omitnan');

    figure
    ax = figdesign(2,2,'type','usletter','orient','landscape');
    axes(ax(1))
    imagesc(g1_mean)
    cx = climscale;
    colormap(ax(1),gouldian);
    title('Group 1')

    axes(ax(2))
    imagesc(g2_mean)
    caxis(cx);
    colormap(ax(2),gouldian);
    title('Group 2')

    axes(ax(3))
    hold all

    imagesc(g1_mean - g2_mean);
    [r,c] = find(isnan(p_adj));
    if ~isempty(r)
        plot(c,r,'x')
    end

    cx = climscale;
    caxis(max(abs(cx))*[-1 1]);
    colormap(ax(3),redblue_equalized);
    colorbar

    %Plot significant regions as a contour
    if(any(sigbins(:)))
        [~,hc] = contour(sigbins,[1 1],'color','k','linewidth',1.5);
        legend(hc,'Regions of Significance')
    end

    axis tight;
    title('Group 1 - Group 2')

    axes(ax(4))
    imagesc(p_adj)
    hold on
    if ~isempty(r)
        plot(c,r,'x')
    end

    caxis([0 FDR*2])
    axis xy
    colormap(ax(4),'parula');
    title('Adjusted p-values')
    colorbar
end

function demo_func(FDR,paired,nonparam)

%Define dataset
N1 = 20;
N2 = 30;
N = N1 + N2;

%Set peaks resolution
T = 30;

group1 = zeros(T,T,N1);
group2 = zeros(T,T,N2);

%Create functions
f1 = @(x)peaks(x)*rand*10 + randn(T)*5+50;
f2 = @(x)fliplr(peaks(x))*rand*10 + randn(T)*5+50;

%Generate data
for ii = 1:N
    if ii<=N1
        data = f1(T);
        data(:,end-6:end) = nan;
        group1(:,:,ii) = data;
    else
        data = f2(T);
        data(:,end-2:end) = nan;
        group2(:,:,ii-N1) = data;
    end
end

fdr_bhtest2_gen('group1',group1,'group2',group2,'FDR',FDR,'paired', paired,'nonparam',nonparam);


