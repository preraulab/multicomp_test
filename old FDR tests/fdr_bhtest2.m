function [sigbins, p_adj, p_values] =  fdr_bhtest2(varargin)
%FDRtest2 Computes FDR regions of significance for a 2d matrix without independence assumption
%
%   Usage:
%   FDRtest2() RUNS DEMO
%   [sigbins, p_adj, p_values] =  FDRtest2(group1, group2, 'FDR', FDR)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required
%   FDR:  level for FDR acceptance (Default: 0.1)
%   use_mattest: logical flag for whether to use mattest (default: false)
%   mattest_options: Name-Value argument cell array for mattest function
%   calls (default: {'permute',false}) 
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
%      FDRtest2(); %RUNS DEMO
%
%   Copyright 2021 Michael J. Prerau, Ph.D.
%
%   Last modified 12/08/2021
%********************************************************************

%Call the examples for no input
if nargin ==0
    demo;
    return;
end

%%
p = inputParser;
addRequired(p,'group1',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addRequired(p,'group2',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addOptional(p,'FDR',0.1,@(x)validateattributes(x,{'numeric','1d'},{'nonempty','positive','<=',1}));
addOptional(p,'use_mattest',false,@islogical);
addOptional(p,'mattest_options',{'permute',false},@iscell);
addOptional(p,'paired',true,@islogical);
addOptional(p,'ploton',true,@islogical);

parse(p,varargin{:});

input_arguments = struct2cell(p.Results);
input_flags = fieldnames(p.Results);
eval(['[', sprintf('%s ', input_flags{:}), '] = deal(input_arguments{:});']);

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

%Compute linear perm test with global bounds
[linear_sigbins, linear_p_adj, linear_p_values] = fdr_bhtest(g1_redim, g2_redim,  FDR,...
    use_mattest, mattest_options, paired, false);

%Reshape output
sigbins = reshape(linear_sigbins, R,C);
p_adj = reshape(linear_p_adj, R,C);
p_values = reshape(linear_p_values, R,C);

%Plot results
if ploton
    figure
    hold all
  
    imagesc(nanmean(reshape(g1_redim, R, C, N1),3)-nanmean(reshape(g2_redim, R, C, N2),3));
    
    cx = climscale;
    caxis(max(abs(cx))*[-1 1]);
    colormap(redbluemap);
    colorbar
    
    %Plot significant regions as a contour
    if(any(sigbins(:)))
        [~, h_sigregions] = contour(sigbins,[1 1],'color','k','linewidth',1.5);
        legend(h_sigregions,'Regions of Significance');
    end
   
    axis tight;
end

function demo

%Define dataset
N1 = 20;
N2 = 20;
N = N1 + N2;

%Set peaks resolution
T = 50;

group1 = zeros(T,T,N1);
group2 = zeros(T,T,N2);

%Create functions
f1 =@(x)peaks(x)*rand*10 + randn(T)*5;
f2 =@(x)fliplr(peaks(x))*rand*10 + randn(T)*5;

%Generate data
for ii = 1:N
    if ii<=N1
        group1(:,:,ii) = f1(T);
    else
        group2(:,:,ii-N1) = f2(T);
    end
end

group1(:,30,1) = inf;

FDRtest2(group1,group2, 'paired', true);


