function [sigbins, acceptance_bounds, true_stat] = gpermtest2(varargin)
%GPERMTEST2 Computes global acceptance bounds and regions of significance
%for a given statistic for two sets of multidimensional observations
%
%   Usage:
%   gpermtest2() RUNS DEMO
%   [sigbins, acceptance_bounds, true_stat] = gpermtest2(group1, group2, alpha_level)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required
%   alpha_level:  global acceptance alpha (Default: 0.05)
%   statfcn: handle statistic to compute across trials (Default: @(x)nanmean(x,2) MUST COMPUTE ACROSS DIM 2)
%   iterations: numeric, number of iterations to generate null distribution
%   ploton: boolean for plot output traces and significance regions (default: true)
%
%   Output:
%   sigbins: A matrix of signifiance bins
%   sig_regions: A cell array of contiguous significant dimensions
%   acceptance_bounds: the 2 x <dimensions> global acceptance bounds at a level defined by alpha_level
%   true_stat: difference in mean activity between the two groups
%
%   Example:
%      gpermtest2(); %RUNS DEMO
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

%Parse inputs
p = inputParser;
addRequired(p,'group1',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addRequired(p,'group2',@(x)validateattributes(x,{'numeric','2d'},{'nonempty'}));
addOptional(p,'alpha_level',0.05,@(x)validateattributes(x,{'numeric','1d'},{'nonempty','positive','<=',1}));
addOptional(p,'statfcn', @(x)nanmean(x,2),@(x)isa(x,'function_handle'));
addOptional(p,'iterations',10000,@(x)validateattributes(x,{'numeric','1d'},{'nonempty','positive'}));
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

%Compute linear perm test with global bounds
[linear_sigbins, ~, linear_bounds, linear_true_stat] = gpermtest(g1_redim, g2_redim, alpha_level, statfcn, iterations, false);

%Reshape output
sigbins = reshape(linear_sigbins, R,C);
acceptance_bounds = reshape(linear_bounds, R,C);
true_stat = reshape(linear_true_stat, R,C);

if ploton
    figure
    
    hold all
    imagesc(nanmean(reshape(g1_redim, R, C, N1),3)-nanmean(reshape(g2_redim, R, C, N2),3));
    
    cx = climscale;
    caxis(max(abs(cx))*[-1 1]);
    colormap(redbluemap);
    colorbar
    
    if any(sigbins,'all')
        [~,h_sigregions] = contour(sigbins,1,'color','k', 'LineWidth', 1.5);
        legend(h_sigregions,'Significant Regions');
    end
end
end 


function demo

%Define dataset
N1 = 20;
N2 = 33;
N = N1 + N2;

%Set peaks resolution
T = 50;

%Initialize data and null matrices
group1 = zeros(T,T,N1);
group2 = zeros(T,T,N2);

%Create functions
f1 =@(x)peaks(x)*rand*10 + randn(T)*.5;
f2 =@(x)fliplr(peaks(x))*rand*10 + randn(T)*.5;

%Generate data
for ii = 1:N
    if ii<=N1
        group1(:,:,ii) = f1(T);
    else
        group2(:,:,ii-N1) = f2(T);
    end
end

gpermtest2(group1, group2);
end
