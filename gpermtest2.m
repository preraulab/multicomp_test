function [sigbins, acceptance_bounds, true_stat] = gpermtest2(group1, group2,  alpha_level, statfcn, iterations, ploton)
%gpermtest2 Computes global acceptance bounds and regions of significance
%for a given statistic for two sets of multidimensional observations
%
%   Usage:
%   gpermtest2() RUNS DEMO
%   [sig_regions, acceptance_bounds] = gpermtest2(group1, group2, xvals, yvals, alpha_level, statfcn, iterations, ploton)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required
%   alpha_level:  global acceptance alpha (Default: 0.05)
%   statfcn: handle statistic to compute across trials (Default: @(x)nanmean(x,2) MUST COMPUTE ACROSS DIM 2)
%   ploton: (default: true)
%
%   Output:
%   sig_regions: A cell array of contiguous significant dimensions
%   acceptance_bounds: the 2 x <dimensions> global acceptance bounds at a level defined by alpha_level
%
%   Example:
%      gpermtest2(); %RUNS DEMO
%
%   Copyright 2021 Michael J. Prerau, Ph.D.
%
%   Last modified 11/01/2021
%********************************************************************

%Call the examples for no input
if nargin == 0
    
    %Define dataset
    N1 = 20;
    N2 = 33;
    N = N1 + N2;
    
    %Set peaks resolution
    T = 50;
    
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
end

%%
%Get group sizes
[R,C,N1] = size(group1);
[R2,C2,N2] = size(group2);

%Check that they are the same size
assert(R == R2 && C == C2,'Group dims must be the same');

if nargin<3 || isempty(alpha_level)
    alpha_level = 0.05;
end

if nargin<4 || isempty(statfcn)
    statfcn = @(x)nanmean(x,2);
end

if nargin<5 || isempty(iterations)
    iterations = 100000;
end

if nargin <6
    ploton = true;
end

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



