function [sigbins, p_adj, p_values] =  FDRtest2(group1, group2, alpha_level, iterations, ploton)
%FDRtest2 Computes FDR regions of significance for a 2d matrix
%
%   Usage:
%   FDRtest2() RUNS DEMO
%   [sigbins, p_adj, p_values] =  FDRtest2(group1, group2, alpha_level, iterations, ploton)
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
%      FDRtest2(); %RUNS DEMO
%
%   Copyright 2021 Michael J. Prerau, Ph.D.
%
%   Last modified 11/03/2021
%********************************************************************

%Call the examples for no input
if nargin ==0
    demo;
    return;
end

%%
%Get group sizes
[R,C,N1] = size(group1);
[R2,C2,N2] = size(group2);

%Check that they are the same size
assert(R == R2 && C == C2,'Group dims must be the same');

%Set up default parameters
if nargin<3 || isempty(alpha_level)
    alpha_level = 0.1;
end

if nargin<4 || isempty(iterations)
    iterations = false;
end

if nargin <5
    ploton = true;
end

%Reshape to linear
g1_redim = reshape(group1,size(group1,1)*size(group1,2),size(group1,3));
g2_redim = reshape(group2,size(group2,1)*size(group2,2),size(group2,3));

%Compute linear perm test with global bounds
[linear_sigbins, linear_p_adj, linear_p_values] = FDRtest(g1_redim, g2_redim,  alpha_level,  iterations, false);

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
        [~, hc] = contour(sigbins,[1 1],'color','k','linewidth',1.5);
        legend(hc,'Regions of Significance');
    end
   
    axis tight;
end

function demo

%Define dataset
N1 = 20;
N2 = 33;
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

FDRtest2(group1,group2);


