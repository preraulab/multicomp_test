function [sig_regions, acceptance_bounds, true_stat, ax] = gperm2(group1, group2, xvals, yvals, varargin)
%GLOBALPERMTEST Computes global acceptance bounds and regions of significance
%for a given statistic for two sets of multidimensional observations
%
%   Usage:
%   globalpermtest() RUNS DEMO
%   [sig_regions, acceptance_bounds] = gperm2(group1, group2, xvals, yvals, alpha_level, statfcn, iterations, ploton)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required
%   alpha_level:  global acceptance alpha (Default: 0.05)
%   statfcn: handle statistic to compute across trials (Default: @(x)mean(x,2) MUST COMPUTE ACROSS DIM 2)
%   ploton: (default: true) 
%
%   Output:
%   sig_regions: A cell array of contiguous significant dimensions
%   acceptance_bounds: the 2 x <dimensions> global acceptance bounds at a level defined by alpha_level
%
%   Example:
%      globalpermtest(); %RUNS DEMO
%
%   Copyright 2021 Michael J. Prerau, Ph.D.
%
%   Last modified 11/01/2021
%********************************************************************

%Call the examples for no input
if nargin==0
    
    
end


g1_redim = reshape(group1,size(group1,1)*size(group1,2),size(group1,3));
g2_redim = reshape(group2,size(group2,1)*size(group2,2),size(group2,3));

gperm(g1_redim, g2_redim);
