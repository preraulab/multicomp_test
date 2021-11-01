function [sigbins_all, sig_regions, acceptance_bounds, true_stat, ax] = gperm(varargin)
%GLOBALPERMTEST Computes global acceptance bounds and regions of significance
%for a given statistic for two sets of multidimensional observations
%
%   Usage:
%   globalpermtest() RUNS DEMO
%   [sig_regions, acceptance_bounds] = globalpermtest(group1, group2, x_bins, alpha_level, statfcn, iterations, group1_name, group2_name, ploton)
%
%   Input:
%   group1, group2: in form <dimensions> x <trials> -- required
%   alpha_level:  global acceptance alpha (Default: 0.05)
%   statfcn: handle statistic to compute across trials (Default: @(x)mean(x,2) MUST COMPUTE ACROSS DIM 2)
%   ploton: 0 no plot, 1 plot traces (default), 2 plot mean and 2*standard errors     % Previously boolean for plot output (Default: true)
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
    demo;
    return;
end


defaults={'','',[],.05,@(x)mean(x,2),10000,'Group 1','Group 2',1};
%Handle defaults
defaults(setdiff(1:length(varargin), find(cellfun(@isempty,varargin)))) = varargin(~cellfun(@isempty,(varargin)));


p1=defaults{1};
p2=defaults{2};
xvals=defaults{3};

alpha_level=defaults{4};
statfcn=defaults{5};
iterations=defaults{6};
group1_name=defaults{7};
group2_name=defaults{8};
ploton=defaults{9};

locations=size(p1,1);
if isempty(xvals)
    xvals=1:locations;
end

%The difference in mean activity between the two scenarios
true_stat=statfcn(p1)-statfcn(p2);

%Combine both groups
all=[p1 p2];

Np1=size(p1,2);

%Generate null distribution
null_mat=zeros(locations,iterations);

disp(['Computing test with ' num2str(iterations) ' iterations at alpha level ' num2str(alpha_level) '...']);

parfor i=1:iterations
    %Get a random permutation of the labels
    inds=randperm(size(all,2));
    
    %Indices for new pseudo-groups
    inds1=inds(1:Np1);
    inds2=inds((Np1 + 1):end);
    
    %Compute the stat and take absolute value
    null_mat(:,i)=abs((statfcn(all(:,inds1))-statfcn(all(:,inds2))));
end

%Sort the trials for ease of computing the pointwise N-tile
sorted_null=sort(null_mat,2);

%Set initial bounds at highest %ile
cutoff_ind = iterations;

% Number of individual trials falling outside the global bounds
numout=0;

%Compute the citerion for the number out
out_crit = ceil(alpha_level*iterations*2);

gbounds=sorted_null(:,end);

%Shrink the bounds until the number of outside trials is %alpha
if out_crit<iterations && out_crit>0
    while numout<out_crit && cutoff_ind>0
        %Calculate the global bounds
        gbounds=sorted_null(:,cutoff_ind);
        
        %Check the number out of bounds
        numout=sum(any(null_mat>=gbounds));
        
        %Shrink the bounds
        cutoff_ind = cutoff_ind - 1;
    end
elseif out_crit>=iterations
    warning('Global bounds includes all iterations. Increase iterations or reduce p-value');
    gbounds = zeros(size(sorted_null(:,1)));
elseif out_crit == 0
    warning('Global bounds is below precision for number of iterations. Increase iterations or increase p-value');
    gbounds = sorted_null(:,end);
end

%Create symmetrical bounds
hi = gbounds;
lo = -gbounds;

%Find the significant bins
sigbins_all=(true_stat>=hi | true_stat<=lo);
acceptance_bounds=[hi,lo];
[cons_all,sig_regions]=consecutive(sigbins_all);

%Plot the results
if ploton
    figure1=figure('units','normalized','position',[0 0 1 1],'color','w');
    
    ax(1)  = axes('Parent',figure1,...
        'Position',[0.08 0.376666666666667 0.87 0.246666666666667*2]);
    ax(2)  = axes('Parent',figure1,'Position',[0.08 0.05 0.87 0.246666666666667]);
    
    linkaxes(ax,'x');
    
    axes(ax(1))
    hold on;
    if ploton==2
        m1 = mean(p1,2,'omitnan');
        se1 = std(p1,0,2,'omitnan')/sqrt(size(p1,2));
        lo1 = m1-2*squeeze(se1);
        hi1 = m1+2*squeeze(se1);
        shadebounds(xvals, m1, hi1, lo1, 'r', [1 .7 .7], 'none');
        m2 = mean(p2,2,'omitnan');
        se2 = std(p2,0,2,'omitnan')/sqrt(size(p2,2));
        lo2 = m2-2*squeeze(se2);
        hi2 = m2+2*squeeze(se2);
        shadebounds(xvals, m2, hi2, lo2, 'b', [.7 .7 1], 'none');
    else
        plot(xvals,p1,'r');
        plot(xvals,p2,'b');
    end
    title([group1_name ' vs. ' group2_name],'fontsize',15);
    axis tight;
    xlim(xvals([1 end]));
    
    axes(ax(2))
    hold on;
    
    %Plot global acceptance bounds
    plot(xvals,acceptance_bounds,'r','linewidth',2);
    %Plot data mean
    plot(xvals,true_stat,'b','linewidth',2);
    axis tight
    
    %Plot significant regions
    yl=ylim;
    for i=1:length(cons_all)
        inds=sig_regions{i};
        h=fill(xvals([inds(1) inds(1) inds(end) inds(end)]),[yl(1) yl(2) yl(2) yl(1)],'g','edgecolor','none');
        uistack(h,'bottom');
    end
    axis tight;
    
    ylabel('Difference Statistic','fontsize',12);
    title('Global Acceptance Test','fontsize',15);
    
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
        g1(:,ii) = smooth(f1(x)+randn(size(x))*.25);
    else
        g2(:,ii-N1) = smooth(f2(x)+randn(size(x))*.25);
    end
end

gperm(g1, g2)
