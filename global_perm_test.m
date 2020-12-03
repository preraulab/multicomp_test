function [sig_regions, acceptance_bounds, data_stat, ax] = global_perm_test(varargin)
%GLOBALPERMTEST Computes global acceptance bounds and regions of significance
%for a given statistic for two sets of multidimensional observations
%
%   Usage:
%   globalpermtest() RUNS DEMO
%   [sig_regions, acceptance_bounds] = globalpermtest(group1, group2, x-values, alpha_level, statfcn, iterations, ploton)
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
%   Copyright 2015 Michael J. Prerau, Ph.D.
%
%   Last modified 05/21/2018
%********************************************************************

%Call the examples for no input
if nargin==0
    global_perm_test('var');
    title('Global Acceptance For Difference in Variance Across Conditions','fontsize',15);
    
    global_perm_test('mean');
    title('Global Acceptance For Difference in Mean Across Conditions','fontsize',15);
    
    % Create annotations
    annotation(gcf,'textbox',...
        [0.748210023866354 0.0764946921443736 0.112171837708825 0.0307855626326964],...
        'String',{'Significant Regions'},...
        'FontWeight','bold',...
        'FontSize',12,...
        'FitBoxToText','off',...
        'LineStyle','none');
    
    annotation(gcf,'textbox',...
        [0.421551312649173 0.206176220806794 0.16257756563245 0.0307855626326964],...
        'String',{'Null Distribution of Permuted Means'},...
        'FontWeight','bold',...
        'FontSize',12,...
        'FitBoxToText','off',...
        'LineStyle','none');
    
    annotation(gcf,'textbox',...
        [0.102028639618147 0.136112526539278 0.0772792362768496 0.0337855626326964],...
        'String',{'Observed Mean'},...
        'FontWeight','bold',...
        'FontSize',12,...
        'FitBoxToText','off',...
        'LineStyle','none');
    
    annotation(gcf,'textbox',...
        [0.125011933174232 0.26243949044586 0.105011933174224 0.0307855626326964],...
        'String',{'95% Acceptance Bands'},...
        'FontSize',12,...
        'FontWeight','bold',...
        'FitBoxToText','off',...
        'LineStyle','none');
    
    annotation(gcf,'arrow',[0.220167064439141 0.265513126491647],...
        [0.271823779193206 0.266454352441614]);
    
    annotation(gcf,'arrow',[0.220167064439141 0.266109785202864],...
        [0.267515923566879 0.200636942675159]);
    
    annotation(gcf,'arrow',[0.1736276849642 0.223150357995227],...
        [0.155050955414013 0.156050955414013]);
    
    
    return;
end

if nargin==1
    %Number of observations for each case
    trials=100;
    %Number of spatial locations
    locations=150;
    %number of iterations
    iterations=5000;
    %Plot result with traces
    ploton=1;
    %alpha
    alpha_level=.05;
    %x-values
    xvals=1:locations;
    
    if strcmpi(varargin{1},'var') || strcmpi(varargin{1},'variance')
        %Stat function
        statfcn=@(x)var(x,[],2);
        
        %Make constant condition
        p1=randn(locations,trials);
        p2=randn(locations,trials);
        
        varsize=[logspace(0,log(3),locations/2) logspace(log(3),0,locations/2)];
        for i=1:locations
            p2(i,:)=p2(i,:)*varsize(i);
        end
        
        %Smooth across location to simulate dependency
        for i=1:trials
            p1(:,i)=smooth(p1(:,i),20,'sgolay');
            p2(:,i)=smooth(p2(:,i),20,'sgolay');
        end
    else
        %Stat function
        statfcn=@(x)mean(x,2);
        
        %Make constant condition
        p1=randn(locations,trials);
        
        %Make variable condition
        p2=[];
        pos=[linspace(0,1,locations/2) linspace(1,0,locations/2)];
        for i=1:locations
            p2=[p2; randn(1,trials)+pos(i)];
        end
        
        %Smooth across location to simulate dependency
        for i=1:trials
            p1(:,i)=smooth(p1(:,i),20,'sgolay');
            p2(:,i)=smooth(p2(:,i),20,'sgolay');
        end
    end
elseif nargin>9
    error('Too many input arguments');
else
    defaults={'','',[],.05,@(x)mean(x,2),10000,'Group 1','Group 2',1};
    defaults(1:length(varargin))=varargin;
    
    p1=defaults{1};
    p2=defaults{2};
    xvals=defaults{3};
    
    alpha_level=defaults{4};
    statfcn=defaults{5};
    iterations=defaults{6};
    group1_name=defaults{7};
    group2_name=defaults{8};
    ploton=defaults{9};
end


locations=size(p1,1);
if isempty(xvals)
    xvals=1:locations;
end

%The difference in mean activity between the two scenarios
data_stat=statfcn(p1)-statfcn(p2);

%Combine both groups
all=[p1 p2];

Np1=size(p1,2);
% Np2=size(p2,2);

%Generate two groups boostrap samples from the shuffled distribution and compute the
%difference in means between them
diff_all=zeros(locations,iterations);
for i=1:iterations
    inds=randperm(size(all,2));
    inds1=inds(1:Np1);
    inds2=inds((Np1 + 1):end);
    diff_all(:,i)=(statfcn(all(:,inds1))-statfcn(all(:,inds2)));
end

%Sort the trials for ease of computing the pointwise percentile
sdiff_all=sort(diff_all,2);

%Set initial bounds at highest %ile
pbounds=[1 iterations];

% Number of individual trials falling outside the global bounds
numout=0;
if ~(any(any(sdiff_all)))
    gbounds=sdiff_all(:,pbounds);
else
    %Shrink the bounds until the number of outside trials is %5
    while numout<floor(iterations*alpha_level)+1
        %Calculate the global bounds
        gbounds=sdiff_all(:,pbounds);
        
        %Calculate the low and hig bounds
        lo=repmat(gbounds(:,1),1,iterations);
        hi=repmat(gbounds(:,2),1,iterations);
        
        %Check the number out of bounds
        numout=sum(any(~((diff_all>=lo) & (diff_all<=hi))));
        
        %Shrink the bounds
        pbounds(1)= pbounds(1)+1;
        pbounds(2)= pbounds(2)-1;
    end
end

%Find the significant bins
sigbins_all=~(data_stat>=gbounds(:,1) & data_stat<=gbounds(:,2));
acceptance_bounds=gbounds;
[cons_all,sig_regions]=consecutive(sigbins_all);

%Plot the results
if ploton
    figure1=figure('units','normalized','position',[0 0 1 1],'color','w');
    
    % Create axes
    %     ax(1) = axes('Parent',figure1,...
    %         'Position',[0.08 0.703333333333333 0.87 0.246666666666667]);
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
        %     title('Condition 1','fontsize',15);
        %     xlim(xvals([1 end]));
        %
        %     axes(ax(2))
        plot(xvals,p2,'b');
    end
    title([group1_name ' vs. ' group2_name],'fontsize',15);
    axis tight;
    xlim(xvals([1 end]));
    
    axes(ax(2))
    hold on;
    %Plot null distribution
%     plot(xvals,diff_all(1:10:end,:),'color',[1 .6 1]);
    %Plot global acceptance bounds
    plot(xvals,acceptance_bounds,'r','linewidth',2);
    %Plot data mean
    plot(xvals,data_stat,'b','linewidth',2);
    
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

