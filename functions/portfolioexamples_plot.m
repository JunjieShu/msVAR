function portfolioexamples_plot(varargin)
% portfolioexamples_plot - Specialized plotting function for portfolioexamples.
%
% This undocumented support function may be modified or removed in a future release.
% It is used by portfolioexamples.m in the Financial Toolbox.
%
% The function plots multiple line and scatter plots on a single figure.
% Any line plots to be included in the legend must be first arguments into this function.
% No error-checking takes place in this function.
%
% Usage:
%	First argument must be a title for the plot.
%	Remaining arguments are cell arrays for each plot overlay.
%
%	{ 'type', rsk, ret, label, style, line}, ...
%
%	'type'	'line'		line plot
%			'scatter'	scatter plot
%	rsk					risk proxies
%	ret					return proxies
%	label				labels for each series
%	style				style for plot (string to specify color, also linestyle if line plot)
%						default is 'b' for line plots and 'g' for scatter plots
%	line				if a line plot, the width of the line, default is 2

% Copyright 2011 The MathWorks, Inc.

plottitle = varargin{1};
plotlegend = [];

for i = 2:nargin
	plotinfo = varargin{i};
	
	plottype = plotinfo{1};
	plotrsk = plotinfo{2};
	plotret = plotinfo{3};
	if numel(plotinfo) > 3
		plotlabel = plotinfo{4};
	else
		plotlabel = [];
	end
	if numel(plotinfo) > 4
		plotstyle = plotinfo{5};
		if isempty(plotstyle)
			plotstyle = 'b';
		end
	else
		if strcmpi(plottype,'line')
			plotstyle = 'b';
		else
			plotstyle = 'g';
		end
	end
	if numel(plotinfo) > 5
		plotline = plotinfo{6};
		if isempty(plotline)
			plotline = 2;
		end
	else
		if strcmpi(plottype,'line')
			plotline = 2;
		else
			plotline = [];
		end
	end
	
	% line plot
	if strcmpi(plottype,'line')
		for k = 1:size(plotrsk,2)
			plot(plotrsk(:,k), plotret(:,k), plotstyle, 'LineWidth', plotline);
			if i == 2 && k == 1
				hold on
			end
			if ~isempty(plotlabel) && ~isempty(plotlabel{k})
				plotlegend = [ plotlegend, {plotlabel{k}} ];
			end
		end
		
	% scatter plot
	else
		if any(plotstyle == '.')
			scatter(plotrsk, plotret, plotstyle);
		else
			scatter(plotrsk, plotret, plotstyle, 'Filled');
		end			
		if i == 2
			hold on
		end
		if ~isempty(plotlabel)
			for k = 1:numel(plotrsk)
				if ~isempty(plotlabel{k})
					text(plotrsk(k) + 0.005, plotret(k), plotlabel{k}, 'FontSize', 8);
				end
			end
		end
	end
end

if ~isempty(plotlegend)
	legend(plotlegend,'Location','SouthEast');
end

title(['\bf' plottitle ]);
xlabel('Standard Deviation of Returns (Annualized)');
ylabel('Mean of Returns (Annualized)');
grid on

hold off
