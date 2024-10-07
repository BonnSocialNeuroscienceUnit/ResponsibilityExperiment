function [oo,ooo] = meansembar(varargin)
% function hh = meansembar(varargin)
% does barplots with SEM of input matrix.
%
% INPUTS
% if one input: the data. That's the minimum input, rest is optional.
% if >1:
%   input 1: the x positions of the bars or sets of bars
%   input 2: the data matrix (1 to 3 dimensions, second dimension determines number of bars per set of bars)
%   input 3: the width of the bars (default = 0.8)
%   input 4: the color of the bars (default is defaultColorOrder)
%   input 5: the fontsize (default is 14, looks nice)
%   input 6: the width of the errorbars (default is .2, looks nice)
%   input 7: the type of error: 'SEM' or 'CI' (default is SEM)
% OUTPUTS
%   output 1: handles to bars
%   output 2: handles to errorbars
%
% EXAMPLE:
% meansembar([1 2],rand(10,3,2),1,[1 0 0;0 1 0;0 0 1]); % plots 2 sets of 3 bars (red, green, blue) with errorbars nicely centered.
% meansembar(rand(12,4)) %does barplot of mean(rand(12,4)) + SEM of rand(12,4), in first default color

fontSize = 14;
colors = repmat(get(gcf,'defaultAxesColorOrder'),20,1);
width = .8;
errBarWidth = .2;
errorType = 'SEM';

if nargin == 1
  data = varargin{1};
else
  x = varargin{1};
  data = varargin{2};
  try width = varargin{3}; end
  try colors = varargin{4}; end
  try fontSize = varargin{5}; end
  try errBarWidth = varargin{6}; end
  try errorType = varargin{7}; end
end
m = squeeze(nanmean(data))';

if ~exist('x','var')
  x = 1:size(m,1);
end

switch errorType
  case 'SEM'
    err = squeeze(nanstd(data))';
    for aa=1:size(data,2)
      for bb=1:size(data,3)
        for cc=1:size(data,4)
          Ns(aa,bb,cc) = length(find(~isnan(data(:,aa,bb,cc))));
        end
      end
    end
    try err = err ./ sqrt(Ns); catch; err = err ./ sqrt(Ns'); end
    errL = err; errH = err;
  case 'CI'
    CI = ci_t(data,95); % find 95% confidence intervals for the mean, assuming normally-distributed data
    errL = m - squeeze(CI(1,:,:))';
    errH = squeeze(CI(2,:,:))' - m;
end

% plot bar
hh = bar(x,m,width);

pause(0.1); %pause allows the figure to be created

% set bar colors
for r = 1:length(hh)
  set(hh(r),'facecolor',[colors(r,:)])
end

hold on
set(gca,'fontsize',fontSize)

% for each bar, add errorbar
Nbars = size(m,2);
% ErrorbarPos = [(Nbars-1):-2:-(Nbars-1)] * .14;
for r = 1:size(m,2)
  try H(r) = get(hh(r),'children');
    % plot errorbars
    hhh(:,r) = errorbar(hh(r).XData+hh(r).XOffset,m(:,r),errL(:,r),errH(:,r),'.k');
    
    ch = get(hhh(r),'children');     Hbars = reshape(get(ch(2),'xdata'),3,[]);
    barwidth = nanmean(diff(Hbars(:,2)));
    HbarWidth = barwidth/2*errBarWidth/.115;
    %         offsets = kron(ones(length(x)*3,1),[HbarWidth -HbarWidth]);
    Hbars(:,[2:3:end]) = [Hbars(1,[1:3:end])-HbarWidth; Hbars(1,[1:3:end])+HbarWidth; repmat(NaN,1,size(Hbars,2)/3)];
    Hbars(:,[3:3:end]) = [Hbars(1,[1:3:end])-HbarWidth; Hbars(1,[1:3:end])+HbarWidth; repmat(NaN,1,size(Hbars,2)/3)];
    %         newxdata = [repmat(nanmean(reshape(Hbars,3,[]))',1,2)+offsets NaN(size(offsets,1),1)]';
    %         set(ch(2),'linewidth',1,'xdata',newxdata(:)');
    set(ch(2),'linewidth',1,'xdata',Hbars(:)');
    
  catch  
    hhh(r) = errorbar(hh(r).XData+hh(r).XOffset,m(:,r),errL(:,r),errH(:,r),'.k');
    %         hhh(r) = errorbar(hh(r).XData-ErrorbarPos(r),m(:,r),se(:,r),se(:,r),'.k');
  end
  set(hhh,'markersize',1,'linewidth',1)
end

set(gca,'xlim',[0.2 x(end)+0.8])
ylabel(['Mean +- SEM, N = ' num2str(size(data,1))])
hold off

if nargout>0, oo = hh; ooo = hhh; end
