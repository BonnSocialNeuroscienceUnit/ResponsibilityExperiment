function [hdlMeanSem,hdlDots,hdlLeg] = meanErrPlotDots(data,varargin)
% function [hdlMeanSem,hdlDots,hdlLeg] = meanErrPlotDots(data,varargin)
% 
% Plots mean, all datapoints, and error in 3 possible forms: SEM, STD and
% 95% CI.

% Display variables
marker = '__os^v><dph*';
linespec = 'ob';
errorType = 'CI'; % alternatives: 'SEM', 'STD'
markersize = 6; % 12;
markerfacecolor = [1 1 1];
dotSize = 16;
dotColor = [.68 .68 .68];
xoffset = .15;
linewidth = 2;
capSize = 0; % the size of the little horizontal lines at the end of thh error bars

% Get the data
% data = varargin{1};
x = 1:size(data,2);

% set / get potential display variables passed as arguments, defaults are specified above
for n = 1:nargin-2 % because the last argument will be a value, not a variable type
  if ischar(varargin{n}) % if there's a text argument, it's considered one of the display variables above
    eval([varargin{n} ' = varargin{n+1};']) % the variable named in varagin{n} is set to the value coming as next argument
  end
end

% 
% for v = 1:narg
%   if ischar(varargin{v})
%     errorType = varargin{v};
%     idx = find([1:narg]~=v);
%     varargin = varargin(idx);
%   end
% end
% narg = length(varargin);
% 
% nocol = zeros(1,narg);
% for v = 1:narg
%   if ischar(varargin{v})
%     linespec = varargin{v};
%   else
%     nocol(v) = 1;
%   end
% end
% if sum(nocol) ~= narg
%   if narg == 3
%     x = varargin{1};
%     data = varargin{2};
%     % 3rd input is color
%     xoffset = varargin{4};
%   elseif narg == 3
%     x = varargin{1};
%     data = varargin{2};
%     % 3rd input is color
%   elseif narg == 2
%     data = varargin{1};
%     x = 1:size(data,2);
%   end
% else
%   if narg == 2
%     x = varargin{1};
%     data = varargin{2};
%   elseif narg == 1
%     data = varargin{1};
%     x = 1:size(data,2);
%   end
% end

m = squeeze(nanmean(data,1))';

switch errorType
  case 'SEM'
    err = squeeze(nanstd(data))' / sqrt(size(data,1));
    errL = err; errH = err;
    
  case 'CI'
    CI = ci_t(data,95); % find 95% confidence intervals for the mean, assuming normally-distributed data
    errL = m - squeeze(CI(1,:,:))';
    errH = squeeze(CI(2,:,:))' - m;
  case 'STD'
    err = squeeze(nanstd(data))';
    errL = err; errH = err;
end

if size(m,2) == 1
  m = m'; errL = errL'; errH = errH';
end

% s = size(m,2);
%
% approxOffset = 1/(s + 1.5 + s/10 - 1/(s*5));
% approxOffsets = cumsum(repmat(approxOffset,1,size(m,2)));
% approxOffsets = approxOffsets-mean(approxOffsets);


orighold = ishold;
if ~orighold
  cla
end
hold on
if size(m,1) > 1
  xd = linspace(-xoffset,xoffset,size(data,3));
  cols = get(gca,'colororder');
  cols = [cols; cols/2];
  if size(m,1)<=4
    cols = niceColorsDark(size(m,1));
  elseif size(m,1)<=7
    cols = niceColorsWarm2(size(m,1));
  else
    cols = niceColors(size(m,1));
  end
  markerfacecolor = cols;
  for r = 1:size(m,1)
    Xnoise = ( rand(size(data,1), size(m,2))-.5 ) / 10;
    hdlDots(r,:) = plot(repmat(x+xd(r),1*size(data,1),1) + Xnoise,data(:,:,r),'.k','markersize',dotSize,'color',dotColor); hold on
    hdlMeanSem(r) = errorbar(x+xd(r), m(r,:), errL(r,:), errH(r,:), linespec,'markersize',markersize, 'marker',marker(r),'linewidth',linewidth,'markerFaceColor',markerfacecolor(r,:),'color',cols(r,:),'CapSize',capSize);
  end
  hdlLeg = legend([hdlMeanSem(1,:)],num2str([1:r]'));
else
  Xnoise = ( rand(size(data,1),size(data,2))-.5 ) / 8;
  hdlDots = plot(repmat(x,1*size(data,1),1) + Xnoise,data,'.k','markersize',dotSize,'color',dotColor); hold on
  hdlMeanSem = errorbar(x, m, errL, errH, linespec,'markersize',markersize, 'marker',marker(1),'linewidth',linewidth,'markerFaceColor',markerfacecolor,'CapSize',capSize);
end

set(gca,'xlim',[x(1)-0.8 x(end)+0.8],'xtick',x)
ylabel(['Mean +- SEM, N = ' num2str(size(data,1))])
if ~orighold
  hold off
end
% if nargout>0, hdlMeanSem = hhh; end
