function varargout = dotPlot(S)
% function varargout = dotPlot(S)
% dotPlot - Doesn't plot box plots!
%
% based on: function notBoxPlot(y,x,jitter)
%
%
% Purpose
% An alternative to a box plot, where the focus is on showing raw
% data. Plots columns of y as different groups located at points
% along the x axis defined by the optional vector x. Points are
% layed over a 1.96 SEM (95% confidence interval) in red and a 1 SD
% in blue. The user has the option of plotting the SEM and SD as a
% line rather than area. Raw data are jittered along x for clarity. This
% function is suited to displaying data which are normally distributed.
% Since, for instance, the SEM is meaningless if the data are bimodally
% distributed.
%
%
% Inputs: a structure with the following fields (only y is mandatory):
% y - each column of y is one variable/group. If x is missing or empty
%     then each column is plotted in a different x position.
%
% x - optional, x axis points at which y columns should be
%     plotted. This allows more than one set of y values to appear
%     at one x location. Such instances are coloured differently.
% Note that if x and y are both vectors of the same length this function
% behaves like boxplot (see Example 5).
%
% jitter - how much to jitter the data for visualization
%          (optional). The width of the boxes are automatically
%          scaled to the jitter magnitude. If jetter is empty or
%          missing then a default value is used.
%
% groupingVar - 1 number per subject to determine subject groups
% colors - colors for display [0,1] 3 times
% colorContourBlack - 1 makes color of contour of dots black, 0 makes it same as dot color
% groupColors - different colors for each group
% markersize - default is 8
% MUmarkers - markers for the mean of each set of dots; can be 'o' (default), 'bar', 'errorbar', or 'change' (changes with subject group)
% MUmarkerSize - default = markersize*2
% errorBarWidth - default = 5
%
% Outputs
% H - structure of handles for plot objects. Only returned if an
%     output argument is requested.
%
% Example 1 - example with dots as central tendency marker
% figure
% S.y = randn(20,5); dotPlot(S);
%
% Example 2 - example with bars as central tendency marker
% figure
% S.y = randn(20,5); S.MUmarkers = 'bar'; dotPlot(S);
%
% Example 3 - example with additional grouping variable that assigns
% rows to 2 different groups
% figure
% S.y = randn(10,2); S.groupingVar = [1 1 1 2 2 1 2 1 2 2]; dotPlot(S);
%
% Rob Campbell - January 2010
% Johannes Schultz - Oct 2018
%
% also see: boxplot

% Get variables from input structure S, or give them default values
if ~isstruct(S); y=S; clear S; S.y=y; end % in case the user is lazy and just inputs data
if isfield(S,'y'); y = S.y; else; error('Input data needed as S.y'); end
if isfield(S,'x') && ~isempty(S.x); x = S.x; else, x=1:size(y,2); end
if isfield(S,'jitter'); jitter = S.jitter; else
  try jitter = (x(2)-x(1))/3; catch
    jitter=0.3; %larger value means greater amplitude jitter
  end
end
if isfield(S,'groupingVar'); groupingVar = S.groupingVar; end
if isfield(S,'colors'); colors = S.colors; end
if isfield(S,'colContourBlack'); colContourBlack = S.colContourBlack; else, colContourBlack = 0; end
if isfield(S,'groupColors'); groupColors = S.groupColors; end
if isfield(S,'markersize'); markersize = S.markersize; else markersize = 8; end
if isfield(S,'legend'); legendTxt = S.legend; end
if isfield(S,'MUmarkers'); MUmarkers = S.MUmarkers; else MUmarkers = 'errorline'; end % can be 'o', 'bar', 'errorbar', 'errorline', 'errorlineblack', 'change' or the markers you want
if isfield(S,'MUmarkers'); MUmarkers = S.MUmarkers; else MUmarkers = 'errorbar'; end % can be 'o', 'bar', 'errorbar', 'errorline', 'errorlineblack', 'change' or the markers you want
if isfield(S,'MUmarkerSize'); MUmarkerSize = S.MUmarkerSize; else MUmarkerSize = markersize*2; end
if isfield(S,'errorBarWidth'); errorBarWidth = S.errorBarWidth; else errorBarWidth = 5; end
if isfield(S,'errorType'); errorType = S.errorType; else errorType = 'SEM'; end

% For future developments:
% median or mean
% SEM or STD or CI95 or like boxplot notches
% boxplotin front/behind or not

% Check input arguments
if nargin==0
  help(mfilename)
  return
end

if isvector(y), y=y(:); end

%If x is logical then the function fails. So let's make sure it's a double
x=double(x);

if jitter==0
  warning('A zero value for jitter means no patch object visible')
end

NpGrp = size(y,2);

% 3rd dimension: different groups based on grouping variable
if exist('groupingVar','var') && ~isempty(groupingVar)
  grps = unique(groupingVar); Ngrps = length(grps);
  if strcmpi(MUmarkers,'change')
    MUmarkers = 'ooddssvv^^pphh';
  end
  for i = 1:Ngrps
    Ns(i) = length(find(groupingVar==grps(i)));
  end
  data = NaN(max(Ns),NpGrp*Ngrps);
  for i = 1:Ngrps
    idx = [1:NpGrp] + NpGrp * (i-1);
    data(1:Ns(i),idx) = y(find(groupingVar==grps(i)),:);
  end
  y = data;
  origx = x;
  x = repmat(x,1,Ngrps) + kron([0:Ngrps-1]*NpGrp,ones(1,NpGrp));
else Ngrps = 1; if strcmpi(MUmarkers,'change'); MUmarkers = 'o'; end
end

% 3rd dimension: different groups based on 3rd dim of input variable
if size(y,3)>1
  Ngrps = size(y,3);
  grps = [0:Ngrps-1]';
  data = [];
  for i = 1:Ngrps
    data = [data y(:,:,i)];
  end
  y = data;
  origx = x;
  x = repmat(x,1,Ngrps) + kron([0:Ngrps-1]*NpGrp,ones(1,NpGrp));
  if strcmpi(MUmarkers,'change')
    MUmarkers = 'ooddssvv^^pphh';
  end
end

if isvector(y) && isvector(x) && length(x)>1
  x=x(:);
  
  if length(x)~=length(y)
    error('length(x) should equal length(y)')
  end
  
  u=unique(x);
  for ii=1:length(u)
    f=x==u(ii);
    h(ii)=notBoxPlot(y(f),u(ii),jitter);
  end
  
  %Make plot look pretty
  if length(u)>1
    xlim([min(u)-1,max(u)+1])
    set(gca,'XTick',u)
  end
  
  if nargout==1
    varargout{1}=h;
  end
  
  return
  
end


if length(x) ~= size(y,2)
  error('length of x doesn''t match the number of columns in y')
end


%We're going to render points with the same x value in different
%colors so we loop through all unique x values and do the plotting
%with nested functions. No clf in order to give the user more
%flexibility in combining plot elements.
hold on
[uX,~,b]=unique(x);
h=[];

% determine colors
if ~exist('colors','var') || isempty(colors)
  if Ngrps==1
    Ncols = length(unique(x));
  else
    Ncols = length(unique(origx));
  end
  if Ncols<12
    cols = niceColorsMany(Ncols);
  else
    cols = niceColors(Ncols);
  end
else
  cols = colors;
end
if Ngrps>1 % there are groups, either repeat colors or replace with group colors
  if exist('groupColors','var') && ~isempty(groupColors)
    cols = repmat(groupColors,NpGrp,1);
  else
    cols = repmat(cols,Ngrps,1); % repeat colors within each group
  end
end
for ii=1:length(uX)
  f=b==ii;
  if ~any(strncmp(MUmarkers,{'bar';'errorbar';'errorline';'errorlineblack';'change'},10))
    MUmarker = MUmarkers(ceil(ii/NpGrp));
    if rem(ceil(ii/NpGrp),2) % by default mu marker color is white
      MUcol = [1 1 1];
    else
      MUcol = cols(ii,:); % every 2nd group, make mu marker color same as dots
    end
  else
    MUmarker = MUmarkers;
    MUcol = [1 1 1];
  end
  h = [h,myPlotter(x(f),y(:,f),cols(ii,:),colContourBlack,MUmarker,MUmarkerSize,MUcol,errorBarWidth,errorType)];
end

hold off

%Tidy up plot: make it look pretty
if length(x)>1
  set(gca,'XTick',unique(x))
  xlim([min(x)-1,max(x)+1])
end

if nargout==1
  varargout{1}=h;
elseif nargout==2
  varargout{1}=h;
  varargout{2}=y; % only helpful when using grouping variable
end

if exist('legendTxt','var')
  legend([h(:).data],legendTxt,'location','best')
end

% Nested functions follow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function h=myPlotter(X,Y,col,colContourBlack,MUmarker,MUmarkerSize,MUcol,errorBarWidth,errorType)
    
    mu=nanmean(Y); %Requires the stats toolbox
    switch errorType
      case 'SEM'
        % err=SEM_calc(Y); %Supplied external function
        err = squeeze(nanstd(Y))' / sqrt(length(Y));
        errL = err; errH = err;
      case 'CI'
        CI = ci_t(Y,95); % find 95% confidence intervals for the mean, assuming normally-distributed data
        errL = mu - squeeze(CI(1,:,:))';
        errH = squeeze(CI(2,:,:))' - mu;
    end
    % SD=nanstd(Y);  %Requires the stats toolbox
    
    %The plot colors to use for multiple sets of points on the same x
    %location
    %  cols=hsv(length(X)+1)*0.5;
    %  cols(1,:)=0;
    % jitScale=jitter*0.55; %To scale the patch by the width of the jitter
    
    %     if strcmpi(MUmarker,'bar') % plot mean as bar
    %       h_sem = errorbar(X,mu,SEM,'k');
    %     end
    
    % Plot mean and SEM as bar for mean & thinner bar for SEM
    for k=1:length(X)
      if strcmpi(MUmarker,'bar') % plot mean as bar
        h(k).sem=plot([X(k),X(k)],[mu(k)-errL(k),mu(k)+errH(k)],'-r',...
          'linewidth',errorBarWidth,'color',[0 0 0]);
        h(k).mu=bar(X(k),mu(k),'facecolor',col);
      elseif strcmpi(MUmarker,'errorbar') % plot mean as bar
        h(k).mu=bar(X(k),mu(k),'facecolor',col); % plot bar now, errorbar after the dots
        %       elseif strcmpi(MUmarker,'errorline') % plot mean as line, SEM as lines, all in condition-specific color
        %         h(k).mu=line([X(k)-.5 X(k)+.5],[mu(k) mu(k)],'color',col,'linewidth',2);
        %         h(k).sem=errorbar(X(k),mu(k),SEM(k),'color',col,'linewidth',2,'CapSize',18);
        %       elseif strcmpi(MUmarker,'errorlineblack') % plot mean as black line
        %         h(k).mu=line([X(k)-.5 X(k)+.5],[mu(k) mu(k)],'color',[0 0 0],'linewidth',2);
        %         h(k).sem=errorbar(X(k),mu(k),SEM(k),'.k','linewidth',2,'CapSize',18);
      end
      
      thisY=Y(:,k);
      thisY=thisY(~isnan(thisY));
      thisX=repmat(X(k),1,length(thisY));
      
      %Plot jittered raw data
      %       if strcmpi(MUmarker,'errorlineblack') % plot dots in black
      %         C = [.5 .5 .5];
      %       else
      C = col+([1 1 1] - col)/2; % lighten up the color
      if colContourBlack; Cb = [0 0 0]; % marker contour color
      else, Cb = C;
      end
      
      J = (rand(size(thisX))-0.5)*jitter;
      h(k).data=plot(thisX+J, thisY, 'o', 'color', Cb,...
        'markerfacecolor', C,...
        'markersize',markersize);
      
      if strcmpi(MUmarker,'errorline') % plot mean and SEM as lines in condition-specific color
        h(k).mu=line([X(k)-.5 X(k)+.5],[mu(k) mu(k)],'color',col,'linewidth',2); % here instead of col I could use [0 0 0] for black mean errorline
        h(k).sem=errorbar(X(k),mu(k),errL(k),errH(k),'color',col,'linewidth',2,'CapSize',18);
      elseif strcmpi(MUmarker,'errorbar') % plot mean as bar
        h(k).sem=errorbar(X(k),mu(k),errL(k),errH(k),'.k','linewidth',2,'CapSize',18);
      elseif strcmpi(MUmarker,'errorlineblack') % plot mean and SEM as black lines
        h(k).mu=line([X(k)-.5 X(k)+.5],[mu(k) mu(k)],'color',[0 0 0],'linewidth',2);
        h(k).sem=errorbar(X(k),mu(k),errL(k),errH(k),'.k','linewidth',2,'CapSize',18);
      end
    end
    
    % Plot mean and SEM as dot & bar
    if ~any(strncmp(MUmarker,{'bar','errorbar','errorline','errorlineblack'},10))
      for k=1:length(X)
        h(k).sem=plot([X(k),X(k)],[mu(k)-errL(k),mu(k)+errH(k)],'-r',...
          'linewidth',errorBarWidth,'color',col);
        h(k).xAxisLocation=x(k);
        h(k).mu=plot(X(k),mu(k),'o','color',col,...
          'markerfacecolor',MUcol,...
          'markersize',MUmarkerSize,...
          'linewidth',2);
      end
    end
    %     h(:).sem = h_sem;
  end %function myPlotter



end %function notBoxPlot
