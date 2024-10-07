function [str,H,p,STATS,stars] = ttest_nice(a,b,type,plotFlag,fontsize4text,printOutFlag,statsInTitleFlag,inputNames)

if ~exist('plotFlag','var') || isempty(plotFlag)
  plotFlag = 0;
end
if ~exist('printOutFlag','var') || isempty(printOutFlag)
  printOutFlag = 1;
end
if ~exist('statsInTitleFlag','var') || isempty(statsInTitleFlag)
  statsInTitleFlag = 1;
end
if ~exist('inputNames','var')  || isempty(inputNames)
  inputNames = {'A','B'};
end
if ~exist('b','var') || isempty(b)
  [H,p,CI,STATS] = ttest(a);
  [d,r] = cohens_d(a);
  type = '1-sample';
else
  if isempty(b)
    [H,p,CI,STATS] = ttest(a);
    [d,r] = cohens_d(a);
    type = '1-sample';
  end
  if ~exist('type','var')
    type = 'paired';
  end
  switch type
    case 'paired'
      [H,p,CI,STATS] = ttest(a,b);
      [d,r] = cohens_d(a,b,1);
    case '1-sample'
      [H,p,CI,STATS] = ttest(a,b);
      [d,r] = cohens_d(a,b,0);
    case '2-sample'
      [H,p,CI,STATS] = ttest2(a,b);
      [d,r] = cohens_d(a,b,0);
  end
end
STATS.cohensd = d;

% give stars to p values
for i = 1:length(p)
  if p(i) < 0.001; stars{i} = '***';
  elseif p(i) < 0.01; stars{i} = '**';
  elseif p(i) < 0.05; stars{i} = '*';
  elseif p(i) < 0.1; stars{i} = '^';
  else stars{i} = ' ';
  end
end
if ~exist('fontsize4text','var')
  fontsize4text = 14;
end
for i = 1:length(p)
  %     if p<0.001
  %         str{i} = sprintf('t(%d)=%.2f, p<0.001 (%s, 2-tailed)',STATS.df(i),STATS.tstat(i),type);
  %     else
  str{i} = sprintf('t(%d)=%.2f, p=%.4f, d=%.2f (%s, 2-tailed)',STATS.df(i),STATS.tstat(i),p(i),d(i),type);
  %     end
  if nargout==0 & printOutFlag
    disp(str{i})
  end
end
str = char(str{:});

doPlot = 0;
if plotFlag == 1; figure; doPlot = 1; end
if plotFlag ~= 0 && ishandle(plotFlag); axes(gca); doPlot = 1; end
if plotFlag == 0; doPlot = 0; end
if doPlot
  if ~isempty(b)
    if min(size(a)) == 1
      a = a(:); b = b(:);
    end
    switch type
      case 'paired'
        %         meansembarFigure([a b a-b],'fontsize',fontsize4text,'xticklabel',{inputNames{1},inputNames{2},'diff'}); starPos = 3;
        %         dotPlot([a b a-b]); set(gca,'xticklabel',[inputNames 'diff']); starPos = 3;
        hh = distributionPlot([a b],'color',[.8 .8 .8],'showMM',4); hold on;
        hh{2}(1).LineWidth=2; hh{2}(2).LineWidth=2; hh{2}(1).Marker='.';
        plot([a b]','k'); set(gca,'xticklabel',inputNames); starPos = 1.5;
        
      case '1-sample'
        hold on
        plot([0 size(a,2)+1],[b b],'k')
        %         meansembarFigure(a,'fontsize',fontsize4text,'color',[.5 .5 .5]); starPos = 3;
        dotPlot2(a); starPos = 1.5;
        for i = 1:length(p)
          if H(i)
            tt = text(i,.9*max(get(gca,'Ylim')),'*','fontsize',fontsize4text*2);
            set(tt,'horizontalAlignment','center')
          end
        end
      case '2-sample'
        temp = NaN(max([length(a) length(b)]),2);
        temp(1:length(a),1) = a;
        temp(1:length(b),2) = b;
        %         meansembarFigure(temp,'fontsize',fontsize4text,'xticklabel',{inputNames}); starPos = 1.5;
        dotPlot2(temp); set(gca,'xticklabel',inputNames); starPos = 1.5;
    end
    if statsInTitleFlag
      title(str)
    end
    % print a star when significant:
    if H
      if p<0.001
        %                 str{i} = sprintf('t(%d)=%.2f, p<0.001 (%s, 2-tailed)',STATS.df(i),STATS.tstat(i),type);
        tt(1) = text(starPos,.9*max(get(gca,'Ylim')),['<0.001'],'fontsize',fontsize4text);
      else
        tt(1) = text(starPos,.9*max(get(gca,'Ylim')),sprintf('%.3f',p),'fontsize',fontsize4text);
      end
      %             text(2.925,.9*max(get(gca,'Ylim')),'*','fontsize',36);
      tt(2) = text(starPos,.95*max(get(gca,'Ylim')),stars{1},'fontsize',round(fontsize4text*2));
      set(tt,'horizontalAlignment','center')
    end
  elseif isempty(b)
    %     meansembarFigure(a,'fontsize',fontsize4text); % set(gca,'xticklabel',{'A'})
    dotPlot2(a); % set(gca,'xticklabel',inputNames); starPos = 1.5;
    for i = 1:length(p)
      if H(i)
        tt = text(i,.9*max(get(gca,'Ylim')),stars{i},'fontsize',fontsize4text*2);
        set(tt,'horizontalAlignment','center')
        %                 text(1/i,.1,'*','sc')%,'fontsize',36);
      end
    end
    if statsInTitleFlag
      title(str)
    end
  end
end