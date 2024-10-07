function colors = niceColorsWarm2(Ncolors)
% just a simple one-liner to produce up to 4 colors for plotting. Uses
% predefined colors from a nice figure in a paper

if nargin == 0
  figure; meansembarFigure(rand(6),'color',niceColorsWarm2(6))
  return
end

if Ncolors > 6
    error('Only 6 colors implemented')
end

darkBlue = [5 82 147];
green = [45 159 44];
orange = [255 127 14];
red = [213 1 4];
darkgreen = (darkBlue + green) / 2;
redorange = (red + orange) / 2;

niceColors = [...
    darkBlue;...
    green;...
    orange;...
    red;...
    darkgreen;...
    redorange;...
    ]/255;

colors = niceColors(1:Ncolors,:);

% % to see them:
% figure; meansembarFigure(rand(4),'color',niceColorsWarm2(4))
