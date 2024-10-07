function varargout = printSVG(varargin)
% prints current figure to SVG file at (pretty much) same size as on
% screen. If 
% 
% INPUTS: 4 (all optional)
% - 1st argument: name, default = the one from figure
% - 2nd argument: file format, default = 'png'
% - 3rd argument: resolution, default = 300 DPI (irrelevant for eps formats)
% - 4th argument: image width in cm, default = 20. Value defines
%                       width of figure, scales height to keep aspect ratio
% 
% EXAMPLE:
% figure; plot([1:10]);
% printfig('test','jpeg',200,10)  % will create a .jpg file with 200dpi,
%                                 % 10cm width (actually: equivalent pixel 
%                                 % number only), called "test.jpg"
% or: printfig([],'epsc')
%
% js - july 2008

% FOR SCREEN SIZE / PRINTED SIZE PROBLEM, SEE:
% 
% I want to set the size of the figure in cm, but keep all proportions the
% same in print. That's not easy: If I change the paperposition values
% through command line (e.g.: set(gcf,'paperPosition',[.2 .5 .4 3])),
% PaperPositionMode goes to "manual", the printed output is of
% the new specified size, but nothing changes on screen!! And Fontsize and
% linewidth are kept as was, so do not scale with figure size.
% HOW to get scaled fontsizes after rescaling from command line? Maybe it says here:
% http://www.mathworks.com/support/solutions/data/1-P51QW.html?solution=1-P51QW

% first get figure settings for size:
set(gcf,'PaperPositionMode','auto','papertype','A4','paperUnits','centimeters')
sz=get(gcf,'paperposition'); sz=sz(3:4);
inch = 2.5381;

% set defaults:
arg{1} = get(gcf,'name');       % file name
if isempty(arg{1}); f=gcf; arg{1} = ['Figure' num2str(f.Number)];  end
arg{2} = 'svg';                 % file format
arg{3} = 300;                   % resolution (dpi)
arg{4} = 20;                     % figure width in cm

% replace defaults by inputs given:
for n = 1:nargin
    if ~isempty(varargin{n})
        arg{n} = varargin{n};
    end
end

% make name filename-compatible:
arg{1} = strrep(arg{1},':','-');
arg{1} = strrep(arg{1},'\','-');
arg{1} = strrep(arg{1},'/','-');

% use meaningful variable names:
[name,fileformat,fileresolution,imagewidth] = deal(arg{:});

if ~strcmpi(name(end-2:end),fileformat)
  name = [name '.' fileformat];
end

% effective file resolution used by print.m to get the right image size in
% pixels, keeping fontsizes and everything else as on screen. Haven't found
% a way to set the printed file size in cm to what is desired, only pixel
% number is ok. Use GIMP, Photoshop or don't know what to change the image
% resolution / size
fileresolutionEff = imagewidth/sz(1) * fileresolution;

% give file info:
npix = round(([imagewidth imagewidth/sz(1)*sz(2)]/inch)*fileresolution);
disp(sprintf('Creating image "%s" with about %d by %d pixels, for size of %.2f by %.2f cm at resolution of %d dpi',...
    name,npix(1),npix(2),imagewidth,imagewidth/sz(1)*sz(2),fileresolution))

% now just print it
print(gcf,['-d' fileformat],['-r' num2str(fileresolutionEff)],name)
% print(gcf,['-d' fileformat],['-r' num2str(fileresolutionEff)],'-zbuffer',name)
    