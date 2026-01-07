% Copyright 2025 Simon P. de Szoeke and Estefanía Quiñones 
% Meléndez.
% 
% Permission is hereby granted, free of charge, to any person 
% obtaining a copy of this software and associated documentation 
% files (the “Software”), to deal in the Software without 
% restriction, including without limitation the rights to use, 
% copy, modify, merge, publish, distribute, sublicense, and/or 
% sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following 
% conditions:
%
% The above copyright notice and this permission notice shall be 
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.
function cmap = b2rcolormap(clevs)
% Sets the colormap to a high-contrast nonclashing
% blue-red colormap. clevs is the array of contour levels, or the
% scalar number of levels. (Usually, enter the length
% of the colormap + 1.)
%
% Simon de Szoeke
% 1 Dec 2005
% Thanks to Todd Mitchell at JISAO for the blue colormap.

%clevs = 10:10:100; %set your own

if length(clevs)==1
  ncols = round(clevs)-1;
else
  ncols = length(clevs)-1;
end
b = [230 255 255; 160 240 255; 80 180 255; 30 110 250; 10 50 200; 10 50 120 ]/255; 
r = fliplr(b);
n = length(b)+length(r);
cmap = interp2((1:3)',(0:n-1)/(n-1),[flipud(b); r],(1:3)',(0:ncols-1)/(ncols-1),'linear');
colormap(cmap);
