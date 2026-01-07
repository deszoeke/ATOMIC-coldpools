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
% Transform accuracy test
% Test script for colorspace

% Pascal Getreuer 2006-2010


fprintf(['\nTransform accuracy test\n\n',...
      'To verify the invertibility of the color transfomations, this test\n',...
      'transforms sRGB data to a space, inverts, and compares with the\n',...
      'original data.\n']);
N = 1e5;            % Number of points to test
A = rand(N,3);      % Generate points uniformly in the sRGB colorspace

% Include pure black and pure white
A(1,:) = 0;
A(2,:) = 1;

Space = {'YPbPr', 'YCbCr', 'JPEG-YCbCr', 'YDbDr', 'YIQ','YUV', 'HSV', ...
      'HSL', 'HSI', 'XYZ', 'Lab', 'Luv', 'LCH', 'CAT02 LMS'};
fprintf('\n Transform          RMSE Error   Max Error\n\n');

for k = 1:length(Space)
   B = colorspace([Space{k},'<-RGB'],A);  % Convert to Space{k}
   R = colorspace(['RGB<-',Space{k}],B);  % Convert back to sRGB
   RMSE = sqrt(mean((A(:) - R(:)).^2));
   MaxError = max(abs(A(:) - R(:)));
   fprintf(' RGB<->%-10s   %9.2e    %9.2e\n', Space{k}, RMSE, MaxError);
end

fprintf('\n\n');
