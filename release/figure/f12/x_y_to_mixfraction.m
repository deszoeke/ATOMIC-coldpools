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
function [fen, fss, fdd] = x_y_to_mixfraction(yml,xml, yen,xen, yss,xss, ydd,xdd)
% [fen, fss, fdd] = th_q_to_mixfraction(thml,qml, then,qen, thss,qss, thdd,qdd)
% Diagnoses BL air as a 3-part mixture fen + fss + fdd = 1
% of entrained air, air in equilibrium with the sea surface,
% and saturated downdraft air.
% The mixing fractions of these 3 end members are inverted algebraically
% from observed thml, qml as in de Szoeke (2018).
%
% (c) 2020-08-19 Simon P. de Szoeke 

% xen(xen>30)=NaN;
% solves algebraically for mixing fractions fen,fss,fdd at each time
kay=(yss-ydd)./(xdd-xss);
beta=ydd+kay.*xdd;
gamma=yen+kay.*xen-beta;
fen=(yml+kay.*xml-beta)./gamma;
fss=(xml-xdd-(xen-xdd).*fen)./(xss-xdd);
fdd=1-fen-fss;

