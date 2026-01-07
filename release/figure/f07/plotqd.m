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
function plotqd(q_surf_b, qi_surf, q_ent_b, qi_ent, n)
% plotqd(q0,qi0, q1,qi1, n)
% q - deltaD plot mixing line between end members 0,1.

    Rvsmow = 155.76e-6; % unitless, deuterium

    % isotope plotting helper functions
    delta = @(qi,q) qi./(Rvsmow*q) - 1;
    mixt = @(q1,q2,n) q1:(q2-q1)/n:q2;

    qx = mixt(q_ent_b,q_surf_b,20);
    dy = 1e3*delta( mixt(qi_ent,qi_surf, n), mixt(q_ent_b,q_surf_b, n) );
    
    plot( qx, dy, 'k-')
    plot( qx(1), dy(1), 'ko', 'MarkerEdgeColor','k', 'MarkerFaceColor','k') % entrainment
    plot( qx(end), dy(end), 'ko', 'MarkerEdgeColor','k', 'MarkerFaceColor',1+[0, 0, 0]) % surface
end