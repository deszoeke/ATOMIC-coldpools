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
function inset_plot(i,tcold,Taf,cp_matrix,qairf,dDf,rrf)
    axes('Position',[.69 .55 .15 .15])
    box on
    plot(tcold(i,1:61),Taf(cp_matrix(i,1:61)),'.');
    datetick('x','HH:MM','keeplimits','keepticks')
    ylabel('Ta [C]')
    title(['CP #',num2str(i),'; onset on ',datestr(tcold(i,1))])
%     ylim([22.5 28])

    axes('Position',[.69 .35 .15 .15])
    box on
    plot(tcold(i,1:61),qairf(cp_matrix(i,1:61)),'.');
    datetick('x','HH:MM','keeplimits','keepticks')
    ylabel('q [g/kg]')
%     ylim([11 18])

    axes('Position',[.69 .15 .15 .15])
    box on
    plot(tcold(i,1:61),dDf(cp_matrix(i,1:61)),'.');
    datetick('x','HH:MM','keeplimits','keepticks')
    ylabel('\deltaD [permil]')
    yyaxis right
    plot(tcold(i,1:61),rrf(cp_matrix(i,1:61)),'.r')
    ylabel('RR [mm/hr]')

%     ylim([-80 -64])