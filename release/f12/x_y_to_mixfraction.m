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

