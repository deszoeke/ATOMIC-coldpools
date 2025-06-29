function [fen, fss, fdd] = th_q_to_mixfraction(thml,qml, then,qen, thss,qss, thdd,qdd)
% [fen, fss, fdd] = th_q_to_mixfraction(thml,qml, then,qen, thss,qss, thdd,qdd)
% Diagnoses BL air as a 3-part mixture fen + fss + fdd = 1
% of entrained air, air in equilibrium with the sea surface,
% and saturated downdraft air.
% The mixing fractions of these 3 end members are inverted algebraically
% from observed thml, qml as in de Szoeke (2018).
%
% (c) 2020-08-19 Simon P. de Szoeke 

qen(qen>30)=NaN;
% solves algebraically for mixing fractions fen,fss,fdd at each time
kay=(thss-thdd)./(qdd-qss);
beta=thdd+kay.*qdd;
gamma=then+kay.*qen-beta;
fen=(thml+kay.*qml-beta)./gamma;
fss=(qml-qdd-(qen-qdd).*fen)./(qss-qdd);
fdd=1-fen-fss;

