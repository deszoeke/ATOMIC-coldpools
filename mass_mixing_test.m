% Valid for 2-member mixing %
function [fw1,fw2] = mass_mixing_test(ind,fev,fss,q_ev,q_surf)    
    fa1 = fev(ind);
    fa2 = fss(ind);
    q1 = q_ev(ind);
    q2 = q_surf(ind);
    q = fa1.*q1 + fa2.*q2;
    fw1 = (fa1.*q1)./q;
    fw2 = (fa2.*q2)./q;
