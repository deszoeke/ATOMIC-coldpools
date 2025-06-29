% Valid for 3-member mixing %
function [fw1,fw2,fw3,q] = mass_mixing(ind,fev,fss,fen,q_ev,q_surf,q_en)    
    fa1 = fev(ind);
    fa2 = fss(ind);
    fa3 = fen(ind);
    q1  = q_ev(ind);
    q2  = q_surf(ind);
    q3  = q_en(ind);
    q   = fa1.*q1 + fa2.*q2 + fa3.*q3;
    fw1 = (fa1.*q1)./q;
    fw2 = (fa2.*q2)./q;
    fw3 = (fa3.*q3)./q;
