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