function dt = adj_t(star,L, zm,zrf)
% dt = adj_t(star,L, zm,zrf)
% adjust thermal variable from measurement height zm to reference height zm
% Trf = Tm + adj_t(Tstar,L, zm,zrf)
    k=0.4;
    dt = (star./k).*(log(zrf/zm)-psi_T(zrf./L)+psi_T(zm./L));
end