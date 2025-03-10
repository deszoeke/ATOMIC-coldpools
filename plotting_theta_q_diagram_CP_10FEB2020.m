%% Plotting theta-q diagram for CP 10 FEB 2020 ONLY
% code runs with output from "q_dD_FRAMEWORK_COLD_POOLS.m" script
% requires running "q_dD_FRAMEWORK_COLD_POOLS.m" script first
% and loading 'conserved_properties_for_isotopes_1min_1.1to1.3km_nearest_interp.mat'
load 'conserved_properties_for_isotopes_1min_1.1to1.3km_nearest_interp.mat' t1min th_ob q_ob th_surf q_surf th_ent q_ent th_dd q_dd
%% Identifying cold pool times in larger matrix
for k = [2:11,13:16]
    dummy  = find(t1min >= time(t_max_ind(k)));
    dummy2 = find(t1min >= time(t_min_ind(k))); % min temp of theta-q space (for each cold pools)
    dummy3 = find(t1min >= time(t_end_ind(k)));

    max_ind(k) = dummy(1);
    min_ind(k) = dummy2(1);
    end_ind(k) = dummy3(1);
end
%% For individual cold pool CP2
i = 16; % #2 strong cp
% entire CP2 data in theta-q
th_dd(max_ind(i):end_ind(i));
th_ent(max_ind(i):end_ind(i));
th_surf(max_ind(i):end_ind(i));
th_ob(max_ind(i):end_ind(i));

q_dd(max_ind(i):end_ind(i));
q_ent(max_ind(i):end_ind(i));
q_surf(max_ind(i):end_ind(i));
q_ob(max_ind(i):end_ind(i));

% to highlight time of min temp
th_dd(min_ind(i));
th_ent(min_ind(i));
th_surf(min_ind(i));
th_ob(min_ind(i));

q_dd(min_ind(i));
q_ent(min_ind(i));
q_surf(min_ind(i));
q_ob(min_ind(i));

%% Background data in theta-q preceding CP2
% 60-min (background) theta-q space
% theta-q data is in 1-min averages
i = 16; % #2 strong cp
win = 61; % averaging window;
dummy  = max_ind(i)-win:1:max_ind(i)-1;
th_bg = th_ob(dummy);
q_bg  = q_ob(dummy);

%% Plotting cold pool progression in theta-q space
yellow = [0.9290 0.6940 0.1250];
orange = [0.8500 0.3250 0.0980];
blue = [0 0.4470 0.7410];
i = 16; % #2 strong cp
% dummy = max_ind(i):1:max_ind(i)+win-1; % represents first 60 minutes of cp
% dummy = max_ind(i):1:end_ind(i); % represents entire cold pool
% figure; 
hold on;
dummy = min_ind(i); % represents time of min temp
plot(th_ob(dummy),q_ob(dummy),'ok')
scatter(th_surf(dummy),q_surf(dummy),55,yellow,'filled')
scatter(th_dd(dummy),q_dd(dummy),55,blue,'filled')
scatter(th_ent(dummy),q_ent(dummy),55,orange,'filled')
plot(th_bg,q_bg,'o','Color',[.5 .5 .5]) % bg
% scatter(th_bg,q_bg,25,[.5 .5 .5])

%% Plotting cold pool progression in theta-q mixing fractions space
% find it in 'table_values_mixing_fractions.m' script