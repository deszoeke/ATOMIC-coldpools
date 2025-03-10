%% Computing values for summary table (iso) mixing fractions
% section run twice; once for each cold pool data set
% Using CP#8 and CP#16

% for actual peak d
    [~,peak_ind] = max(dDcold,[],2,'omitnan'); % peak: ~-70 permil
% for "chosen" peak dD (GIVES POSITIVE f_en FRACTIONS)
for k = [8,16]
    dummy = [fen(k,1:61); dDcold(k,1:61); 1:61];
    dummy2 = sortrows(dummy');
    ind = find(dummy2(:,1)>=0);
    peak_ind(k) = dummy2(ind(1),3);
end

% for T_min "valley" [MIGHT NOT BE WORKING CORRECTLY]
    % peak_ind = find(t_cp == t(t_min_ind(i)))

% for individual cold pools
i =  8; % #1 strong cp
i = 16; % #2 strong cp

value2 = fss(i,peak_ind(i))   % peak value
value3 = nanmean(fss_bg(i,:)) % bg value

value2 = fen(i,peak_ind(i))
value3 = nanmean(fen_bg(i,:))

value2 = fee(i,peak_ind(i))
value3 = nanmean(fee_bg(i,:))

% mean for all cold pools (14)
    % for T_min "valley"
    % for k = 2:length(t_min_ind)
    %     min_ind(k) = find(t_cp(i,:) == t(t_min_ind(k)));
    % end
c=0;
for k = [2:11,13:16]
    c = c+1;
    value4(c) = fss(k,peak_ind(k)); % peak value
    value5(c) = fen(k,peak_ind(k)); % peak value
    value6(c) = fee(k,peak_ind(k)); % peak value
end
    value4 = nanmean(value4(2:end))
    value5 = nanmean(value5(2:end))
    value6 = nanmean(value6(2:end))

value3 = nanmean(nanmean(fss_bg([2:11,13:16],:)))     % bg value
value3 = nanmean(nanmean(fen_bg([2:11,13:16],:)))     % bg value
value3 = nanmean(nanmean(fee_bg([2:11,13:16],:)))     % bg value

%% Computing values for summary table (theta-q) mixing fractions
load('conserved_properties_for_isotopes_1min_1.1to1.3km_linear_interp.mat')
% mean of theta-q space mixing fractions
    value1 = nanmean(fss(iso_cold_pool_flag==0)) % non-cold pools
    value2 = nanmean(fss(iso_cold_pool_flag==1)) % cold pools
    
    value1 = nanmean(fen(iso_cold_pool_flag==0)) % non-cold pools
    value2 = nanmean(fen(iso_cold_pool_flag==1)) % cold pools
    
    value1 = nanmean(fdd(iso_cold_pool_flag==0)) % non-cold pools
    value2 = nanmean(fdd(iso_cold_pool_flag==1)) % cold pools
    
    % percentage difference for theta-q space mixing fractions
    per_dif = (value1-value2)/((value1+value2)/2)*100; %  4 for fss
    per_dif = (value1-value2)/((value1+value2)/2)*100; % 21 for fen
    per_dif = (value1-value2)/((value1+value2)/2)*100; % 38 for fdd
    % for ss fractions it is 0.27 for cold pools and 0.26 for non-cold pools.
    % Meanwhile, for dd fractions (0.28 for cold pools and 0.19 for non-cold pools).  
    % for en fractions (0.45 for cold pools and 0.55 for non-cold pools). 

% valley (min temp) of theta-q space mixing fractions (for all cold pools)
    for k = [2:11,13:16]
        dummy = find(t1min >= time(t_max_ind(k)));
%         dummy2 = find(t1min >= time(t_min_ind(k)))
        f_max_ind(k) = dummy(1);
%         f_min_ind(k) = dummy2(1)
    end
    
%     value3 = nanmean(fss(f_min_ind(2:end)))
%     value3 = nanmean(fen(f_min_ind(2:end)))
%     value3 = nanmean(fdd(f_min_ind(2:end)))

% peak (max dD) of theta-q space mixing fractions (for all cold pools)   
    for k = [2:11,13:16]
        dummy = find(t1min >= tcold(k,peak_ind(k)));
        f_peak_ind(k) = dummy(1);
    end
    
    value3 = nanmean(fss(f_peak_ind([2:11,13:16])))
    value3 = nanmean(fen(f_peak_ind([2:11,13:16])))
    value3 = nanmean(fdd(f_peak_ind([2:11,13:16])))

% 60-min (background) theta-q space mixing fractions (for all cold pools)
% fraction data is in 1-min averages
    win = 61; % averaging window;
    nanfss = fss;
    nanfss(iso_cold_pool_flag==1) = NaN;
    nanfdd = fdd;
    nanfdd(iso_cold_pool_flag==1) = NaN;
    nanfen = fen;
    nanfen(iso_cold_pool_flag==1) = NaN;
    
    value4 = [];
    value5 = [];
    value6 = [];
    
    for k = [2:11,13:16]
        dummy = f_max_ind(k)-win:1:f_max_ind(k)-1;
        value4(k) = nanmean(nanfss(dummy));
        value5(k) = nanmean(nanfen(dummy));
        value6(k) = nanmean(nanfdd(dummy));
    end
    
    value4 = nanmean(value4([2:11,13:16]))
    value5 = nanmean(value5([2:11,13:16]))
    value6 = nanmean(value6([2:11,13:16]))

% for individual cold pools
i =  8; % #1 strong cp
i = 16; % #2 strong cp
% % valley (min temp)
% value3 = nanmean(fss(f_min_ind(i)))
% value3 = nanmean(fen(f_min_ind(i)))
% value3 = nanmean(fdd(f_min_ind(i)))

% peak (max dD)
value3 = nanmean(fss(f_peak_ind(i)))
value3 = nanmean(fen(f_peak_ind(i)))
value3 = nanmean(fdd(f_peak_ind(i)))

% background
dummy  = f_max_ind(i)-win:1:f_max_ind(i)-1;
value4 = nanmean(nanfss(dummy))
value5 = nanmean(nanfen(dummy))
value6 = nanmean(nanfdd(dummy))

%% Mean values for each individual cp [don't make any physical sense]
% Only after running section for both cold pools
value1  = 0.4053 % CP 8  FROM nanmean(fss)
value2  = 0.5741 % CP 16 FROM nanmean(fss)

value1  = 0.4563 % CP 8  FROM nanmean(fen)
value2  = 0.2496 % CP 16 FROM nanmean(fen)

value1  = 0.1384 % CP 8  FROM nanmean(fee)
value2  = 0.1763 % CP 16 FROM nanmean(fee)

% percentage difference for theta-q space mixing fractions
per_dif = (value1-value2)/((value1+value2)/2)*100

%% Plotting cold pool progression in theta-q mixing fractions space
hold on;
i =  8; % #1 strong cp
dummy = f_max_ind(i):1:f_max_ind(i)+win-1;
scatter(fss(dummy),fdd(dummy),55,'k')
scatter(fss(dummy),fdd(dummy),55,th_ob(dummy),'filled')
dummy  = f_max_ind(i)-win:1:f_max_ind(i)-1;
plot(fss(dummy),fdd(dummy),'.','Color',[.5 .5 .5]) % bg

i = 16; % #2 strong cp
dummy = f_max_ind(i):1:f_max_ind(i)+win-1;
scatter(fss(dummy),fdd(dummy),55,'k')
scatter(fss(dummy),fdd(dummy),55,th_ob(dummy),'filled')
dummy  = f_max_ind(i)-win:1:f_max_ind(i)-1;
plot(fss(dummy),fdd(dummy),'.','Color',[.5 .5 .5]) % bg
han = colorbar;
han.YLabel.String = ('\theta [K]');
colormap(jet(12))
caxis([294 298])


