function [t_max,t_min,t_max_ind,t_min_ind,t_end,t_end_ind,end_flag,cp_ind,delta_T,T_max,T_min,Taf] = cold_pool_detection_algorithm(t1min,Ta)

    Taf = movmean(Ta,11,'omitnan'); % Ta filtered => 11-min running average
    del_T = Taf(2:end)-Taf(1:end-1);
    cand = find(del_T<-0.05); % candidates positions
    %% Separating candidates [118 possible cold pools]
    c=1;
    cand_ind(c) = cand(1);
    for i = 2:length(cand)-1
        if cand(i)-cand(i-1)>1 && cand(i+1)-cand(i)==1
            c = c+1;
            cand_ind(c) = cand(i);
        end
    end
    % Time passed between possible cold pools
    t_pass = t1min(cand_ind(2:end)) - t1min(cand_ind(1:end-1)); % in days
    t_pass = t_pass*24*60; % in minutes
    t_max_offset = 20; % in minutes
    find(t_pass<t_max_offset);
    %% Identifying t_max: time of the cold pool front onset
    c=1;
    for k = 1:length(cand_ind)
        for ii = cand_ind(k):-1:1
            if del_T(ii)>0 && ((t1min(cand_ind(k))-t1min(ii+1))*24*60)<t_max_offset
                t_max(c) = t1min(ii+1);
                t_max_ind(c) = ii+1;
                c = c+1;
                break
            elseif del_T(ii)>0 && ((t1min(cand_ind(k))-t1min(ii+1))*24*60)>=t_max_offset
                t_max(c) = t1min(cand_ind(k)-t_max_offset);
                t_max_ind(c) = (cand_ind(k))-t_max_offset;
                c = c+1;
                break
            end
        end
    end
    % Time passed between t_max times
    % t_p = t1min(cand_ind) - t_max; % in days
    % t_p = t_p*24*60; % in min
    t_p = cand_ind - t_max_ind; % in min
    find(t_p>20);
    %% Identifying t_min: time of end of cold pool front
    c=1;
    for k = 1:length(cand_ind)
        for ii = cand_ind(k):1:length(del_T)
            if del_T(ii)>0
                t_min(c) = t1min(ii);
                t_min_ind(c) = ii;
                c = c+1;
                break
            end
        end
    end
    t_min_offset = 20; % within 20min of the previous minimum
    %% Identifying subsequent cold pools that can be merged
    % Add a flag for merged with previous!!!
    % 17 merged cold pools identified; length(merge_flag(merge_flag==0)) = 17
    merge_flag(1) = 1;
    for ii = 2:length(cand_ind)
        % checking whether subsequent cold pools are less than 20 min apart
        if cand_ind(ii) - t_min_ind(ii-1) <= t_min_offset
            chunk = find(Ta(t_min_ind(ii-1):cand_ind(ii))-Ta(t_min_ind(ii-1)) > 0.5); % should I be using Taf instead?
            if isempty(chunk) == 1
                % subsequent cold pools can be merged
                merge_flag(ii) = 0; % 1 = stand alone cold pool
                                 % 0 = can be merged with previous
            else
                merge_flag(ii) = 1;
            end
        else
            merge_flag(ii) = 1;
        end
    end
    %% Merging subsequent cold pools identified in previous section
    t_max_ind = t_max_ind(merge_flag==1);
    t_max = t_max(merge_flag==1);
    t_min_ind = t_min_ind(merge_flag==1);
    t_min = t_min(merge_flag==1);
    %% Identifying t_end: end of cold pool
    % 22 cold pools end with next cold pool!!!
    % find(end_flag==0)
    % ans =
    %      2     3     14    22    30    40    45    46    47    50    52    
    %      57    66*   69*   73*   74*   88    91    92    93    99    101
    % *Only 4 cold pools did not recover during the iso data time frame!!!
    cp_ind = cand_ind(merge_flag==1);
    delta_T = Taf(t_max_ind) - Taf(t_min_ind);
    T_min = Taf(t_min_ind);
    T_max = Taf(t_max_ind);

    for k = 1:length(cp_ind)-1
        for ii = t_min_ind(k):1:cp_ind(k+1)
    %         if Taf(ii) > T_min(k) + delta_T(k)/exp(1) % for recovery defined as when ~1/3 of delta_T is reached; 1/3 = 0.3333; 1/exp(1) = 0.3679 = e^(-1); 
            if Taf(ii) > T_min(k) + delta_T(k)*(1-(1/exp(1))) % for recovery defined as when ~2/3 of delta_T is reached; 2/3 = 0.6667; 1-(1/exp(1)) = 0.6321 = 1 - e^(-1);
                t_end(k) = t1min(ii);
                t_end_ind(k) = ii;
                end_flag(k) = 1; % 1 = cold pool ends on its own; recovered cold pool
                if t_end_ind(k) - t_max_ind(k) > 120
                    t_end_ind(k) = t_max_ind(k) + 120;
                    t_end(k) = t1min(t_end_ind(k));
                end
                break
            end
            t_end(k) = t1min(ii-1);
            t_end_ind(k) = ii-1;
            end_flag(k) = 0; % cold pool ends at the onset of the next cold pool
            if t_end_ind(k) - t_max_ind(k) > 120
                t_end_ind(k) = t_max_ind(k) + 120;
                t_end(k) = t1min(t_end_ind(k));
            end
        end
    end

    k = length(cp_ind);
    for ii = t_min_ind(k):1:length(t1min)
    %     if Taf(ii) > T_min(k) + delta_T(k)/exp(1) % for recovery defined as when ~1/3 of delta_T is reached; 1/3 = 0.3333; 1/exp(1) = 0.3679 = e^(-1); 
        if Taf(ii) > T_min(k) + delta_T(k)*(1-(1/exp(1))) % for recovery defined as when ~2/3 of delta_T is reached; 2/3 = 0.6667; 1-(1/exp(1)) = 0.6321 = 1 - e^(-1);
            t_end(k) = t1min(ii);
            t_end_ind(k) = ii;
            end_flag(k) = 1; % 1 = cold pool ends on its own; recovered cold pool
            if t_end_ind(k) - t_max_ind(k) > 120
                t_end_ind(k) = t_max_ind(k) + 120;
                t_end(k) = t1min(t_end_ind(k));
            end
            break
        end
        t_end(k) = t1min(ii-1);
        t_end_ind(k) = ii-1;
        end_flag(k) = 0; % cold pool ends at the onset of the next cold pool
        if t_end_ind(k) - t_max_ind(k) > 120
            t_end_ind(k) = t_max_ind(k) + 120;
            t_end(k) = t1min(t_end_ind(k));
        end
    end
    % Time passed between t_max and t_end times [total cold pool times]
    total_p = t_end_ind - t_max_ind; % in min
    t_end_offset = 1; % 5 minutes
    %% REPEAT: Identifying subsequent cold pools that can be merged
%     % Add a flag for merged with previous!!!
%     % 17 merged cold pools identified; length(merge_flag(merge_flag==0)) = 17
%     merge_flag2(1) = 1;
%     for ii = 2:length(t_max_ind)
%         % checking whether subsequent cold pools are less than 1 min apart
%         if t_max_ind(ii) - t_end_ind(ii-1) <= t_end_offset
%                 % subsequent cold pools can be merged
%                 merge_flag2(ii) = 0; % 1 = stand alone cold pool
%                                      % 0 = can be merged with previous
%         else
%             merge_flag2(ii) = 1;
%         end
%     end
    %% REPEAT: Merging subsequent cold pools identified in previous section
%     t_max_ind = t_max_ind(merge_flag2==1);
%     t_max = t_max(merge_flag2==1);
%     t_min_ind = t_min_ind(merge_flag2==1);
%     t_min = t_min(merge_flag2==1);
%     t_end_ind = t_end_ind(merge_flag2==1);
%     t_end = t_end(merge_flag2==1);
%     T_min = T_min(merge_flag2==1);
%     T_max = T_max(merge_flag2==1);
%     total_p = t_end_ind - t_max_ind; % in min
%     cp_ind = cp_ind(merge_flag2==1);
    %%
    cold_pool_flag_1min = zeros(size(t1min));
    for k = 1:length(t_max)
        cold_pool_flag_1min(t_max_ind(k):t_end_ind(k)) = 1;
    end
    % for plotting purposes
    % variables for shading
    yvalues = [22 22 28 28];
    figure;
    hold on;
    for k = 1:length(cp_ind)
    xvalues = [t1min(t_max_ind(k)) t1min(t_end_ind(k)) t1min(t_end_ind(k)) t1min(t_max_ind(k))];
    patch(xvalues, yvalues, [0.8 0.8 0.8]);%,'LineStyle','none')
    clearvars xvalues
    end
    plot(t1min,Taf,'-k','Color','k')
    yyaxis right
    plot(t1min,cold_pool_flag_1min,'-r','Color','r')
    datetick('x','mm/dd','keeplimits','keepticks')
    title('cold pool flag'); 
    % ylim([-1 2])
end