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
function [t_end,t_end_ind,end_flag] = t_end_w_max_length(delta_T,T_min,max_len,cand_ind,t_min_ind,Taf,time)
%% Identifying t_end: end of cold pool
for k = 1:length(cand_ind)-1
    for ii = t_min_ind(k):1:cand_ind(k+1)
%       if Taf(ii) > T_min(k) + delta_T(k)/exp(1)         % for recovery defined as when ~1/3 of delta_T is reached; 1/3 = 0.3333; 1/exp(1) = 0.3679 = e^(-1); 
        if Taf(ii) > T_min(k) + delta_T(k)*(1-(1/exp(1))) % for recovery defined as when ~2/3 of delta_T is reached; 2/3 = 0.6667; 1-(1/exp(1)) = 0.6321 = 1 - e^(-1);
            if ii < t_min_ind(k)+max_len
                t_end(k) = time(ii);
                t_end_ind(k) = ii;
                end_flag(k) = 1; % 1 = cold pool ends on its own; recovered cold pool
                break
            end
        end
        t_end(k) = time(t_min_ind(k)+max_len);
        t_end_ind(k) = t_min_ind(k)+max_len;
        end_flag(k) = 0; % cold pool ends after max_len amount of time has passed
    end
end
k = length(cand_ind);
for ii = t_min_ind(k):1:length(time)
%   if Taf(ii) > T_min(k) + delta_T(k)/exp(1)         % for recovery defined as when ~1/3 of delta_T is reached; 1/3 = 0.3333; 1/exp(1) = 0.3679 = e^(-1); 
    if Taf(ii) > T_min(k) + delta_T(k)*(1-(1/exp(1))) % for recovery defined as when ~2/3 of delta_T is reached; 2/3 = 0.6667; 1-(1/exp(1)) = 0.6321 = 1 - e^(-1);
        t_end(k) = time(ii);
        t_end_ind(k) = ii;
        end_flag(k) = 1; % 1 = cold pool ends on its own; recovered cold pool
        break
    end
    t_end(k) = time(t_min_ind(k)+max_len);
    t_end_ind(k) = t_min_ind(k)+max_len;
    end_flag(k) = 0; % cold pool ends after max_len amount of time has passed
end
