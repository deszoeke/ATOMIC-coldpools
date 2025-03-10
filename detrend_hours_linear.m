function [dTA] = detrend_hours_linear(TA,time,hours)
    for k = 1:hours*60:length(time)-hours*60
        x = time(k:k+hours*60-1);
        y = TA(k:k+hours*60-1);
        p = polyfit(x,y,1);
        f = polyval(p,x);
        dTA(k:k+hours*60-1) = y - f;
    end
end