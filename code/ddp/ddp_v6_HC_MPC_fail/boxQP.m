function [x_opt, clamped] = boxQP(h, g, x_lo, x_up)
    clamped = false;
    
    if h > 0
        x_star = -g / h;
        if x_lo <= x_star && x_star <= x_up
            x_opt = x_star;
        elseif x_star < x_lo
            x_opt = x_lo;
            clamped = true;
        else
            x_opt = x_up;
            clamped = true;
        end
    elseif h < 0
        f_x_lo = 0.5 * h * x_lo^2 + g * x_lo;
        f_x_up = 0.5 * h * x_up^2 + g * x_up;
        if f_x_lo <= f_x_up
            x_opt = x_lo;
            clamped = true;
        else
            x_opt = x_up;
            clamped = true;
        end
    else
        if g > 0
            x_opt = x_lo;
            clamped = true;
        elseif g < 0
            x_opt = x_up;
            clamped = true;
        else
            x_opt = (x_lo + x_up) / 2;
            clamped = false;
        end
    end
end