function [outcome] = kifTipDettach(motor_off_rate,time_step)
% function to determine if Kif18b at tip of microtubule will dettach

    a = 0.43;
    b = 70.8;
    c = 10.0;
    x = motor_off_rate * time_step;

    % calc probability of detaching
    pr_dettach = 1 - a*exp((-x)/b)-(1-a)*exp((-x)/c);

    % generate p ~ Unif(0,1)
    p = rand;
    
    if p < pr_dettach
        outcome = true;
    else
        outcome = false;
    end
end

