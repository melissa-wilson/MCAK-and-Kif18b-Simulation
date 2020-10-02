function [outcome] = motorAttach(motor_rate,time_step)
% function to determine if motor will attach or dettach based on on/off
% rate and time passed

    % determine probability of motor j attaching/dettaching
    pr_attach = 1-exp(-motor_rate*time_step);
            
    % generate p ~ Unif(0,1)
    p = rand;
    
    if p < pr_attach
        outcome = true;
    else
        outcome = false;
    end
end