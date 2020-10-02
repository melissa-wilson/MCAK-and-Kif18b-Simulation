% Melissa Wilson
% Final Project

clear;
clc;
close all;

% constants and general stuff
number_kif = 200;  % number of Kif18b motors
number_mcak_vector = [4 50 100 200]; % number of MCAK motors
initial_MT_Length_vector = [2000 5000 10000];  % MT Length in nm
duration = 20*60; % duration in sec
time_step = 0.01; % simulation time step in sec
total_steps = round(duration/time_step); % number of time steps
both_depol = 3.620; % nm/s; depolymerization rate of MT when both MCAK and Kif18b at end
motor_off_rate = 0.03; % s^-1

% stuff associated with Kif18b
v_kif_mean = 349.3; % nm/s; mean Kif18b velocity
v_kif_sd = 102.4; % nm/s; sd of Kif18b velocity
kif_on = 1; % on-rate constant for the motor in microM^-1 s^-1

% stuff associated with MCAK
v_mcak_mean = 0; % nm/s; mean MCAK velocity
v_mcak_sd = 380; % nm/s; sd of  MCAK velocity
mcak_on = 1; % on-rate constant for the motor in microM^-1 s^-1
mcak_depol = 7.130; % nm/s; depolymerization rate of MT when only MCAK at end

% empty 3d matrices for storing # of interaction events and MCAK alone
% depolymerization events
association_mat = zeros(3,4,10);
MCAK_alone_mat = zeros(3,4,10);

% repeat whole simulation 10 times
for x = 1:10
    
    % loop through each starting MCAK number
    for m = 1:4

        number_mcak = number_mcak_vector(m); % number of MCAK motors

        % loop through each starting microtubule length
        for n = 1:3

            initial_MT_Length = initial_MT_Length_vector(n);  % MT Length in nm
            mt_conc = 2e-7 * initial_MT_Length; % concentration of microtubule polymer microM
            mcak_on_rate = mcak_on * mt_conc; % on rate of motor in s^-1
            kif_on_rate = kif_on * mt_conc; % on rate of motor in s^-1

            time = zeros(1, total_steps);  % initialize time array
            MT_length_counter = zeros(1,number_kif)+initial_MT_Length; % each MT gets own length
            MT_lengths = zeros(number_kif, total_steps); % length of each MT over time
            time_counter = 0;
            kif_position = zeros(number_kif,total_steps); % position of each kif over time
            kif_position_counter = zeros(1,number_kif); % current position of each kif
            mcak_position = zeros(number_mcak,total_steps); % position of each MCAK over time
            mcak_position_counter = zeros(1,number_mcak); % current position of each MCAK
            mcak_mt = zeros(1,number_mcak); % which microtubule MCAK is bound to (vector of MT indicies)
            mt_mcak = zeros(1,number_kif); % which MCAK on current microtubule (vector of MCAK indicies)
            number_associations = 0; % keeps track of how many times Kif18b and MCAK have associated w each other
            number_MCAK_depol = 0; % keeps track of how many times MCAK alone depolymerized

            % loop through each time step
            for i = 1: total_steps

                time(i) = time_counter; % assign time array
                time_counter = time_counter + time_step;  % increment time counter

                % loop thru each kif18b
                for j = 1:number_kif

                    % if jth kif18b dettached
                    if kif_position_counter(j) == 0

                        % simulate attachment
                        outcome = motorAttach(kif_on_rate,time_step);

                        if outcome == true

                            % motor will attach at a random spot on the microtubule
                            kif_position_counter(j) = rand*MT_length_counter(j);

                        end

                    % if jth kif18b attached
                    else

                        % if no MCAK on same MT
                        if mt_mcak(j) == 0

                            % if at end
                            if kif_position_counter(j) >= MT_length_counter(j)

                                % Kif18b doesn't  depolymerize MTs
                                % simulate dettachment
                                outcome = kifTipDettach(motor_off_rate,time_step);

                                if outcome == true

                                    % kif will dettach
                                    kif_position_counter(j) = 0;

                                end

                            % if not at end
                            else

                                % simulate dettachment
                                outcome = motorAttach(motor_off_rate,time_step);

                                if outcome == true

                                    % kif will dettach
                                    kif_position_counter(j) = 0;

                                else

                                    % kif will not dettach, move kif
                                    kif_vel = v_kif_mean + v_kif_sd*randn;
                                    kif_position_counter(j) = kif_position_counter(j) + kif_vel*time_step;

                                end
                            end

                        % if MCAK on same MT
                        else

                            % if very close to each other
                            if abs(kif_position_counter(j) - mcak_position_counter(mt_mcak(j))) < 10

                                if kif_position_counter(j) ~= mcak_position_counter(mt_mcak(j))

                                    % set equal to each other for simplicity
                                    kif_position_counter(j) = mcak_position_counter(mt_mcak(j));

                                    % increase counter that tracks kif and MCAK
                                    % association
                                    number_associations = number_associations + 1;

                                % if at end
                                elseif kif_position_counter(j) >= MT_length_counter(j)

                                    % simulate kif dettachment
                                    outcome = kifTipDettach(motor_off_rate,time_step);

                                    if outcome == true

                                        % both will dettach
                                        kif_position_counter(j) = 0;
                                        mcak_position_counter(mt_mcak(j)) = 0;
                                        mcak_mt(mt_mcak(j)) = 0;
                                        mt_mcak(j) = 0;

                                    % will not dettach
                                    elseif outcome == false && MT_length_counter(j) < 1

                                        % MT too short, reset length
                                        MT_length_counter(j) = 1;

                                    else

                                        % depolymerize at MCAK + Kif18b rate
                                        MT_length_counter(j) = MT_length_counter(j) - both_depol*time_step;

                                        % both motors move back too
                                        kif_position_counter(j) = kif_position_counter(j) - both_depol*time_step;
                                        mcak_position_counter(mt_mcak(j)) = mcak_position_counter(mt_mcak(j)) - both_depol*time_step;

                                    end

                                % if not at end
                                else

                                    % simulate dettachment
                                    outcome = motorAttach(motor_off_rate,time_step);

                                    if outcome == true

                                        % both will dettach
                                        kif_position_counter(j) = 0;
                                        mcak_position_counter(mt_mcak(j)) = 0;

                                    else

                                        % won't dettach
                                        % both motors move to end at kif rate
                                        velocity = v_kif_mean + v_kif_sd*randn;
                                        kif_position_counter(j) = kif_position_counter(j) + velocity*time_step;
                                        mcak_position_counter(mt_mcak(j)) = mcak_position_counter(mt_mcak(j)) + velocity*time_step;

                                    end
                                end

                            % if not in same location    
                            else

                                % move each motor independently

                                % if kif at end
                                if kif_position_counter(j) >= MT_length_counter(j)

                                    % set kif to be at end
                                    kif_position_counter(j) = MT_length_counter(j);

                                    % simulate kif dettachment
                                    kif_outcome = kifTipDettach(motor_off_rate,time_step);

                                    if kif_outcome == true

                                        % kif dettachs
                                        kif_position_counter(j) = 0;

                                    end

                                % if kif not at end
                                else

                                    % simulate kif dettachment
                                    kif_outcome = motorAttach(motor_off_rate,time_step);

                                    if kif_outcome == true

                                        % kif dettachs
                                        kif_position_counter(j) = 0;

                                    else

                                        % kif moves
                                        kif_vel = v_kif_sd*randn + v_kif_mean;
                                        kif_position_counter(j) = kif_position_counter(j) + kif_vel*time_step;

                                    end
                                end
                            end
                        end
                    end
                end

                % loop thru each MCAK
                for k = 1:number_mcak

                    % if kth MCAK not bound
                    if mcak_mt(k) == 0

                        % simulate attachment
                        outcome = motorAttach(mcak_on_rate,time_step);

                        if outcome == true

                            % will attach

                            % pick which MT to bind to
                            mt = randi([1 number_kif]);

                            % while the chosen MT already has an MCAK
                            while(sum(mt == mcak_mt) ~= 0)

                                % pick another MT
                                 mt = randi([1 number_kif]);

                            end

                            % index corresponds to MCAK, value what MT the kth MCAK is
                            % bound to
                            mcak_mt(k) = mt;

                            % index corresponds to microtubule, value what MCAK is on it
                            mt_mcak(mt) = k;

                            % motor will attach at a random spot on the microtubule
                            mcak_position_counter(k) = rand*MT_length_counter(mt);

                        end

                    % if kth MCAK bound to a MT
                    else

                        % if kth MCAK not associated with a Kif18b
                        if abs(kif_position_counter(mcak_mt(k)) - mcak_position_counter(k)) > 10

                            % simulate mcak dettachment
                            mcak_outcome = motorAttach(motor_off_rate,time_step);

                            if mcak_outcome == true

                                % mcak dettachs
                                mcak_position_counter(k) = 0;
                                mt_mcak(mcak_mt(k)) = 0;
                                mcak_mt(k) = 0;

                            else

                                % mcak doesn't dettach

                                if MT_length_counter(mcak_mt(k)) < 1

                                    % too short, reset length
                                    MT_length_counter(mcak_mt(k)) = 1;

                                % if MCAK past plus end
                                elseif mcak_position_counter(k) > MT_length_counter(mcak_mt(k))

                                    % set MCAK to same length as end
                                    mcak_position_counter(k) = MT_length_counter(mcak_mt(k));

                                    % record MCAK alone depolymerization event
                                    number_MCAK_depol = number_MCAK_depol + 1;

                                % if MCAK at plus end
                                elseif mcak_position_counter(k) == MT_length_counter(mcak_mt(k))

                                    % depolymerize at MCAK alone rate
                                    MT_length_counter(mcak_mt(k)) = MT_length_counter(mcak_mt(k)) - mcak_depol*time_step;

                                    % MCAK moves back too
                                    mcak_position_counter(k) = mcak_position_counter(k) - mcak_depol*time_step;

                                % if MCAK at minus end
                                elseif mcak_position_counter(k) <= 0

                                    % depolymerize at MCAK alone rate
                                    MT_length_counter(mcak_mt(k)) = MT_length_counter(mcak_mt(k)) - mcak_depol*time_step;

                                    % make sure not past end
                                    mcak_position_counter(k) = 0;

                                % if MCAK not on ends  
                                else

                                    % mcak follows a random walk on the MT
                                    mcak_vel = v_mcak_sd*randn + v_mcak_mean;
                                    mcak_position_counter(k) = mcak_position_counter(k) + mcak_vel*time_step;

                                end
                            end
                        end
                    end
                end

                % save positions
                 mcak_position(:,i) = mcak_position_counter;
                 kif_position(:,i) = kif_position_counter;
                 MT_lengths(:,i) = MT_length_counter;
            end
            
            % save info on associations and depolymerization events
            association_mat(n,m,x) = number_associations;
            MCAK_alone_mat(n,m,x) = number_MCAK_depol;
            
            % only plot if first iteration
            if x == 1
                % determine position of subplot
                if n == 1 && m == 1
                    p = 1;
                elseif n == 1 && m == 2
                    p = 2;
                elseif n == 1 && m == 3
                    p = 3;
                elseif n == 1 && m == 4
                    p = 4;
                elseif n == 2 && m == 1
                    p = 5;
                elseif n == 2 && m == 2
                    p = 6;
                elseif n == 2 && m == 3
                    p = 7;
                elseif n == 2 && m == 4
                    p = 8;
                elseif n == 3 && m == 1
                    p = 9;
                elseif n == 3 && m == 2
                    p = 10;
                elseif n == 3 && m == 3
                    p = 11;
                elseif n == 3 && m == 4
                    p = 12;
                end
                
                max_mcak_position = max(mcak_position,[],2);  % calculate max position of each MCAK
                figure(1)
                subplot(3,4,p)
                hist(max_mcak_position)
                xlabel('Maximum MCAK Position (nm)')
                ylabel('Frequency')
                title([num2str(initial_MT_Length), ' nm MT, ', num2str(number_mcak), ' MCAK'])
                xlim([0 initial_MT_Length])

                max_kif_position = max(kif_position,[],2);  % calculate max position of each Kif18b
                figure(2)
                subplot(3,4,p)
                hist(max_kif_position)
                xlabel('Maximum Kif18b Position (nm)')
                ylabel('Frequency')
                title([num2str(initial_MT_Length), ' nm MT, ', num2str(number_mcak), ' MCAK'])
                xlim([0 initial_MT_Length])

                figure(3)
                subplot(3,4,p)
                plot(mcak_position, -time, '.', 'MarkerSize', 6)
                xlabel('MCAK position (nm)')
                ylabel('time (s)')
                title([num2str(initial_MT_Length), ' nm MT, ', num2str(number_mcak), ' MCAK'])
                axis([0 initial_MT_Length -duration 0])

                figure(4)
                subplot(3,4,p)
                plot(kif_position, -time, '.', 'MarkerSize', 6)
                xlabel('Kif18b position (nm)')
                ylabel('time (s)')
                title([num2str(initial_MT_Length), ' nm MT, ', num2str(number_mcak), ' MCAK'])
                axis([0 initial_MT_Length -duration 0])

                figure(5)
                hold all
                plot(time,MT_lengths);
                xlabel('time (s)')
                ylabel('Microtubule length (nm)')
                axis([0 duration 0 max(initial_MT_Length_vector)]);
                hold off
            end
        end
    end
end

% print mean number of association events and MCAK alone depolymerization
% events for each combo of MT length and MCAK number over the 10 iterations
mean(association_mat,3)
mean(MCAK_alone_mat,3)