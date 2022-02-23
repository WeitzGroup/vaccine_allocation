function state_solution = state_solver(parsM, parsC,...
    parsT, parsS, initial_state)
state_solution.A = zeros(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 7); % [t, S, E, I, R, V, D]
state_solution.B = zeros(round((parsT.tf - parsT.t0)/parsT.dt) + 1, 7); % [t, S, E, I, R, V, D]

state_solution.A(1,:) = [parsT.t0, initial_state.A];
state_solution.B(1,:) = [parsT.t0, initial_state.B];

for i = 2:round((parsT.tf - parsT.t0)/parsT.dt) + 1
    state_prev.A = state_solution.A(i - 1,2:end);
    state_prev.B = state_solution.B(i - 1,2:end);
    
    if parsS.VA - state_prev.A(end-1)> 0
        if parsS.VA - state_prev.A(end-1)- ...
                parsS.vaccination_rate_baseline*parsT.dt > 0
            vac_rate.A = parsS.vaccination_rate_baseline; 
        else    % not enough vaccines to sustain nominal rate for dt interval
            vac_rate.A = (parsS.VA-state_prev.A(end-1))/parsT.dt; 
        end    
    else
        vac_rate.A = 0;
    end  
    
    if parsS.VB - state_prev.B(end-1)> 0
        if parsS.VB - state_prev.B(end-1)- ...
                parsS.vaccination_rate_baseline*parsT.dt > 0
            vac_rate.B = parsS.vaccination_rate_baseline; 
        else    % not enough vaccines to sustain nominal rate for dt interval
            vac_rate.B = (parsS.VB-state_prev.B(end-1))/parsT.dt; 
        end
    else
        vac_rate.B = 0;
    end 
    
    parsS.home = 1;
    der.A = state_f(parsM, parsS, state_prev, vac_rate.A);
    
    parsS.home = 0;
    der.B = state_f(parsM, parsS, state_prev, vac_rate.B);
    
    state_solution.A(i,2:end) = state_prev.A + parsT.dt.*der.A';
    state_solution.B(i,2:end) = state_prev.B + parsT.dt.*der.B';
    
    % if any state is below 0, set to 0
    
    
    state_solution.A(i,1) = parsT.t0 + (i - 1)*parsT.dt;
    state_solution.B(i,1) = state_solution.A(i,1);
    
end

end