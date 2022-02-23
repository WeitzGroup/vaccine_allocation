function f = state_f(parsM, parsS, state, vac_rate)
% vector field f(x) - 6 by 1 vector - S E I R V D 

if parsS.home == 1 % home is A
    contact_vec_home = parsM.cA;
    contact_vec_away = parsM.cB;
    state_vec_home = state.A;
    state_vec_away = state.B;
else
    contact_vec_home = parsM.cB;
    contact_vec_away = parsM.cA;
    state_vec_home = state.B;
    state_vec_away = state.A;
end

F = force_of_infection(parsM,contact_vec_home,...
    contact_vec_away, state_vec_home, state_vec_away);

S = state_vec_home(1);
E = state_vec_home(2);
I = state_vec_home(3);
R = state_vec_home(4);
V = state_vec_home(5);
D = state_vec_home(6);

if S > 0
    vac_rate_true = vac_rate;
else
    vac_rate_true = 0;
end

f = [-F*S - vac_rate_true; F*S - E/parsM.Tinc; E/parsM.Tinc - I/parsM.Tinf;...
    (1 - parsM.mu)*I/parsM.Tinf; vac_rate_true; parsM.mu*I/parsM.Tinf];

end