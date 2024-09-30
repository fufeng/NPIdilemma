function Fv=SIR_SD_FC(t,Y,beta,beta_fc,gamma,omega,kappa,C_SD,C_I,C_FC)
%I_star = -ln(1 - C_d/C_i)/beta;
Fv(1,1)= -beta*Y(4)*Y(1)*Y(2) - beta_fc*(1-Y(4)-Y(3))*Y(1)*Y(2); % Susceptible
Fv(2,1)= beta*Y(4)*Y(1)*Y(2) + beta_fc*(1-Y(4)-Y(3))*Y(1)*Y(2) -gamma*Y(2); % Infected
Fv(3,1)= omega*Y(3)*Y(4)*tanh(kappa/2*(-C_SD+C_I*(1-exp(-beta*Y(2))))) ...
    + omega*Y(3)*(1-Y(4)-Y(3))*tanh(kappa/2*(-C_SD+C_FC+C_I*(1-exp(-beta_fc*Y(2))))); %sd
Fv(4,1)= omega*Y(3)*Y(4)*tanh(kappa/2*(C_SD-C_I*(1-exp(-beta*Y(2))))) ...
    + omega*Y(4)*(1-Y(4)-Y(3))*tanh(kappa/2*(C_FC+C_I*(exp(-beta*Y(2))-exp(-beta_fc*Y(2))))); %nsd
