function [] = lk_mc()
% Schlogl model:
% A  --> X
% K1 = k1/k2;
% B --> X
% Chain of reactions:
% A0 -> A1 -> A2 ->A3 ->A4 +2X -> 3X
% X -> B4 ->B3 -> B2 -> B1 -> B0
% A0 and B0 are the resevoirs

% Variables:
V = 1.00;
R = 0.00198858775; 
T = 298.15;
RT = R*T;

% kinetic parameters
k_boundary_rate = 1.00;
scale1 = 1;
scale2 = 1;
k1 = 1*3.0/V * scale1* k_boundary_rate;
k2 = 1*0.6/V * scale1* k_boundary_rate;
k3 = 1*0.25/V * scale2* k_boundary_rate;
k4 = k3
%k4=2.95;
K1 = k1/k2
K2 = k3/k4
k_boundary_rate = k_boundary_rate/V;

% Overall free energy change:
%DG = -3.85 %Gives a A0 =~ 10
%DG = -5.00
DG = -2.00
%DG = -3.85
%DG = -0.00
L = exp(-DG/RT)
A2B = L * K2/K1;
B0 = 100
A0 = round(A2B * B0)
pause;
X0 = 100;
X = X0;

L = K1/K2 * A0/B0
L1 = K1*A0/X
L2 = (1/K2) * X/B0
L1L2 = L1*L2
G_target = -RT* log(K1/K2 * A0/B0)

A1 = A0;
A2 = A0;
A3 = A0;
A4 = A0;
B1 = B0;
B2 = B0;
B3 = B0;
B4 = B0;
nsteps = 10000;
nsteps = 1000;
choices = rand(nsteps,1);

x_states = [];
a_states = [];
b_states = [];
Gstates = [];
deltaGstates = [];
xflux = [];
newtime = 0;
isteps = 0;
nreactions = 10;
L_thermo = zeros(nreactions, nsteps);
L_kinetic = zeros(nreactions, nsteps);

%Thermodynamic states - equilibrate the system:
for i=1:1:nsteps
   isteps = isteps + 1;

%Calculate the free energy changes of the reactions:
[f1,r1,f2,r2,f3,r3,f4,r4,f5,r5,f6,r6,fb3,rb3,fb4,rb4,fb5,rb5,fb6,rb6] = free_energy_changes(X, A0, A1, A2, A3, A4, B0, B1, B2, B3, B4, K1, K2);
   fall = f1+f2+f3+f4+f5+f6+fb3+fb4+fb5+fb6;
   rall = r1+r2+r3+r4+r5+r6+rb3+rb4+rb5+rb6;
   vall = 1 + rall + fall;

   %Check to see if any of the reactants went to zero:
   if((X < 0) || (B4 < 0))
      pause_sim = 1;
   else
      pause_sim = 0;
   end;
   % Random number:
   choice = choices(i);
if(X >= 1)
   % Choose a reaction:
   if (choice <= f1/vall)
      X = X +1;
      A4 = A4-1;
   elseif (choice <= (f1+f2)/vall)
      X = X +1;
      B4 = B4 -1;
   elseif (choice <= (f1+f2+f3)/vall)
      A1 = A1 +1;
   elseif (choice <= (f1+f2+f3+f4)/vall)
      A2 = A2 +1;
      A1 = A1 -1;
   elseif (choice <= (f1+f2+f3+f4+f5)/vall)
      A3 = A3 +1;
      A2 = A2 -1;
   elseif (choice <= (f1+f2+f3+f4+f5+f6)/vall)
      A4 = A4 +1;
      A3 = A3 -1;
   elseif (choice <= (f1+f2+f3+f4+f5+f6+fb3)/vall)
      B1 = B1 -1;
      %B0 is a constant
   elseif (choice <= (f1+f2+f3+f4+f5+f6+fb3+fb4)/vall)
      B1 = B1 + 1;
      B2 = B2 - 1;
   elseif (choice <= (f1+f2+f3+f4+f5+f6+fb3+fb4+fb5)/vall)
      B2 = B2 + 1;
      B3 = B3 - 1;
   elseif (choice <= (f1+f2+f3+f4+f5+f6+fb3+fb4+fb5+fb6)/vall)
      B3 = B3 + 1;
      B4 = B4 - 1;
%========================================
   elseif (choice <= (fall+r1)/vall)
      X = X -1;
      A4 = A4 +1;
   elseif (choice <= (fall+r1+r2)/vall)
      X = X -1;
      B4 = B4 +1;
   elseif (choice <= (r1+r2+r3)/vall)
      A1 = A1 -1;
   elseif (choice <= (r1+r2+r3+r4)/vall)
      A2 = A2 -1;
      A1 = A1 +1;
   elseif (choice <= (r1+r2+r3+r4+r5)/vall)
      A3 = A3 -1;
      A2 = A2 +1;
   elseif (choice <= (r1+r2+r3+r4+r5+r6)/vall)
      A4 = A4 -1;
      A3 = A3 +1;
   elseif (choice <= (r1+r2+r3+r4+r5+r6+rb3)/vall)
      B1 = B1 +1;
      %B0 is a constant
   elseif (choice <= (r1+r2+r3+r4+r5+r6+rb3+rb4)/vall)
      B1 = B1 - 1;
      B2 = B2 + 1;
   elseif (choice <= (r1+r2+r3+r4+r5+r6+rb3+rb4+rb5)/vall)
      B2 = B2 - 1;
      B3 = B3 + 1;
   elseif (choice <= (r1+r2+r3+r4+r5+r6+rb3+rb4+rb5+rb6)/vall)
      B3 = B3 - 1;
      B4 = B4 + 1;

   else
      % do nothing - the most likely state is the current state;
   end
   if (pause_sim ==1)
      X = X
      B4 = B4
      pause;
   end
   x_states = [x_states X];
   a_states = [a_states A2];
   b_states = [b_states B2];
else
   % By deduction, if the no X is present, the only possible reaction is that of B --> X:
   % In this case, the free energy of this state is simply -RT ln Keq * A/B where A/B = constant.
   %if (X < 1)
     % This is equivalent to assuming that an B--> X reaction must occur and is the only possibility
   if (pause_sim ==1)
      X = X
      B4 = B4
   end
     X = 1;
     %nsteps = isteps;
     %break;
   x_states = [x_states X];
   end;

end;

x_states = [];
a_states = [];
b_states = [];
Gstates = [];
deltaGstates = [];
xflux = [];
newtime = 0;
isteps = 0;

figure;
hold on;
%Thermodynamic states - Data collection:
for i=1:1:nsteps
   isteps = isteps + 1;
   % likelihood of forward reaction 
   % Free energy *change* requires incrementing products from current population
   % Current Free energy requires uses products from current population
%   if(i < 20000) X = 100; end;
      
%Calculate the free energy changes of the reactions:
[f1,r1,f2,r2,f3,r3,f4,r4,f5,r5,f6,r6,fb3,rb3,fb4,rb4,fb5,rb5,fb6,rb6] = free_energy_changes(X, A0, A1, A2, A3, A4, B0, B1, B2, B3, B4, K1, K2);
   fall = f1+f2+f3+f4+f5+f6+fb3+fb4+fb5+fb6;
   rall = r1+r2+r3+r4+r5+r6+rb3+rb4+rb5+rb6;
   vall = 1 + rall + fall;
   G1 = f1; 
   G2 = f2; 
   G3 = r1; 
   G4 = r2;

   if (f1 == 0) 
      G1 = 1; 
   end
   if (f2 == 0) 
      G2 = 1; 
   end
   if (r1 == 0) 
      G3 = 1; 
   end
   if (r2 == 0) 
      G4 = 1; 
   end
   if((X < 0) || (B4 < 0))
      pause_sim = 1;
   else
      pause_sim = 0;
   end;
 
   % Random number
   choice = choices(i);
if(X >= 1)
   %Choose a reaction:
   if (choice <= f1/vall)
      X = X +1;
      A4 = A4-1;
   elseif (choice <= (f1+f2)/vall)
      X = X +1;
      B4 = B4 -1;
   elseif (choice <= (f1+f2+f3)/vall)
      A1 = A1 +1;
   elseif (choice <= (f1+f2+f3+f4)/vall)
      A2 = A2 +1;
      A1 = A1 -1;
   elseif (choice <= (f1+f2+f3+f4+f5)/vall)
      A3 = A3 +1;
      A2 = A2 -1;
   elseif (choice <= (f1+f2+f3+f4+f5+f6)/vall)
      A4 = A4 +1;
      A3 = A3 -1;
   elseif (choice <= (f1+f2+f3+f4+f5+f6+fb3)/vall)
      B1 = B1 -1;
      %B0 is a constant
   elseif (choice <= (f1+f2+f3+f4+f5+f6+fb3+fb4)/vall)
      B1 = B1 + 1;
      B2 = B2 - 1;
   elseif (choice <= (f1+f2+f3+f4+f5+f6+fb3+fb4+fb5)/vall)
      B2 = B2 + 1;
      B3 = B3 - 1;
   elseif (choice <= (f1+f2+f3+f4+f5+f6+fb3+fb4+fb5+fb6)/vall)
      B3 = B3 + 1;
      B4 = B4 - 1;
%========================================
   elseif (choice <= (fall+r1)/vall)
      X = X -1;
      A4 = A4 +1;
   elseif (choice <= (fall+r1+r2)/vall)
      X = X -1;
      B4 = B4 +1;
   elseif (choice <= (r1+r2+r3)/vall)
      A1 = A1 -1;
   elseif (choice <= (r1+r2+r3+r4)/vall)
      A2 = A2 -1;
      A1 = A1 +1;
   elseif (choice <= (r1+r2+r3+r4+r5)/vall)
      A3 = A3 -1;
      A2 = A2 +1;
   elseif (choice <= (r1+r2+r3+r4+r5+r6)/vall)
      A4 = A4 -1;
      A3 = A3 +1;
   elseif (choice <= (r1+r2+r3+r4+r5+r6+rb3)/vall)
      B1 = B1 +1;
      %B0 is a constant
   elseif (choice <= (r1+r2+r3+r4+r5+r6+rb3+rb4)/vall)
      B1 = B1 - 1;
      B2 = B2 + 1;
   elseif (choice <= (r1+r2+r3+r4+r5+r6+rb3+rb4+rb5)/vall)
      B2 = B2 - 1;
      B3 = B3 + 1;
   elseif (choice <= (r1+r2+r3+r4+r5+r6+rb3+rb4+rb5+rb6)/vall)
      B3 = B3 - 1;
      B4 = B4 + 1;

   else
      % do nothing - the most likely state is the current state;
   end
   if (pause_sim ==1)
      X = X
      B4 = B4
      pause;
   end
   x_states = [x_states X];
   a_states = [a_states A2];
   b_states = [b_states B2];
else
   % By deduction, if the no X is present, the only possible reaction is that of B --> X:
   % In this case, the free energy of this state is simply -RT ln Keq * A/B where A/B = constant.
   %if (X < 1)
     % This is equivalent to assuming that an B--> X reaction must occur and is the only possibility
   if (pause_sim ==1)
      X = X
      B4 = B4
   end
     X = 1;
     %nsteps = isteps;
     %break;
   X = X
   x_states = [x_states X];
end;
   %fitness = -RT*(log(G1)+log(G2))
   [f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1, B2, B3, B4, K1, K2);
   DG_forward = -RT*log(f1*r2*f3*f4*f5*f6*fb3*fb4*fb5*fb6);
   Gstates = [Gstates DG_forward];
   
   L_thermo(1,i) = 1;
   L_thermo(2,i) = f3;
   L_thermo(3,i) = f4;
   L_thermo(4,i) = f5;
   L_thermo(5,i) = f6;
   L_thermo(6,i) = f1;
   L_thermo(7,i) = r2;
   L_thermo(8,i) = fb6;
   L_thermo(9,i) = fb5;
   L_thermo(10,i) = fb4;
   L_thermo(11,i) = fb3;
      
   %plot([0:10], cumsum(-RT*log(L_thermo(:,i))));
   %plot([0:10], (-RT*log(L_thermo(:,i))));
   plot([0:10], [A0 A1 A2 A3 A4 X B4 B3 B2 B1 B0], 'b-');
   title('Product Counts');
   %pause;

   %Constants:
%   if(acetyl_coA== 0) acetyl_coA = 1; end

   %G = G_equil + (log(x/a) + log(y/x) + log (1/y));  
   %deltaG = (log(x/a) + log(y/x) + log (1/y));  
   %Gstates = [Gstates G];
   %deltaGstates = [deltaGstates deltaG];
end
   
isteps = isteps
nx_states = size(x_states)
mean_x = mean(x_states)
mean_a = mean(a_states)
mean_b = mean(b_states)
arith_avg_g = mean(Gstates)
   
figure
plot([1:nsteps] , L_thermo(1,:),'b-', ...
     [1:nsteps] , L_thermo(2,:),'g-', ...
     [1:nsteps] , L_thermo(3,:),'r-', ...
     [1:nsteps] , L_thermo(4,:),'c-', ...
     [1:nsteps] , L_thermo(5,:),'m-', ...
     [1:nsteps] , L_thermo(6,:),'y-', ...
     [1:nsteps] , L_thermo(7,:),'k-', ...
     [1:nsteps] , L_thermo(8,:),'b--', ...
     [1:nsteps] , L_thermo(9,:),'g--', ...
     [1:nsteps] , L_thermo(10,:),'r--');
title('Likelihood of Rxn per Simulation Step');

figure
thermo_L = [mean(L_thermo(1,:)),mean(L_thermo(2,:)),mean(L_thermo(3,:)),mean(L_thermo(4,:)),mean(L_thermo(5,:)),mean(L_thermo(6,:)),mean(L_thermo(7,:)),mean(L_thermo(8,:)),mean(L_thermo(9,:)),mean(L_thermo(10,:))];
scale = 1/log(2);
plot([0:10],scale*log(L_thermo),'b:');
title('Likelihood of Rxn');

%G_thermo = -RT* cumsum(log(thermo_L))
%G_kinetic = -RT* cumsum(log(kinetic_L))
figure
G_thermo = -RT* cumsum(log(L_thermo));
plot([0:10],G_thermo,'b');
title('Cumulative Free Energy Change');
figure
G_thermo = -RT* (log(L_thermo));
plot([0:10],G_thermo,'b');
title('Free Energy Change');


function [f1,r1,f2,r2,f3,r3,f4,r4,f5,r5,f6,r6,fb3,rb3,fb4,rb4,fb5,rb5,fb6,rb6] = free_energy_changes(X, A0, A1, A2, A3, A4, B0, B1, B2, B3, B4, K1, K2)

   [f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1, B2, B3, B4, K1, K2);
   current_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);

   % A4 + 2X --> 3X
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X+1, A0, A1, A2, A3, A4-1, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   f1 = exp(sum(next_log_likelihood - current_log_likelihood));
   % B --> X
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X+1, A0, A1, A2, A3, A4, B0, B1, B2, B3, B4-1, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   f2 = exp(sum(next_log_likelihood - current_log_likelihood));
   % likelihood of reverse reaction
   % 3X --> A + 2X
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X-1, A0, A1, A2, A3, A4+1, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   r1 = exp(sum(next_log_likelihood - current_log_likelihood));
        
   % X --> B
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X-1, A0, A1, A2, A3, A4, B0, B1, B2, B3, B4+1, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   r2 = exp(sum(next_log_likelihood - current_log_likelihood));
   %A0 -> A1
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1+1, A2, A3, A4, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   f3 = exp(sum(next_log_likelihood - current_log_likelihood));
   %A1 -> A0
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1-1, A2, A3, A4, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   r3 = exp(sum(next_log_likelihood - current_log_likelihood));
        
   % A1 -> A2
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1-1, A2+1, A3, A4, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   f4 = exp(sum(next_log_likelihood - current_log_likelihood));
   % A2 -> A1
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1+1, A2-1, A3, A4, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   r4 = exp(sum(next_log_likelihood - current_log_likelihood));
   % A2 -> A3
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2-1, A3+1, A4, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   f5 = exp(sum(next_log_likelihood - current_log_likelihood));
   % A3 -> A2
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2+1, A3-1, A4, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   r5 = exp(sum(next_log_likelihood - current_log_likelihood));
   % A3 -> A4
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3-1, A4+1, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   f6 = exp(sum(next_log_likelihood - current_log_likelihood));
   % A4 -> A3
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3+1, A4-1, B0, B1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   r6 = exp(sum(next_log_likelihood - current_log_likelihood));

   % B1 -> B0 (B0 constant)
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1-1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   fb3 =  exp(sum(next_log_likelihood - current_log_likelihood));
   % B0 -> B1
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1+1, B2, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   rb3 =  exp(sum(next_log_likelihood - current_log_likelihood));
   % B2 -> B1
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1+1, B2-1, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   fb4 =  exp(sum(next_log_likelihood - current_log_likelihood));
   % B1 -> B2
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1-1, B2+1, B3, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   rb4 =  exp(sum(next_log_likelihood - current_log_likelihood));
   % B3 -> B2
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1, B2+1, B3-1, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   fb5 =  exp(sum(next_log_likelihood - current_log_likelihood));
   % B2 -> B3
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1, B2-1, B3+1, B4, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   rb5 =  exp(sum(next_log_likelihood - current_log_likelihood));
   % B4 -> B3
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1, B2, B3+1, B4-1, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   fb6 =  exp(sum(next_log_likelihood - current_log_likelihood));
   % B3 -> B4
   [f1 f3 f4 f5 f6 r2 fb3 fb4 fb5 fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1, B2, B3-1, B4+1, K1, K2);
   next_log_likelihood = log([f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6]);
   rb6 =  exp(sum(next_log_likelihood - current_log_likelihood));
return;

function [f1,f3,f4,f5,f6,r2,fb3,fb4,fb5,fb6] = current_free_energy(X, A0, A1, A2, A3, A4, B0, B1, B2, B3, B4, K1, K2)

   counts = [X A0 A1 A2 A3 A4 B0 B1 B2 B3 B4];
   index = find(counts < 1);
   if(isempty(index) ~=1)
      counts(index) = 0.0000000001;
      X = counts(1);
      A0 = counts(2);
      A1 = counts(3);
      A2 = counts(4);
      A3 = counts(5);
      A4 = counts(6);
      B0 = counts(7);
      B1 = counts(8);
      B2 = counts(9);
      B3 = counts(10);
      B4 = counts(11);
   end;

   % The likelihood of any reaction is K*Q, where Q is written as reactants/products
   % A4 +X^2 --> X^3
   f1 = K1*(A4)/X;
   % X --> B
   r2 = (1/K2)*X/B4;

   %A0 --> A1
   f3 = A0/A1;
   f4 = (A1)/A2;
   f5 = (A2)/(A3);
   f6 = (A3)/(A4);

   % B1 --> B0
   fb3 = B1/(B0);
   fb4 = (B2)/(B1);
   fb5 = (B3)/(B2);
   fb6 = (B4)/(B3);
return;
