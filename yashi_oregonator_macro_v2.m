clear all
close all 

% Author: Ya Shi (Andrew) Zhang - yashi.zhang@nyu.edu
% Simulate the Belousov-Zhabotinsky Reactions through
% the Oregonator scheme
% i.e. reactions of the form:
%   (I) A + Y -> X       | Rate K1
%  (II) X + Y -> P       | Rate K2
% (III) B + X -> 2X + Z  | Rate K3
%  (IV) 2X -> Q          | Rate K4
%   (V) Z -> fY          | Rate K5
% Where P, Q are products of reaction, f is the stoichiometric
% factor (defaulted to 1 in literature), X is HBrO2, Y is Br(-),
% Z is Ce(4+), and both A and B are BrO3(-). P and Q do not need 
% to be modeled as they are products of the reaction.

% Reference Volume
refV = 1e16;
% Avogadro's Constant
NA    = 6.02214076e23;
% Time 
tspan = [0 500];
% Initial Concentration of [X Y Z] measured in #Molecules/mm^3
alpha = 500;
eta   = 1;
row   = 500;
y0    = [5.025e-11*alpha 3e-7*eta 2.412e-8*row] * (NA/refV);
%[alpha eta row] = [488.68 0.99796 488.68] is an unstable steady state

[t,y] = ode15s(@Oregonator,tspan,y0); % Dynamic Eq. is stiff

% Normal log concentration plot
labels = ["log(X)", "log(Y)", "log(Z)"];
stackedplot(t,log(y), "Title", "Macroscopic: log(#Molecules) vs. Time",...
    "DisplayLabels", labels, "xlabel", "Time")
saveas(gcf,'macro.png');

%figure
%logy = log(y);
% 3d trajectory plot
%plot3(logy(:,1),logy(:,2),logy(:,3));


function dydt = Oregonator(t,y)
    % Reference Volume
    refV = 1e16;
    % Avogadro's Constant
    NA = 6.02214076e23;
    % Parameter values taken from Oscillations in chemical systems
    % Vector of macro rate constants (1/(#Molecules*Sec))
    k  = [1.34 1.6e9 8e3 4e7 1*(NA/refV)] / (NA/refV);
    % Concentrations of [A] and [B] in #Molecules/mm^3
    A  = 0.06 * (NA/refV);
    B  = 0.06 * (NA/refV);
    % Stoichiometric constant of Y in reaction V
    f  = 1;

    dydt    = zeros(3,1);
    dydt(1) = k(1)*A*y(2) - k(2)*y(1)*y(2) + k(3)*B*y(1) - 2*k(4)*y(1)*y(1);
    dydt(2) = -k(1)*A*y(2) - k(2)*y(1)*y(2) + f*k(5)*y(3);
    dydt(3) = k(3)*B*y(1) - k(5)*y(3);
end

% Sources: Refer to paper