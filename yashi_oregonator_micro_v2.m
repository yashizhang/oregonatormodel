clear all
close all 

% Author: Ya Shi (Andrew) Zhang - yashi.zhang@nyu.edu
% Simulate the Belousov-Zhabotinsky Reaction through
% the Oregonator scheme for chemical reactions
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

refV = 1e17; % Reference Volume
tmax = 50; % Maximum Time
f    = 1; % Stoichiometric Factor
NA   = 6.02214076e23; % Avogadro's Constant
K    = [1.34 1.6e9 8e3 4e7 1*(NA/refV)] / (NA/refV); % Vector of Macro Rate Constants [K1 K2 K3 K4 K5] (1/(#Molecules*Sec) and 1/sec for K5)
V    = 4^3; % Volume of System  
N    = zeros(1,5); % State of system
% Initial Number of Chemicals (Converted from dimensionless variables in
% [1]). [alpha eta row] = [488.68 0.99796 488.68] is unstable steady state.
alpha = 500;
eta   = 1;
row   = 500;
N(1) =         0.06      * (NA/refV) * V; N(1) = round(N(1)); % #Molecules of A 
N(2) =         0.06      * (NA/refV) * V; N(2) = round(N(2)); % #Molecules of B
N(3) = alpha * 5.025e-11 * (NA/refV) * V; N(3) = round(N(3)); % #Molecules of X
N(4) = eta   * 3e-7      * (NA/refV) * V; N(4) = round(N(4)); % #Molecules of Y
N(5) = row   * 2.412e-8  * (NA/refV) * V; N(5) = round(N(5)); % #Molecules of Z

% Change in N when reaction j happens (from I to V)
dN      = zeros(5,5); 
dN(:,1) = [-1 0 1 -1 0];
dN(:,2) = [0 0 -1 -1 0];
dN(:,3) = [0 -1 1 0 1];
dN(:,4) = [0 0 -2 0 0];
dN(:,5) = [0 0 0 f -1];

timestate   = zeros(1e8,1);
systemstate = zeros(1e8,5);

% Main Loop
t = 0;
eventcounter = 1;
while t < tmax
    rate = rates(N,K,V);
    T = (-log(rand(1,5)))./rates(N,K,V);
    [Tmin, Kmin] = min(T);
    N = N + dN(:,Kmin)';
    systemstate(eventcounter,:) = N;
    t = t + Tmin;
    timestate(eventcounter) = t;
    eventcounter = eventcounter + 1;
end

labels     = ["logNX", "logNY", "logNZ"];
stackedplot(timestate(1:eventcounter-1),...
    log(systemstate(1:eventcounter-1,3:5)),...
    "Title", "Microscopic: log(#Molecules) vs. Time",...
    "DisplayLabels", labels, "xlabel", "Time")
saveas(gcf,'micro.png');

%figure
%plot(timestate(1:eventcounter-1), systemstate(1:eventcounter-1,3), ...
%    'DisplayName', 'X', 'Color', 'r');
%title('NX vs. Time');
%saveas(gcf,'microX.png');
%figure
%plot(timestate(1:eventcounter-1), systemstate(1:eventcounter-1,4), ...
%    'DisplayName', 'Y', 'Color', 'g');
%title('NY vs. Time');
%saveas(gcf,'microY.png');
%figure
%plot(timestate(1:eventcounter-1), systemstate(1:eventcounter-1,5), ...
%    'DisplayName', 'Z', 'Color', 'b');
%title('NZ vs. Time');
%saveas(gcf,'microZ.png');

function rates = rates(N,K,V)
% Takes in the current state of the system (N)
% as well as the vector of reaction coefficients (K)
% to output the probability per unit time that each 
% reaction will occur.
% This specific rate function is for the Oregonator 
% class of oscillating chemical reactions 

rate1 = abs(K(1) * N(1) *  N(4)   ) / V; % (K1/V)*NA*NY
rate2 = abs(K(2) * N(3) *  N(4)   ) / V; % (K2/V)*NX*NY
rate3 = abs(K(3) * N(2) *  N(3)   ) / V; % (K3/V)*NB*NX
rate4 = abs(K(4) * N(3) * (N(3)-1)) / V; % (K4/V)*(NX)*(NX-1)
rate5 = abs(K(5) * N(5));                % K5*NZ

rates = [rate1 rate2 rate3 rate4 rate5];
end

% Sources: Refer to paper