%{
Title: Assignment 4- Eigenfunctions and Finite-Difference modelling
Course Code: AESM1511
Name:Sargun Kaur
Studentnumber: 5833728
%}

% Initial settings
clear variables %remove items from workspace, freeing up system memory;
close all       %close any open figures;
clc             %clear the command window

% Define constant and given parameters
c0 = 299792458;     % [m/s], velocity in free space
mu0 = 4*pi*10^(-7); % magnetic permeability in free space
e0 = 1/(mu0*c0^2);  % electric permittivity in free space
tw = 12*10^(-9);    % [s], t+--he given time window of computation
fc = [500*10^6, 1*10^9, 2*10^9];  % [Hz], given centre frequencies
max_factor = 4;     % relates the factor between centre frequency and maximum frequency fmax
fmax = max_factor*fc;  
z0 = sqrt(mu0/e0);

% Initialising layer parameters [Upper half space, Layer1 - Layer4 (Domain D2), Lower half space]
e_rel = [1, 6, 2, 16, 6, 9]; % relative electric permittivity
mu_rel = [1, 1, 1, 1, 1, 1]; % relative magnetic permeability
e_abs = e_rel*e0; % multiplication of the relative electrical permittivity and permittivity of free space gives absolute electrical permittivity
mu_abs = mu_rel*mu0;  % multiplication of the relative magnetic permeability with permittivity of free space gives absolute magnetic permeability   
s = [0, 10^(-3), 10^(-4), 10^(-2), 10^(-3), 10^(-3)]; % [S/m] electric conductivity given by sigma 
d = [0.24, 0.10, 0.003, 0.05, 0.15];  % thickness of layers (first element is the height of the source above the ground), thickness of lower halfspace is infinite
c_rel = c0./sqrt(mu_rel.*e_rel) % relative velocity c in each layer
c_abs= c0./sqrt(mu_abs.*e_abs)  % absolute velocity c in each layer

% Plotting parameters
scalefactor = 0.25;  % color bar scales to 0.25 of the maximum and minimum amplitude
fontsize = 16;       % fontsize used for titles    
fontstyle ='Arial'; % fontstyle used for titles

%% Task 2: Discretizing Domain, computing maximum velocity and maximum step size
tstep= (1/2*max(fmax))         % time step is max when t=(1/2fmax), we take max(fmax) to have the greatest time step
stepsize= max(c_abs)*tstep % maximum step size is when time step is max
max_stepsize =1/3*min(d)*0.99     % we multiply by 1/3 to make sure that each layer is sampled for atleast 3 points
if stepsize < min(d)
    msg= 'This step is smaller than the smallest layer thickness';
    error(msg);
else 
    disp("Maximum time step:")
    disp(tstep)
    disp("Maximum step size:")
    disp(max_stepsize)
end 

% Make a discrete model and define values of epsilon, sigma and mu for all depth points:
d_sum = cumsum(d)-d(1); % depth of all interfaces
t = linspace(0,tw,(tw/tstep)); % creating a time vector
z = linspace(0,d_sum(end), (d_sum(end)./stepsize)); % creating a depth vector
n_t = numel(t); % calculate the number of time points in the model
n_z = numel(z); % calculate the number of depth points in the model

%% Task 3: Precompute the incident electric field for all times

% Definition of useful constants/variables:

t_arg = t - d(1)/c_abs(1); % This is the argument of the wavelet
tau = sqrt(2)./fc;
Z_0 = sqrt(mu0/e0);

derv_Elect = zeros(1,numel(t), numel(fc));

  % We calculated by hand the (at z = 0 for all t) dW(t)/dt =
  % ((t-tau)/tau^2)exp(-2*pi^2*(t/tau-1)^2)*(-8*pi^2-4*pi^2*(1-4*pi^2)*(t/tau-1)^2)

  %By definition, the derivative of Electirc field as a function of t is 
  %derv_Elect = (-Z_0/2)*dW(t - d(1)/`c_abs(1))/dt
  % so we get for all fcs: 

for i = 1:numel(fc)
    derv_Elect(1,:,i) =(-Z_0/2).*((t-tau(i))./tau(i).^2).*exp(-2*pi^2*(t./tau(i)-1).^2).*(-8*pi^2-4*pi^2*(1-4*pi^2)*(t./tau(i)-1).^2)
end