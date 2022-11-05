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
c_rel = c0./sqrt(mu_rel.*e_rel); % relative velocity c in each layer
c_abs= c0./sqrt(mu_abs.*e_abs);  % absolute velocity c in each layer

% Plotting parameters
scalefactor = 0.25;  % color bar scales to 0.25 of the maximum and minimum amplitude
fontsize = 16;       % fontsize used for titles    
fontstyle ='Arial';  % fontstyle used for titles

%% Task 2: Discretizing Domain, computing maximum velocity and maximum step size
tstep= (1/2*max(fmax));         % time step is max when t=(1/2fmax)
stepsize= max(c_rel)*tstep;     % maximum step size is when time step is max
max_stepsize =1/3*min(d)*0.99;  % we multiply by 1/3 to make sure that each layer is sampled for atleast 3 points


% It is then checked whether the z_step is small enough and if not the
    % z_step and t_step are adjusted
        if stepsize < max_stepsize
            disp('The step size of z is small enough so that also the smallest layer is sampled with three or more points')
        else
            stepsize = max_stepsize;
            disp('The step size of z was too large and was therefore decreased to the maximum possible value')  
            tstep = stepsize/max(c_rel); % when the z_step is decreased we also need to decrease the t_step, because otherwise the condition for z_step>= c_max*t_step is not fulfilled anymore.
        end

% Make a discrete model and define values of epsilon, sigma and mu for all depth points:
d_sum = cumsum(d)-d(1); % depth of all interfaces
t = linspace(0,tw,(tw/tstep)); % creating a time vector
z = linspace(0,d_sum(end),(d_sum(end)/stepsize)); % creating a depth vector
n_t = numel(t); % calculate the number of time points in the model
n_z = numel(z); % calculate the number of depth points in the model

% create empty vectors for e, s and mu:
e_z = zeros(n_z,1); 
s_z = zeros(n_z,1);
mu_z = zeros(n_z,1);

for i = 2:length(d_sum) % loop over all layers of the domain (corresponds to the elements 2:end of the depth intervals)
        for j = 1:n_z
            if z(j) <= d_sum(i) && z(j) >= d_sum(i-1) 
                e_z(j) = e_abs(i);
                mu_z(j) = mu_abs(i);
                s_z(j) = s(i);
            end
        end
end
%% Task 3: Precompute the incident electric field for all times

% Definition of useful constants/variables:

t_arg = t - d(1)/c_abs(1); % This is the argument of the wavelet
tau = sqrt(2)./fc;

derv_Elect = zeros(1,numel(t), numel(fc));

for i = 1:numel(fc)
    derv_Elect(1,:,i) = z0*(2*pi^2).*(1./tau(i))*((t_arg./tau(i))-1).*exp(-2*pi^2*(t_arg./tau(i)-1).^2).*(3-4*pi^2*(t_arg./tau(i)-1).^2);
end

%% Task 4: 
% Write the finite-difference code to calculate the electric and magnetic
% fields at all grid points for each time step
 
%Intialise the electric and magnetic field
E = zeros(numel(z),numel(t)+2,numel(fc)); %Also electric field at time is 0
H = zeros(numel(z),numel(t)+2,numel(fc));

%Creating useful constants needed for equation 37: 
denom = (mu_abs(1)/mu_abs(2)) + (stepsize/(c_abs(1)*tstep));
const_1 = (mu_abs(1)/mu_abs(2))/3;
const_2 = (stepsize/(c_abs(1)*tstep))/3;
const_3 = (2/3)*stepsize/c_abs(1);
const_4 = tstep/stepsize;


for f = 1:numel(fc) % all three center frequencies

    for time = 3:numel(t) % for all times but equation 37 requires a value for time-2 so we start at 3 for numel to work
        % Electric field for 1st depth, defined by boundary condition 37:

        E(1, time, f) = (const_1*(4*E(1+1,time,f)- E(1+2,time,f)) + const_2*(4*E(1,time-1,f) - E(1,time-2,f))+ const_3*derv_Elect(1,time,f))/denom;
        for depth = 1:numel(z)-1 %for all depth points
            %Calculate Magnetic field H based on equation 27, and using E(1) defined above for initial value of H(1):
            H(depth, time, f) = H(depth, time-1, f) + (const_4/mu_z(depth))*(E(depth+1,time,f) - E(depth,time,f));

            %Calculating Electric field for all depths based on equation 26

            %Constants needed as a function of depth
            const_5 = e_z(depth) - s_z(depth)*tstep/2;
            const_6 = e_z(depth) + s_z(depth)*tstep/2;

            %Equation 26
             if depth >= 2
                %Use Equation 26 to calculate all electric fields inside
                %the boundaries
                E(depth,time+1,f) = (const_5/const_6)*E(depth,time,f) + (const_4/const_6)*(H(depth,time,f) - H(depth-1,time,f));
             end
             
        end
        

        %The electric field on the boundary (last depth step, numel(z)) is defined with equation 34. So we re-define it here: 
        const_7 = (mu_abs(6)/mu_abs(5));
        const_8 = (1/3)*const_7;
        const_9 = (stepsize/(c_abs(6)*tstep));
        const_10 = (1/3)*const_9;
        const_11 = s(6)*sqrt(mu_abs(6)/e_abs(6))*stepsize/3;
   
        E(numel(z), time+1, f) = (const_8*(4*E(numel(z)-1,time+1,f)-E(numel(z)-2,time+1,f)) + const_9*(4*E(numel(z),time,f) -E(numel(z),time-1,f)))/(const_7 + const_9 + const_11);
        
    end
    
end




%% Task 5: No explicit loops over depth (ie. use array operations)
%%We repeat Task 4 (iterations of center of frequencies and time) until the
%%definition of the 1st electric field boundary condition.
%%We won't re-define the constant used outside of the loops - see Task 4 for explanations.

%Intialise the electric and magnetic field using different variable name to
%compare with task 4.
E_T5 = zeros(numel(z),numel(t)+2,numel(fc)); %Also electric field at time is 0 (causal)
H_T5 = zeros(numel(z),numel(t)+2,numel(fc));

for f = 1:numel(fc) 
    for time = 3:numel(t) 
        E_T5(1, time, f) = (const_1*(4*E_T5(1+1,time,f)- E_T5(1+2,time,f)) + const_2*(4*E(1,time-1,f) - E_T5(1,time-2,f))+ const_3*derv_Elect(1,time,f))/denom; %from equation 37

        %H at depth 1 has for equation: 
        H_T5(1, time, f) = H_T5(1, time-1, f) + (const_4/mu_z(1))*(E_T5(1+1,time,f) - E_T5(1,time,f)); %from equation 27 
        
        %H for further depths (less than boundary d):
        idx_depth = (2:z);
        H_T5(idx_depth, time, f) = H_T5(idx_depth, time-1, f) + (const_4/mu_z(idx_depth))*(E(idx_depth+1,time,f) - E(idx_depth,time,f));
        %For depth not boundaries
        E(idx_depth,time+1,f) = (const_5/const_6)*E(idx_depth,time,f) + (const_4/const_6)*(H(idx_depth,time,f) - H(idx_depth-1,time,f));

    end
    %for electric field on final boundary:
    %The electric field on the boundary (last depth step, numel(z)) is defined with equation 34. So we re-define it here: 
    const_7 = (mu_abs(6)/mu_abs(5)); %index 5 corresponds to parameters of the scattering domain & index 6 corresponds to parameters of the Lower half space
    const_8 = (1/3)*const_7;
    const_9 = (stepsize/(c_abs(6)*tstep));
    const_10 = (1/3)*const_9;
    const_11 = s(6)*sqrt(mu_abs(6)/e_abs(6))*stepsize/3; 
   
    E_T5(numel(z), time+1, f) = (const_8*(4*E_T5(numel(z)-1,time+1,f)-E_T5(numel(z)-2,time+1,f)) + const_9*(4*E_T5(numel(z),time,f) -E_T5(numel(z),time-1,f)))/(const_7 + const_9 + const_11);
end



