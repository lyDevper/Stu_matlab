%{  
This script for prepare data and parameters for parameter estimator.
1. Load your collected data to MATLAB workspace.
2. Run this script.
3. Follow parameter estimator instruction.
%}

% R and L from experiment
motor_R = 1.12727; %1.58679;
motor_L = 0.0291360334; %0.000538056;
% Optimization's parameters
motor_Eff = 0.5;
motor_Ke = 0.05;
motor_J = 0.001;
motor_B = 0.0001;

% Extract collected dataclear

%pulse2s
a = 1;
Input1 = data{a}{1}.Values.Data;
Time1 = data{a}{1}.Values.Time;
Velo1 = double(data{a}{3}.Values.Data);

%sin1
a = 5;
Input2 = data{a}{1}.Values.Data;
Time2 = data{a}{1}.Values.Time;
Velo2 = double(data{a}{3}.Values.Data);

%stair1
a = 8;
Input3 = data{a}{1}.Values.Data;
Time3 = data{a}{1}.Values.Time;
Velo3 = double(data{a}{3}.Values.Data);

% Plot 
figure(Name='Motor velocity response')
plot(Time3,Velo3,Time3,Input3)