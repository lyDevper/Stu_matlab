%{  
This script for prepare data and parameters for parameter estimator.
1. Load your collected data to MATLAB workspace.
2. Run this script.
3. Follow parameter estimator instruction.
%}

% R and L from experiment
motor_R = 0.92384; %1.58679;
motor_L = 0.0877797954; %0.0291360334; %0.000538056;
% Optimization's parameters
%motor_Km = 0.0527;
motor_Efficiency = 0.8;
motor_Ke = 0.05;
motor_J = 0.001;
motor_B = 0.0001;

% Extract collected data
a = 5;
Input = data{a}{1}.Values.Data;
Time = data{a}{1}.Values.Time;
Velo = double(data{a}{3}.Values.Data);

% down size
Input = Input(1:10:end);
Time = Time(1:10:end);
Velo = Velo(1:10:end);

% Plot 
figure(Name='Motor velocity response')
plot(Time,Velo,Time,Input)