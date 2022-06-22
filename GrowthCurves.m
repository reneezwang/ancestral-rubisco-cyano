%% Growth Curves for Wang et al (2022) PNAS
%%%%% This file:
%%%%% 1) Plots growth curves
%%%%% 2) Calculates best-fit exponential growth constant


%% Get workspace ready
%%%%% Clear workspace if necessary
clear all
clf
close all
%%%%% CHANGE ME to your current working directory
cd '/Volumes/ReneeWang/PHD_ALL/Rubisco_ALL/ANCRubisco_forGithub'
%%%%% Import growth curve data
growthcurves = importdata('allCurves_smoothed.xlsx');


%% Set plotting niceties 
%%%%% Set colors
% Set blues
blues = ["#4D80E6" "#99B3FF" "#004DE6" "#0073BD"];
% Set greens
greens = ["#00B366" "#00804D" "#00994D" "#00B300"];
% Set grays
grays = ["#666666" "#000000"];
% Marker size
sz = 10;
% Line width
w = 1;


%% Plot all growth curves (Fig. S15 in Supplemental)
figure
sgtitle('Figure S15: Supplemental Growth Curve Plots');

subplot(2,2,1)
hold on
%%%%% WT Reference Condition in blue
for i = 1:4
    plot(growthcurves.data(:,1),growthcurves.data(:,(i+1)),...
        'Color',blues(i),...
        'LineStyle','-',...
        'LineWidth',2,...
        'Marker','none');
end
%%%%% WT High CO2 Condition in green
for i = 1:4
    plot(growthcurves.data(:,1),growthcurves.data(:,(i+8)),...
        'Color',greens(i),...
        'LineStyle','-',...
        'LineWidth',2,...
        'Marker','none');
end
%%%%% WT High Light Condition in gray
hold on
for i = 1:2
    plot(growthcurves.data(:,1),growthcurves.data(:,(i+16)),...
        'Color',grays(i),...
        'LineStyle','-',...
        'LineWidth',2,...
        'Marker','none');
end
title('WT')
ylabel('Abs at 720 nm');
xlabel('Time (hrs)');
ylim([0 1]);
box on
hold off

subplot(2,2,2)
hold on
%%%%% ANC Reference Condition in blue
for i = 1:3
    plot(growthcurves.data(:,1),growthcurves.data(:,(i+5)),...
        'Color',blues(i),...
        'LineStyle','-',...
        'LineWidth',2,...
        'Marker','none');
end
%%%%% ANC High CO2 Condition in green
for i = 1:4
    plot(growthcurves.data(:,1),growthcurves.data(:,(i+12)),...
        'Color',greens(i),...
        'LineStyle','-',...
        'LineWidth',2,...
        'Marker','none');
end
%%%%% ANC High Light Condition in gray
for i = 1:2
    plot(growthcurves.data(:,1),growthcurves.data(:,(i+18)),...
        'Color',grays(i),...
        'LineStyle','-',...
        'LineWidth',2,...
        'Marker','none');
end
title('ANC')
xlabel('Time (hrs)');
ylabel('Abs at 720 nm')
ylim([0 1]);
box on
hold off


%% Get best fit for exponential phase of growth cuve
%%%%% This uses the GetExpFit_MCMC.m function

%%%%% Instrumental error of spectrometer
s = 0.00002*(ones(length(growthcurves.data),1));

%%%%% Initial parameters (cannot be 0)
params0.a = 0.05; % pre-exponential factor 'a'
params0.loss_frequency = 0.02; % 1/seconds; also known as 'k'
params0.left_bound = 0.01; % left bound
params0.right_bound = 50; % right bound
params0.offset = 0.01; % offset from 0

step_size = 1e-2; % step size

params0 = assembleParamVector(params0); % convert the structure into a vector for function 

step = step_size * params0; % define step size vector to be step_size percent of initial guess

%%%%%% Bacterial Growth Curve of interest
B = growthcurves.data(:,10);
B_name = growthcurves.textdata(1,10);

n_steps = 1000000; % number of steps; should run 100,000 to 1,000,000 steps

num = length(B(~isnan(B))); % exclude NaN, length of y-data
x = growthcurves.data(1:num,1); % define x-values dependent on length of y-data
y = B(1:num); % define y-values 

%%%%% Run MCMC
results = GetExpFit_MCMC(params0,step,x,y,s,n_steps);

params = vec2ParamsStruct(results.params); % convert to struct


%% Plot Parameter Results
figure
% subplot(6,1,1)
subplot(3,2,1)
hold on
plot(results.chi2);
xlabel('Iteration number');
ylabel('Calculated Chi2');
title('Chi2 values over iterations','FontSize',14)
hold off

% subplot(6,1,2)
subplot(3,2,2)
hold on
histfit(results.loc(:,1));
PDFa = fitdist(results.loc(:,1),'Normal');
title('Pre-exponential value (a)','FontSize',14);
xlabel('Frequency');
ylabel('Value');
hold off

% subplot(6,1,3)
subplot(3,2,3)
hold on
histfit(results.loc(:,2));
PDFloss_frequency = fitdist(results.loc(:,2),'Normal');
title('Exponential term (k)','FontSize',14);
xlabel('Frequency');
ylabel('Value');
hold off

% subplot(6,1,4)
subplot(3,2,4)
hold on
histfit(results.loc(:,3));
PDFleft_bound = fitdist(results.loc(:,3),'Normal');
title('Left-bound','FontSize',14);
xlabel('Frequency');
ylabel('Value');
hold off

% subplot(6,1,5)
subplot(3,2,5)
hold on
histfit(results.loc(:,4));
PDFright_bound = fitdist(results.loc(:,4),'Normal');
title('Right-bound','FontSize',14);
xlabel('Frequency');
ylabel('Value');
hold off

% subplot(6,1,6)
subplot(3,2,6)
hold on
histfit(results.loc(:,5));
PDFoffset = fitdist(results.loc(:,5),'Normal');
title('Offset (b)','FontSize',14);
xlabel('Frequency');
ylabel('Value');
hold off


%% Plot best-fit line
% Equation of best fit line
yBest = @(x) (PDFa.mu*exp(PDFloss_frequency.mu*x)+PDFoffset.mu);
yMax = @(x) ((PDFa.mu+PDFa.sigma)*exp((PDFloss_frequency.mu+PDFloss_frequency.sigma)*x)+(PDFoffset.mu+PDFoffset.sigma));
yMin = @(x) ((PDFa.mu-PDFa.sigma)*exp((PDFloss_frequency.mu-PDFloss_frequency.sigma)*x)+(PDFoffset.mu-PDFoffset.sigma));

figure
hold on

ax = gca;
ax.FontSize = 16;
title(B_name,'with best-fit model','FontSize',18);
xlabel('Time (hours)','FontSize',16);
ylabel('Absorbance at 750 nm','FontSize',16);

% Plot original growth curve
scatter(x,y,sz,'MarkerEdgeColor',greens(1));

% Plot best-fit line
fplot(yBest,'k');
fplot(yMax,':k');
fplot(yMin,':k');

% Plot left-bounds
xline(PDFleft_bound.mu,'b');
xline(PDFleft_bound.mu+PDFleft_bound.sigma,':b');
xline(PDFleft_bound.mu-PDFleft_bound.sigma,':b');

% Plot right-bounds
xline(PDFright_bound.mu,'r');
xline(PDFright_bound.mu+PDFright_bound.sigma,':r');
xline(PDFright_bound.mu-PDFright_bound.sigma,':r');
hold off


%% Define helper functions
function params_vector = assembleParamVector(params_struct)
%%% Converts the params struct to a vector
% input
% 1. params_struct: a struct containing all the fields to be fitted in the optimization code
%
% output:
% a vector to be fitted

params_vector = [params_struct.a, params_struct.loss_frequency, params_struct.left_bound, params_struct.right_bound, params_struct.offset];
end

function params_struct = vec2ParamsStruct(params_vector)
params_struct.a = params_vector(1);
params_struct.loss_frequency = params_vector(2);
params_struct.left_bound = params_vector(3);
params_struct.right_bound = params_vector(4);
params_struct.offset = params_vector(5);
end

