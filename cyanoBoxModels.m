%% Traditional and Proposed Box Model for Wang et al (2022) PNAS
%%%%% For full discussion of model solution, see Supplemental

%% Experimental Results and Model Parameters
%%%%% 1. Experimental results of epsilon P across 3 conditions for WT and ANC
%%%%% Condition 1: Reference Condition
%%%%% Condition 2: High CO2
%%%%% Condition 3: High Light
%%%%% For clarity, errors are not shown because they are smaller than the
%%%%% marker size. 

epsP_WTCondition1 = 7.453;
epsP_WTCondition2 = 18.304;
epsP_WTCondition3 = 7.812;

epsP_ANCCondition1 = 18.247;
epsP_ANCCondition2 = 19.467;
epsP_ANCCondition3 = 24.299;

%%%%% 2. Model parameters
%%%%% a) Fractionations used for both strains
epsIn = 1; % Diffusion
epsPCA = 30; % Powered Carbonic Anhydrase
epsL1 = 0; % Loss 1; Assume no fractionation during loss process
epsL2 = 0; % Loss 2; Assume no fractionation during loss process

%%%%% Epsilon Rubisco specific to each strain
epsRub_WT = 25.2; %WT = 25.2 +/- -.3 std err
epsRub_ANC = 17.2; %ANC = 17.2 +/- 0.6 std err

%% Traditional box model

%%%%% Define model for each strain (They have different Rubisco fractionations)
epsP_WT_traditional = @(f) (1-f)*epsIn + f*epsRub_WT;
epsP_ANC_traditional = @(f) (1-f)*epsIn + f*epsRub_ANC;

%%%%% Solve for f given the epsilon P values measured experimentally
f_WTCondition1 = (epsP_WTCondition1-epsIn)/(epsRub_WT-epsIn);
f_WTCondition2 = (epsP_WTCondition2-epsIn)/(epsRub_WT-epsIn);
f_WTCondition3 = (epsP_WTCondition3-epsIn)/(epsRub_WT-epsIn);
f_ANCCondition1 = (epsP_ANCCondition1-epsIn)/(epsRub_ANC-epsIn);
f_ANCCondition2 = (epsP_ANCCondition2-epsIn)/(epsRub_ANC-epsIn);
f_ANCCondition3 = (epsP_ANCCondition3-epsIn)/(epsRub_ANC-epsIn);

%% Proposed box model

%%%%% Create values of f to input
f1_ANC = linspace(0,1,100);
f2_ANC = linspace(0,1,100);

f1_WT = linspace(0,1,100);
f2_WT = linspace(0,1,100);

%%%%% Solve for possible values of Epsilon P for ANC
j1 = 1;
while j1 <= length(f1_ANC)
    for k = 1:length(f2_ANC)
        epsP_ANC(j1,k) = epsL2 - epsL1 - epsIn + f1_ANC(j1)*(epsPCA-epsL1) + f2_ANC(k)*(epsRub_ANC-epsL2);
    end
    j1 = j1+1;
end

%%%%% Solve for possible values of Epsilon P for WT
j2 = 1;
while j2 <= length(f1_WT)
    for k = 1:length(f2_WT)
        epsP_WT(j2,k) = epsL2 - epsL1 - epsIn + f1_WT(j2)*(epsPCA-epsL1) + f2_WT(k)*(epsRub_WT-epsL2);
    end
    j2 = j2+1;
end

%% Visualize Traditional and Proposed Box Model
%%%%% This is the figure seen in the main text

figure
sgtitle('Main Text: Box Model Outputs');

%%%%% Traditional Box Model Outputs
subplot(1,3,1)
hold on

%%%%% Plot possible model outputs
fplot(epsP_WT_traditional,[0 1],'k');
fplot(epsP_ANC_traditional,'--k');

%%%%% Plot experimental results
sz = 70; % Size of markers
%%%%% WT results
scatter(f_WTCondition1,epsP_WTCondition1,sz,...
    'o',...
    'MarkerFaceColor','blue',...
    'MarkerEdgeColor','blue',...
    'MarkerFaceAlpha',0.5);
scatter(f_WTCondition2,epsP_WTCondition2,sz,...
    'o',...
    'MarkerFaceColor',[0.3,0.5,0.11],...
    'MarkerEdgeColor',[0.3,0.5,0.11],...
    'MarkerFaceAlpha',0.5);
scatter(f_WTCondition3,epsP_WTCondition3,sz,...
    'o',...
    'MarkerFaceColor','black',...
    'MarkerEdgeColor','black',...
    'MarkerFaceAlpha',0.5);
%%%%% ANC results
scatter(f_ANCCondition1,epsP_ANCCondition1,sz,...
    '+',...
    'MarkerEdgeColor','blue',...
    'LineWidth',1.5);
scatter(f_ANCCondition2,epsP_ANCCondition2,sz,...
    '+',...
    'MarkerEdgeColor',[0.3,0.5,0.11],...
    'LineWidth',1.5);
scatter(f_ANCCondition3,epsP_ANCCondition3,sz,...
    '+',...
    'MarkerEdgeColor','black',...
    'LineWidth',1.5);

%%%%% Figure formatting
title('Traditional box model','FontSize',11);
xlabel('f = \phi_L_o_s_s / \phi_I_n','FontSize',11);
ylabel('\epsilon_P','FontSize',11);
xlim([0 1.5]);
ylim([0 30]);
xline(1);
axis square
box on
hold off



%%%%% Place-holder for cartoon of proposed box model
subplot(1,3,2)
axis square
box on



%%%%% Proposed Box Model
%%%%% Assume f = 0.1
ANC_f1_10pct = epsP_ANC(:,10); %10th column is 0.1
WT_f1_10pct = epsP_WT(:,10); %10th column is 0.1
x_newPlot = linspace(0,1,100);

%%%%% Calculate f2 values for epsP results
%%%%% polyfit gives m & b of line
fits_WT_f1_10pct = polyfit(x_newPlot,WT_f1_10pct,1);
fits_ANC_f1_10pct = polyfit(x_newPlot,ANC_f1_10pct,1);
%%%%% Solve for f2 (x=(y-b)/m))
WT_newX1 = (epsP_WTCondition1-fits_WT_f1_10pct(2))/fits_WT_f1_10pct(1);
WT_newX2 = (epsP_WTCondition2-fits_WT_f1_10pct(2))/fits_WT_f1_10pct(1);
WT_newX3 = (epsP_WTCondition3-fits_WT_f1_10pct(2))/fits_WT_f1_10pct(1);
ANC_newX1 = (epsP_ANCCondition1-fits_ANC_f1_10pct(2))/fits_ANC_f1_10pct(1);
ANC_newX2 = (epsP_ANCCondition2-fits_ANC_f1_10pct(2))/fits_ANC_f1_10pct(1);
ANC_newX3 = (epsP_ANCCondition3-fits_ANC_f1_10pct(2))/fits_ANC_f1_10pct(1);

subplot(1,3,3)
hold on
% Plot model results
plot(x_newPlot,WT_f1_10pct,'k');
plot(x_newPlot,ANC_f1_10pct,'--k');
% Plot data on top
scatter(WT_newX1,epsP_WTCondition1,sz,...
    'o',...
    'MarkerFaceColor','blue',...
    'MarkerEdgeColor','blue',...
    'MarkerFaceAlpha',0.5);
scatter(WT_newX2,epsP_WTCondition2,sz,...
    'o',...
    'MarkerFaceColor',[0.3,0.5,0.11],...
    'MarkerEdgeColor',[0.3,0.5,0.11],...
    'MarkerFaceAlpha',0.5);
scatter(WT_newX3,epsP_WTCondition3,sz,...
    'o',...
    'MarkerFaceColor','black',...
    'MarkerEdgeColor','black',...
    'MarkerFaceAlpha',0.5);

scatter(ANC_newX1,epsP_ANCCondition1,sz,...
    '+',...
    'MarkerEdgeColor','blue',...
    'LineWidth',1.5);
scatter(ANC_newX2,epsP_ANCCondition2,sz,...
    '+',...
    'MarkerEdgeColor',[0.3,0.5,0.11],...
    'LineWidth',1.5);
scatter(ANC_newX3,epsP_ANCCondition3,sz,...
    '+',...
    'MarkerEdgeColor','black',...
    'LineWidth',1.5);

title({'Proposed box model','(f_1=0.1)'},'FontSize',11);
xlabel('f_2 = \phi_L_o_s_s_2 / \phi_N_D_H','FontSize',11);
ylabel('\epsilon_P','FontSize',11);
xlim([0 1]);
axis square
box on
hold off


%% Visualize Full Model Outputs for WT

figure
hold on

%%%%% Plot WT model outputs
contourf(f1_WT,f2_WT,epsP_WT,'LineColor','none');

%%%%% Plot WT experimental results onto model outputs
contour(f1_WT,f2_WT,epsP_WT,[epsP_WTCondition1 epsP_WTCondition1],'--w','LineWidth',1.5);
contour(f1_WT,f2_WT,epsP_WT,[epsP_WTCondition2 epsP_WTCondition2],':w','LineWidth',1.5);
contour(f1_WT,f2_WT,epsP_WT,[epsP_WTCondition3 epsP_WTCondition3],'-w','LineWidth',1.5);

%%%%% Show f1 = 0.1
xl = xline(0.1,'-y','f_1 = 0.1');
xl.LineWidth = 2;
xl.FontSize = 14;

%%%%% Figure formatting
title('WT Proposed Box Model: Full Outputs','FontSize',14);
xlabel('f_1 = \phi_l_o_s_s_1/\phi_i_n','FontSize',14);
ylabel('f_2 = \phi_l_o_s_s_2/\phi_N_D_H','FontSize',14);
c = colorbar;
c.Label.String = '\epsilon_P';
c.Label.FontSize = 14;
colormap(winter)
box on
ax = gca;
ax.LineWidth = 2;
axis square
hold off


%% Visualize Full Model Outputs for ANC

figure
hold on

%%%%% Plot ANC model outputs
contourf(f1_ANC,f2_ANC,epsP_ANC,'LineColor','none');

%%%%% Plot ANC experimental results onto model outputs
contour(f1_ANC,f2_ANC,epsP_ANC,[epsP_ANCCondition1 epsP_ANCCondition1],'--k','LineWidth',1.5);
contour(f1_ANC,f2_ANC,epsP_ANC,[epsP_ANCCondition2 epsP_ANCCondition2],':k','LineWidth',1.5);
contour(f1_ANC,f2_ANC,epsP_ANC,[epsP_ANCCondition3 epsP_ANCCondition3],'-k','LineWidth',1.5);

%%%%% Show f1 = 0.1
xl = xline(0.1,'-r','f_1 = 0.1');
xl.LineWidth = 2;
xl.FontSize = 14;

% Figure formatting
title('ANC Proposed Box Model: Full Outputs','FontSize',14);
xlabel('f_1 = \phi_l_o_s_s_1/\phi_i_n','FontSize',14);
ylabel('f_2 = \phi_l_o_s_s_2/\phi_N_D_H','FontSize',14);
c = colorbar;
c.Label.String = '\epsilon_P';
c.Label.FontSize = 14;
colormap(summer)
box on
ax = gca;
ax.LineWidth = 2;
axis square
hold off



