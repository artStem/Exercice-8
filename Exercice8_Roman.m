%% Générale : executer une fois avant les simulations ====================%
clear all
close all
clc

repertoire = ''; 
executable = 'Exercice8';

% Plot parametters ==================================
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 15);
set(groot, 'defaultLineLinewidth', 1.5);

% Constantes ========================================

pi = 3.14159265358979323846264338327950288419716939937510582097494459230;
% g = 9.80665;
% exp = 2.718281828459045;
% e =  1.602176634e-19;
% c = 299792458;
% epsilon0 = 8.85418782e-12;
% mu0 = 4*pi*1e-7;
% G = 6.674e-11;
%% Etude de convergence ================================================= %
input = 'configuration2i.in';
doSimulation = 0;           % 1 = TRUE    0 = FALSE
loadData = 0;               % Pour grosse simulation déjà chargée dans MATLAB
%=========================================================================%
Nsimul = 15;
dt = logspace(-3,-1,Nsimul);
%=========================================================================%
nomFigure{1} = '8_2i_convergence_x_moy';
nomFigure{2} = '8_2i_convergence_p_moy';
nomFigure{3} = '8_2i_convergence_incert_x';
nomFigure{4} = '8_2i_convergence_incert_y';

%=========================================================================%
%=========================================================================%


if loadData == 1
    for i = size(dt')
        nomfile{i,1} = ['dt=', num2str(dt(i)), '_pot'];
        nomfile{i,2} = ['dt=', num2str(dt(i)), '_sqwave'];
        nomfile{i,3} = ['dt=', num2str(dt(i)), '_obs'];
        if doSimulation == 1
            % Simulation
            cmd = sprintf('%s%s %s output_potential=%s output_squared_wave=%s output_observables=%s dt=%.15g', ...
                repertoire, executable, input, nomfile{i,1}, nomfile{i,2}, nomfile{i,3}, dt(i));     % Pour rajouter parametre numérique : %s=%.15g
            disp(cmd)
            system(cmd);
        end
    end
end

% p = polyfit(log10(CFL(1:end-1)),log10(erreur(1:end-1)),1);
% xfit = linspace(1e-2, 1e0);
% yfit = (10^(p(2)))*xfit.^(p(1));
% legende = ['Fit : $\log (\epsilon)$ = ', num2str(p(1),'%.1f'), '$\log(\beta_{CFL})$ + ', num2str(p(2),'%.1f')];

% fig{1} = figure('Name',nomFigure{5},'NumberTitle','off');
% loglog(CFL(1:end-1), erreur(1:end-1), 'linestyle', 'none', 'marker', '+', 'Displayname', 'Numerical values');
% hold on
% loglog(xfit,yfit, 'k--', 'Displayname', legende);
% xlabel('$\beta_{CFL}$');
% ylabel('$\epsilon$ [m$^2$]');
% legend('location', 'best');
% grid on

%% Save figures ==========================================================%
saveFigure = 0;         %1 = TRUE, 0 = FALSE
sauver = 0;%[1 2 3 4 5];

if saveFigure == true
    for i=sauver
        exportgraphics(fig{i},[nomFigure{i}, '.pdf'],'ContentType','vector');
    end
end


%% Fonctions utilisées


