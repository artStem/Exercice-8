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
%% Ordre des valeur dans les fichiers output
% output_potential    : x, V
% Num colonne         : 1, 2
%
% output_squared_wave : [0, x    ]
%                       [t, PSI^2]
% 
%
% output_observable   : t, P<0, P>0, P, E, <x>, <x²>, <p>, <p²>, Heisenberg, Deltax, Deltap
% Num colonne         : 1,  2 ,  3 , 4, 5,  6 ,   7 ,  8 ,   9 ,      10   ,   11  ,  12
%% Etude de convergence ================================================= %
input = 'configuration2i.in';
doSimulation = 0;           % 1 = TRUE    0 = FALSE
loadData = 0;               % Pour grosse simulation déjà chargée dans MATLAB
%=========================================================================%
tfin = 2000;
Nsimul = 30;
N = round(2*logspace(3,4,Nsimul));
%=========================================================================%
nomFigure{1} = '8_2i_convergence';

%=========================================================================%
%=========================================================================%

dt = tfin./N;
if loadData == 1
    for i = 1:size(dt')
        % Création des noms de fichier
        nomfile{i,1} = ['dt=', num2str(dt(i)), '_pot.out'];
        nomfile{i,2} = ['dt=', num2str(dt(i)), '_sqwave.out'];
        nomfile{i,3} = ['dt=', num2str(dt(i)), '_obs.out'];
        % Simulations
        if doSimulation == 1
            cmd = sprintf('%s%s %s output_potential=%s output_squared_wave=%s output_observables=%s dt=%.15g tfin=%.15g', ...
                repertoire, executable, input, nomfile{i,1}, nomfile{i,2}, nomfile{i,3}, dt(i), tfin);     % Pour rajouter parametre numérique : %s=%.15g
            disp(cmd)
            system(cmd);
        end
        % Chargement des données
        data = load(nomfile{i,3});
        % param_fin (colonnes) : <x> fin, <p> fin, Deltax fin, Deltay fin
        param_fin(i,1) = data(end,6);
        param_fin(i,2) = data(end,8);
        param_fin(i,3) = data(end,11);
        param_fin(i,4) = data(end,12);
        tf(i) = data(end,1);
    end
end

% Fit et regression linéaire sur param_fin(:,j)
% Recuperer les valeurs convergées et calculer erreur
beg = 10;
for j = 1:4
    % Convergence ordre 2
    p = polyfit(dt(beg:end).^2, param_fin(beg:end,j), 1);
    val_conv(j) = p(2);
    clear p
    % Calculer erreur
    err(:,j) = abs(1 - param_fin(:,j)./val_conv(j));
end

% Faire un fit en loglog
xfit = linspace(min(dt), max(dt));
for j = 1:4
    p = polyfit(log10(dt),log10(err(:,j)),1);
    yfit(:,j) = (10^(p(2)))*xfit.^(p(1));
    legende{j} = ['Fit slope $a$ = ', num2str(p(1),'%.2f')];
    % legende{j} = ['Fit : $\log (\epsilon)$ = ', num2str(p(1),'%.1f'), '$\log(\Delta t)$ + ', num2str(p(2),'%.1f')];
end

% Plot ===================================
nom_param = {'$\epsilon_{\langle x \rangle}$', '$\epsilon_{\langle p \rangle}$',...
    '$\epsilon_{\langle \Delta x \rangle}$', '$\epsilon_{\langle \Delta p \rangle}$'};
fig{1} = figure('Name',nomFigure{1},'NumberTitle','off');
for j=1:4
   loglog(dt, err(:,j), 'Displayname', nom_param{j},'linestyle', 'none', 'marker', '*');
   hold on
   grid on
   loglog(xfit,yfit(:,j), 'k--', 'Displayname', legende{j});
   xlabel('$\Delta t$ [s]');
   ylabel('$\epsilon$');
   legend('location', 'best','NumColumns',2);
end
set(fig{1},'Position',[100 100 900 426]);


%% Graphes de quantités conservées (avoir run l'étude de convergence !) = %
%=========================================================================%
nomFigure{2} = '8_2i_proba';
nomFigure{3} = '8_2i_err_P_E';
nomFigure{4} = '8_2i_Heisenberg';
%=========================================================================%
%=========================================================================%

% Charger les données
data = load(nomfile{end,3});
t = data(:,1);
Pgauche = data(:,2);
Pdroite = data(:,3);
Ptot = data(:,4);
E = data(:,5);
DeltaxDeltap = data(:,10); %(Delta x)(Delta p)
Deltax = data(:,11);
Deltap = data(:,12);

% Proba
fig{2} = figure('Name',nomFigure{2},'NumberTitle','off');
hold on; grid on;
plot(t,Pgauche, 'Displayname', '$P_{x<0}(t)$');
plot(t,Pdroite, 'Displayname', '$P_{x>0}(t)$');
plot(t,Ptot, 'Displayname', '$P_{tot}(t)$');
xlabel('$t$ [s]');
ylabel('$P$');
legend('location', 'best');

% Erreur sur Proba total et Energie
fig{3} = figure('Name',nomFigure{3},'NumberTitle','off');
semilogy(t,abs(Ptot - 1), 'Displayname', 'Probability');
hold on; grid on;
semilogy(t,abs(E-E(1))./E(1), 'Displayname', 'Energy');
xlabel('$t$ [s]');
ylabel('$\epsilon$');
legend('location', 'best');

% incertitude heisenberg
fig{4} = figure('Name',nomFigure{4},'NumberTitle','off');
plot(t, DeltaxDeltap);
hold on; grid on;
yline(0.5, 'k--','$\frac{\hbar}{2}$', 'fontsize', 20,'linewidth', 2,'LabelHorizontalAlignment', 'left','LabelVerticalAlignment','bottom', 'Interpreter', 'latex');
xlabel('$t$ [s]');
ylabel('$\langle \Delta x \rangle \cdot \langle \Delta p \rangle$');
ylim([0.4 0.9]);

%% Save figures ========================================================= %
saveFigure = 0;         %1 = TRUE, 0 = FALSE
sauver = 0;%[1 2 3 4 5];

if saveFigure == true
    for i=sauver
        exportgraphics(fig{i},[nomFigure{i}, '.pdf'],'ContentType','vector');
    end
end


%% Fonctions utilisées



