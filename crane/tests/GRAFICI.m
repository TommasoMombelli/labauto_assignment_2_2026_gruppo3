% Assicuriamoci che 'time' sia un vettore colonna (Nx1) come gli altri dati,
% dato che nel tuo workspace è caricato come riga (1xN).
close all; clear all; clc; 
% load("C:\Users\sbarr\LabAuto\labauto_assignment_2_2026_gruppo3\crane\tests\test_trj1_shaped_ZVDD_20260520114420.mat")
data=load ("test_trj1_shaped_ZVDD_20260609002348.mat")
data_Base=load("test_trj1_baseline")
t = data.time(:);
t_ref=data_Base.time(:);
name=data.name; 

%% FIGURA 1: Tracking (Posizione, Velocità, Sforzo di Controllo)
figure('Name', ['Tracking - ' name], 'NumberTitle', 'off');

% 1. Posizione (Reale vs Riferimento)
subplot(3, 1, 1);
plot(t_ref, data_Base.reference_position, '--r', 'LineWidth', 1.5);
hold on;
plot(t, data.joint_position, 'b', 'LineWidth', 1.5);
grid on;
title('Inseguimento di Posizione');
xlabel('Tempo [s]');
ylabel('Posizione');
legend('Reference', 'Actual', 'Location', 'best');

% 2. Velocità (Reale vs Riferimento)
subplot(3, 1, 2);
plot(t_ref, data_Base.reference_velocity, '--r', 'LineWidth', 1.5);
hold on;
plot(t, data.joint_velocity, 'b', 'LineWidth', 1.5);
grid on;
title('Inseguimento di Velocità');
xlabel('Tempo [s]');
ylabel('Velocità');
legend('Reference', 'Actual', 'Location', 'best');

% 3. Sforzo di Controllo (Joint Torque)
subplot(3, 1, 3);
plot(t, data.joint_torque, 'k', 'LineWidth', 1.5);
grid on;
title('Sforzo di Controllo (Actuator Force)');
xlabel('Tempo [s]');
ylabel('Coppia / Forza');
legend('Control Action', 'Location', 'best');

%% FIGURA 2: Analisi degli Errori
figure('Name', ['Errori - ' name], 'NumberTitle', 'off');

% Calcolo degli errori
error_pos = data_Base.reference_position - data.joint_position;
error_vel = data_Base.reference_velocity - data.joint_velocity;

% Calcolo MAE (Mean Absolute Error)
mae_pos = mean(abs(error_pos));
mae_vel = mean(abs(error_vel));

% 1. Errore di Posizione
subplot(2, 1, 1);
plot(t, error_pos, 'b', 'LineWidth', 1.5);
grid on;
title(sprintf('Errore di Posizione (MAE = %.4f)', mae_pos));
xlabel('Tempo [s]');
ylabel('Errore Posizione');

% 2. Errore di Velocità
subplot(2, 1, 2);
plot(t, error_vel, 'r', 'LineWidth', 1.5);
grid on;
title(sprintf('Errore di Velocità (MAE = %.4f)', mae_vel));
xlabel('Tempo [s]');
ylabel('Errore Velocità');

% =========================================================================
% IMPOSTAZIONI SALVATAGGIO
% =========================================================================
suffix = '_ZVDD'; % <--- CAMBIA QUESTO SUFFISSO OGNI VOLTA 
% =========================================================================

% Assicuriamoci che 'time' sia un vettore colonna (Nx1)
% t = data.time(:); 

%% FIGURA 1: Tracking (Posizione, Velocità, Sforzo di Controllo)
% Assegniamo la figura alla variabile 'fig1'
fig1 = figure('Name', ['Tracking - ' name], 'NumberTitle', 'off');

% 1. Posizione (Reale vs Riferimento)
subplot(3, 1, 1);
plot(t_ref, data_Base.reference_position, '--r', 'LineWidth', 1.5);
hold on;
plot(t, data.joint_position, 'b', 'LineWidth', 1.5);
grid on;
title('Inseguimento di Posizione');
xlabel('Tempo [s]');
ylabel('Posizione');
legend('Reference', 'Actual', 'Location', 'best');

% 2. Velocità (Reale vs Riferimento)
subplot(3, 1, 2);
plot(t_ref, data_Base.reference_velocity, '--r', 'LineWidth', 1.5);
hold on;
plot(t, data.joint_velocity, 'b', 'LineWidth', 1.5);
grid on;
title('Inseguimento di Velocità');
xlabel('Tempo [s]');
ylabel('Velocità');
legend('Reference', 'Actual', 'Location', 'best');

% 3. Sforzo di Controllo (Joint Torque)
subplot(3, 1, 3);
plot(t, data.joint_torque, 'k', 'LineWidth', 1.5);
grid on;
title('Sforzo di Controllo (Actuator Force)');
xlabel('Tempo [s]');
ylabel('Coppia / Forza');
legend('Control Action', 'Location', 'best');

%% FIGURA 2: Analisi degli Errori
% Assegniamo la figura alla variabile 'fig2'
fig2 = figure('Name', ['Errori - ' name], 'NumberTitle', 'off');

% Calcolo degli errori
error_pos = data_Base.reference_position - data.joint_position;
error_vel = data_Base.reference_velocity - data.joint_velocity;

% Calcolo MAE (Mean Absolute Error)
mae_pos = mean(abs(error_pos));
mae_vel = mean(abs(error_vel));

% 1. Errore di Posizione
subplot(2, 1, 1);
plot(t, error_pos, 'b', 'LineWidth', 1.5);
grid on;
title(sprintf('Errore di Posizione (MAE = %.4f)', mae_pos));
xlabel('Tempo [s]');
ylabel('Errore Posizione');

% 2. Errore di Velocità
subplot(2, 1, 2);
plot(t, error_vel, 'r', 'LineWidth', 1.5);
grid on;
title(sprintf('Errore di Velocità (MAE = %.4f)', mae_vel));
xlabel('Tempo [s]');
ylabel('Errore Velocità');

%% =========================================================================
% SALVATAGGIO DELLE IMMAGINI PNG
% =========================================================================

% Creazione dinamica dei nomi dei file
filename_fig1 = ['Tracking', suffix, '.png'];
filename_fig2 = ['Errori', suffix, '.png'];

% Salvataggio ad alta risoluzione (300 DPI) per il report in LaTeX
exportgraphics(fig1, filename_fig1, 'Resolution', 300);
exportgraphics(fig2, filename_fig2, 'Resolution', 300);

% Messaggio a schermo di conferma
fprintf('\nImmagini salvate con successo nella cartella corrente:\n - %s\n - %s\n', filename_fig1, filename_fig2);


% =========================================================================
% CALCOLO DEI TEMPI DI FINE MOTO E ASSESTAMENTO
% =========================================================================

% 1. Calcolo del tempo di fine riferimento (quando la velocità del ref va a zero)
% Cerchiamo l'ultimo indice in cui la velocità di riferimento è significativa
idx_end_ref = find(abs(data_Base.reference_velocity) > 1e-4, 1, 'last');
if isempty(idx_end_ref)
    t_end_ref = t_ref(end); % Usiamo t_ref per i dati di base
else
    t_end_ref = t_ref(idx_end_ref);
end

% 2. Calcolo del tempo di assestamento (Settling Time al 5%)
% Spostamento totale (AGGIUNTO data_Base.)
spostamento_totale = abs(data_Base.reference_position(end) - data_Base.reference_position(1));
banda_tolleranza = 0.05 * spostamento_totale; % Tolleranza del 5%

% L'errore rispetto alla posizione FINALE (AGGIUNTO data. e data_Base.)
errore_da_target = abs(data.joint_position - data_Base.reference_position(end));

% Cerchiamo l'ultimo istante in cui l'errore è FUORI dalla banda di tolleranza
idx_settling = find(errore_da_target > banda_tolleranza, 1, 'last');

if isempty(idx_settling)
    t_settling = t(end);
else
    t_settling = t(idx_settling);
end

% Stampiamo i risultati in console
fprintf('\n--- ANALISI DEI TEMPI ---\n');
fprintf('Tempo di fine riferimento (Comando): %.3f s\n', t_end_ref);
fprintf('Tempo di assestamento (Reale al 5%%): %.3f s\n', t_settling);