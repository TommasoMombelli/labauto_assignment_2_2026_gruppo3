
% Assicuriamoci che 'time' sia un vettore colonna (Nx1) come gli altri dati,
% dato che nel tuo workspace è caricato come riga (1xN).
close all; clear all; %clc; 

% Ottenimento dinamico della cartella in cui si trova lo script
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd; % Se eseguito riga per riga, usa la directory corrente
end

% load("C:\Users\sbarr\LabAuto\labauto_assignment_2_2026_gruppo3\crane\tests\test_trj1_shaped_ZVDD_20260520114420.mat")
data = load("test_trj1_shaped_EI_20260609002518.mat");
data_Base = load("test_trj1_baseline");

t = data.time(:);
t_ref = data_Base.time(:);
name = data.name; 

% =========================================================================
% IMPOSTAZIONI SALVATAGGIO
% =========================================================================
suffix = '_EI'; % <--- CAMBIA QUESTO SUFFISSO OGNI VOLTA 

%% FIGURA 1: Tracking (Posizione, Velocità, Sforzo di Controllo)
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
filename_fig1 = fullfile(script_dir, ['Tracking', suffix, '.png']);
filename_fig2 = fullfile(script_dir, ['Errori', suffix, '.png']);

exportgraphics(fig1, filename_fig1, 'Resolution', 300);
exportgraphics(fig2, filename_fig2, 'Resolution', 300);

fprintf('\nImmagini salvate con successo nella cartella dello script:\n - Tracking%s.png\n - Errori%s.png\n', suffix, suffix);

%% =========================================================================
% CALCOLO DEI TEMPI DI FINE MOTO E ASSESTAMENTO
% =========================================================================

% 1. Calcolo del tempo finale della Baseline 
idx_end_base = find(abs(data_Base.reference_velocity) > 1e-4, 1, 'last');
if isempty(idx_end_base)
    t_end_base = t_ref(end); 
else
    t_end_base = t_ref(idx_end_base);
end

% 2. Calcolo del tempo finale del Filtro
% Cerchiamo prima se c'è la reference_velocity salvata per il filtro, altrimenti usiamo la velocità reale
if isfield(data, 'reference_velocity')
    idx_end_filter = find(abs(data.reference_velocity) > 1e-4, 1, 'last');
else
    idx_end_filter = find(abs(data.joint_velocity) > 1e-4, 1, 'last');
end

if isempty(idx_end_filter)
    t_end_filter = t(end);
else
    t_end_filter = t(idx_end_filter);
end

% 3. Calcolo DT: Tempo Finale Filtro - Tempo Finale Baseline - 
dt = t_end_filter - t_end_base;

% 4. Calcolo del tempo di assestamento (Settling Time al 5%) PER IL FILTRO
spostamento_totale = abs(data_Base.reference_position(end) - data_Base.reference_position(1));
banda_tolleranza = 0.05 * spostamento_totale; % Tolleranza del 5%

errore_da_target = abs(data.joint_position - data_Base.reference_position(end));
idx_settling = find(errore_da_target > banda_tolleranza, 1, 'last');

if isempty(idx_settling)
    t_settling = t(end);
else
    t_settling = t(idx_settling);
end

% Stampiamo i risultati in console
fprintf('\n--- ANALISI DEI TEMPI ---\n');
fprintf('Tempo finale Baseline              : %.3f s\n', t_end_base);
fprintf('Tempo finale Filtro                : %.3f s\n', t_end_filter);
fprintf('dt (Baseline - Filtro)             : %.3f s\n', dt);
fprintf('Tempo di assestamento (Filtro 5%%)  : %.3f s\n', t_settling);

%% =========================================================================
% CALCOLO SFORZO ED ENERGIA MECCANICA
% =========================================================================
picco_max_posizione = max(abs(data.joint_position));  
picco_max_sforzo = max(abs(data.joint_torque));       

sforzo_quadratico = trapz(t, data.joint_torque.^2);   

potenza_meccanica = abs(data.joint_torque .* data.joint_velocity);
energia_joule = trapz(t, potenza_meccanica);

fprintf('\n--- ANALISI ENERGETICA E PRESTAZIONI ---\n');
fprintf('Picco di Posizione     : %.4f\n', picco_max_posizione);
fprintf('Picco di Sforzo        : %.3f Nm\n', picco_max_sforzo);
fprintf('Sforzo Quadratico Medio: %.2f (Nm)^2*s\n', sforzo_quadratico);
fprintf('Energia Meccanica Pura : %.2f Joule\n', energia_joule);

%% =========================================================================
% ESPORTAZIONE IN EXCEL
% =========================================================================
nome_filtro = strrep(suffix, '_', ''); 
nome_file_excel = fullfile(script_dir, ['Risultati_Ottimali', suffix, '.xlsx']);

% Compilazione della tabella con il nuovo ordine esatto
NuoviDati = table({nome_filtro}, mae_pos, mae_vel, t_end_base, t_end_filter, dt, t_settling, ...
    picco_max_posizione, picco_max_sforzo, sforzo_quadratico, energia_joule, ...
    'VariableNames', {'Filtro', 'MAE_Posizione', 'MAE_Velocita', 'Tempo_Finale_Baseline_s', 'Tempo_Finale_Filtro_s', ...
                      'dt_s', 'Tempo_Assestamento_5perc_s', 'Picco_Max_Posizione', 'Picco_Max_Sforzo_Nm', ...
                      'Sforzo_Quadratico', 'Energia_Joule'});

if isfile(nome_file_excel)
    writetable(NuoviDati, nome_file_excel, 'WriteMode', 'Append', 'WriteVariableNames', false);
    fprintf('\n-> Dati AGGIUNTI al file Excel: %s\n\n', nome_file_excel);
else
    writetable(NuoviDati, nome_file_excel, 'WriteMode', 'Overwrite');
    fprintf('\n-> File Excel CREATO. Dati inseriti: %s\n\n', nome_file_excel);
end
