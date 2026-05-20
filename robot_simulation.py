import os
import time
import yaml
from collections import deque

from scipy.io import savemat, loadmat
from scipy.signal import find_peaks
from datetime import datetime

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from labauto import MuJoCoMechanicalSystem
from labauto import TrapezoidalMotionLaw
from labauto import loadController
from labauto import loadInstructions

# ==============================================================================
# CONFIGURAZIONE FLAG DI SIMULAZIONE
# ==============================================================================
USE_SHAPER = True  # Cambia a True per usare lo Shaper, False per la Baseline

# Inserisci qui il nome ESATTO del file salvato senza shaper (con estensione .mat)
MANUAL_BASELINE_FILE = "test_trj1_baseline.mat"

# ==============================================================================
# FUNZIONE DI ESTRAZIONE PARAMETRI DAI DATI (T, zeta)
# ==============================================================================
def extract_tz_from_data(mat_file_path, Tc):
    try:
        print(f"Estrazione parametri dal file: {mat_file_path}")
        data = loadmat(mat_file_path)
        t = data["time"].flatten()
        # Calcoliamo l'errore di posizione (asse x)
        error_x = (data["reference_position"] - data["joint_position"])[:, 0]
        
        # Analizziamo solo l'ultima metà della simulazione (vibrazione residua libera)
        tail_idx = len(t) // 2
        error_x_tail = error_x[tail_idx:]
        t_tail = t[tail_idx:]
        
        # Cerca i picchi PRINCIPALI (distanti almeno 1.0 sec per evitare falsi positivi)
        min_dist = int(1.0 / Tc)
        
        # 'prominence' forza l'algoritmo a ignorare le microsospensioni e il rumore numerico.
        peaks, _ = find_peaks(error_x_tail, distance=min_dist, prominence=0.01)
        
        if len(peaks) >= 2:
            # 1. Calcolo del Periodo T
            periodi = np.diff(t_tail[peaks])
            T_esatto = np.mean(periodi)
            
            # 2. Calcolo dello Smorzamento zeta (Decremento Logaritmico)
            n = len(peaks) - 1
            x0 = error_x_tail[peaks[0]]
            xn = error_x_tail[peaks[-1]]
            
            if x0 > 0 and xn > 0 and x0 > xn:
                delta = (1 / n) * np.log(x0 / xn)
                zeta_esatto = delta / np.sqrt((2 * np.pi)**2 + delta**2)
            else:
                zeta_esatto = 0.0 # Ampiezza costante o in crescita = smorzamento nullo
                
            print(f"--> Parametri estratti: T = {T_esatto:.4f} s, zeta = {zeta_esatto:.6f}")
            return T_esatto, zeta_esatto
        else:
            print("Non ci sono abbastanza picchi principali nella vibrazione residua. Uso default.")
            return 1.5, 0.0
            
    except Exception as e:
        print(f"Errore durante l'analisi dei dati ({e}). Uso default.")
        return 1.5, 0.0

# ==============================================================================
# CLASSE INPUT SHAPER (ZV - Zero Vibration con Ritardo Frazionario)
# ==============================================================================
class ZVShaper:
    def __init__(self, omega_n, zeta, Tc):
        # Calcolo della frequenza smorzata
        omega_d = omega_n * np.sqrt(1 - zeta**2) if zeta < 1.0 else omega_n
        
        # Calcolo del decadimento K
        K = np.exp(-(zeta * np.pi) / np.sqrt(1 - zeta**2)) if zeta < 1.0 else 1.0
        
        # Calcolo ampiezze ideali degli impulsi
        A1 = 1 / (1 + K)
        A2 = K / (1 + K)
        
        # --- IMPLEMENTAZIONE RITARDO FRAZIONARIO ---
        self.delays = []
        self.weights = []
        
        # 1. Primo impulso (t=0) cade esattamente al campione 0
        self.delays.append(0)
        self.weights.append(A1)
        
        # 2. Secondo impulso (t = pi / omega_d)
        t2 = np.pi / omega_d
        n_esatto = t2 / Tc  # Es: 50.3
        
        N2 = int(np.floor(n_esatto))  # Parte intera (Es: 50)
        alpha = n_esatto - N2         # Parte frazionaria (Es: 0.3)
        
        # Distribuiamo l'ampiezza A2 su due campioni adiacenti
        self.delays.append(N2)
        self.weights.append(A2 * (1 - alpha))  # Quota per il campione N2
        
        self.delays.append(N2 + 1)
        self.weights.append(A2 * alpha)        # Quota per il campione N2+1
        
        # Inizializzazione del buffer storico per i riferimenti
        self.buffer_size = max(self.delays) + 1
        self.history = deque([(0.0, 0.0, 0.0)] * self.buffer_size, maxlen=self.buffer_size)
        self.initialized = False

    def shape(self, q, dq, ddq):
        # Inizializza il buffer con il primo valore reale per evitare "salti" all'avvio
        if not self.initialized:
            self.history = deque([(q, dq, ddq)] * self.buffer_size, maxlen=self.buffer_size)
            self.initialized = True
            
        # Aggiunge i nuovi riferimenti in testa (sinistra) al buffer
        self.history.appendleft((q, dq, ddq))
        
        # Calcola i valori "shapati" combinando pesi e ritardi
        shaped_q = sum(w * self.history[d][0] for w, d in zip(self.weights, self.delays))
        shaped_dq = sum(w * self.history[d][1] for w, d in zip(self.weights, self.delays))
        shaped_ddq = sum(w * self.history[d][2] for w, d in zip(self.weights, self.delays))
        
        return shaped_q, shaped_dq, shaped_ddq

# ==============================================================================
# CONFIGURAZIONE E AVVIO SIMULAZIONE
# ==============================================================================
model_name = "crane"  # folder containing model.xml + control_config.yaml + motion program
program_name = "test_trj1"

with open(f'{model_name}/control_config.yaml', 'r') as file:
    params_yaml = yaml.safe_load(file)
    controller_params = params_yaml['controller']
    dynamic_params = np.array(params_yaml['model_parameters'])

xml_path = f"{model_name}/model.xml"
robot = MuJoCoMechanicalSystem(xml_path=xml_path, motor_actuators=["motor_1"], motor_joints=["joint_1"], spring_joints=[],ee_site="payload")
robot = MuJoCoMechanicalSystem(xml_path=xml_path)
robot.show()
robot.initialize()
dof = robot.get_input_number()

Tc = robot.get_sampling_period()

# ------------------------------------------------------------------------------
# INIZIALIZZAZIONE DELLO SHAPER ESTRAENDO I DATI (SE ATTIVATO)
# ------------------------------------------------------------------------------
shaper = None
if USE_SHAPER:
    mat_baseline_path = f"{model_name}/tests/{MANUAL_BASELINE_FILE}"

    if os.path.exists(mat_baseline_path):
        T_oscillazione, zeta = extract_tz_from_data(mat_baseline_path, Tc)
        omega_n = (2 * np.pi) / T_oscillazione
        shaper = ZVShaper(omega_n, zeta, Tc)
        print(f"[STATUS] Input Shaper ATTIVO: omega_n={omega_n:.3f} rad/s, zeta={zeta:.5f}")
    else:
        print(f"[ERRORE] File baseline non trovato in {mat_baseline_path}. Esecuzione SENZA shaper.")
        USE_SHAPER = False
else:
    print("[STATUS] Esecuzione in modalità BASELINE (Shaper Disattivato)")
# ------------------------------------------------------------------------------

decentralized_ctrl=loadController(Tc,controller_params,dynamic_params,model_name)
decentralized_ctrl.initialize()
decentralized_ctrl.set_umax(robot.get_umax())

measured_output = robot.read_sensor_value()
q0 = measured_output[:dof]
Dq0 = measured_output[dof:]
DDq0 = np.zeros(dof)
initial_reference = np.concatenate((q0, Dq0, DDq0))

max_Dq = np.array([5.5]*dof)
max_DDq = np.array([5.0]*dof)
motion_law_params=dict()
motion_law_params['max_velocity']=max_Dq
motion_law_params['max_acceleration']=max_DDq
ml = TrapezoidalMotionLaw(motion_law_params, Tc)
ml.set_initial_condition(q0)

instructions = loadInstructions(f'{model_name}/{program_name}.txt')
ml.add_instructions(instructions)

joint_torque = robot.read_actuator_value()
feedforward_action = np.array([0.0]*dof)

decentralized_ctrl.starting(initial_reference, measured_output, joint_torque, feedforward_action)

# Simulation loop
t, measured_signal, control_action, reference_signal, link_position = [], [], [], [], []
actual_time = 0.0

while ml.depending_instructions():
    loop_t0 = time.perf_counter()
    target_q, target_Dq, target_DDq = ml.compute_motion_law()

    target_q_is = target_q[0]
    target_Dq_is = target_Dq[0]
    target_DDq_is = target_DDq[0]

    # --- APPLICAZIONE CONDIZIONALE DELL'INPUT SHAPER ---
    if USE_SHAPER and shaper is not None:
        shaped_q, shaped_dq, shaped_ddq = shaper.shape(target_q_is, target_Dq_is, target_DDq_is)
        reference = np.array([shaped_q, shaped_dq, shaped_ddq])
    else:
        reference = np.array([target_q_is, target_Dq_is, target_DDq_is])

    measured_output = robot.read_sensor_value()

    joint_torque = decentralized_ctrl.compute_control_action(reference, measured_output, feedforward_action)
    robot.write_actuator_value(joint_torque)

    t.append(actual_time)
    measured_signal.append(measured_output)
    control_action.append(joint_torque)
    reference_signal.append(reference)
    link_position.append(robot.link_position())
    actual_time += Tc

    robot.simulate()

    computation_time = time.perf_counter() - loop_t0
    time.sleep(max(0.0, Tc - computation_time))

# ==============================================================================
# POST-PROCESSING E SALVATAGGIO DATI
# ==============================================================================
t = np.array(t)
measured_signal = np.array(measured_signal)
control_action = np.array(control_action)
reference_signal = np.array(reference_signal)
link_position = np.array(link_position)

joint_position = measured_signal[:, :dof]
joint_velocity = measured_signal[:, dof:]

reference_position = reference_signal[:, :dof]
reference_velocity = reference_signal[:, dof:2*dof]
reference_acceleration = reference_signal[:, 2*dof:]

timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

# Cambia il suffisso in base alla modalità
mode_suffix = "shaped" if USE_SHAPER else "baseline"
filename = f"{model_name}/tests/{program_name}_{mode_suffix}_{timestamp}.mat" if USE_SHAPER else f"{model_name}/tests/{program_name}_{mode_suffix}.mat"
test_data = {
    "reference_position": reference_position,
    "reference_velocity": reference_velocity,
    "reference_acceleration": reference_acceleration,
    "joint_position": joint_position,
    "joint_velocity": joint_velocity,
    "joint_torque": control_action,
    "link_position": link_position,
    "time": t,
    "name": f"test_{mode_suffix}"
}
savemat(filename, {key: test_data[key] for key in test_data})
print(f"[SALVATAGGIO] File salvato: {filename}")

# ==============================================================================
# PLOTTING DEI RISULTATI
# ==============================================================================
labels = ["x"]

fig1 = make_subplots(
    rows=3, cols=3,
    shared_xaxes=True,
    subplot_titles=[f"Position {a}" for a in labels] +
                   [f"Velocity {a}" for a in labels] +
                   [f"Actuator force {a} (motor-side)" for a in labels]
)

for i, a in enumerate(labels):
    col = i + 1
    fig1.add_trace(go.Scatter(x=t, y=joint_position[:, i], name=f"q_{a}", legendgroup=f"pos_{a}"), row=1, col=col)
    fig1.add_trace(go.Scatter(x=t, y=reference_position[:, i], name=f"qref_{mode_suffix}_{a}", legendgroup=f"pos_{a}", line=dict(dash="dash")), row=1, col=col)
    
    fig1.add_trace(go.Scatter(x=t, y=joint_velocity[:, i], name=f"dq_{a}", legendgroup=f"vel_{a}"), row=2, col=col)
    fig1.add_trace(go.Scatter(x=t, y=reference_velocity[:, i], name=f"dqref_{mode_suffix}_{a}", legendgroup=f"vel_{a}", line=dict(dash="dash")), row=2, col=col)
    
    fig1.add_trace(go.Scatter(x=t, y=control_action[:, i], name=f"F_{a}", legendgroup=f"u_{a}"), row=3, col=col)

fig1.update_xaxes(title_text="Time (s)", row=3, col=2)
fig1.update_layout(title=f"Tracking ({mode_suffix.upper()}): position / velocity / control", height=900, width=1200, legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0))
fig1.update_xaxes(showgrid=True)
fig1.update_yaxes(showgrid=True)

position_error = reference_position - joint_position
velocity_error = reference_velocity - joint_velocity
mae_pos = np.mean(np.abs(position_error), axis=0)
mae_vel = np.mean(np.abs(velocity_error), axis=0)

subplot_titles = (
    [f"Position error {labels[i]} (MAE={mae_pos[i]:4.3f})" for i in range(dof)] +
    [f"Velocity error {labels[i]} (MAE={mae_vel[i]:4.3f})" for i in range(dof)] +
    [f"Actuator force {labels[i]} (motor-side)" for i in range(dof)]
)

fig2 = make_subplots(rows=3, cols=3, shared_xaxes=True, subplot_titles=subplot_titles)

for i, a in enumerate(labels):
    col = i + 1
    fig2.add_trace(go.Scatter(x=t, y=position_error[:, i], name=f"e_q_{a}", legendgroup=f"ep_{a}"), row=1, col=col)
    fig2.add_trace(go.Scatter(x=t, y=velocity_error[:, i], name=f"e_dq_{a}", legendgroup=f"ev_{a}"), row=2, col=col)
    fig2.add_trace(go.Scatter(x=t, y=control_action[:, i], name=f"F_{a}", legendgroup=f"u_{a}"), row=3, col=col)

fig2.update_xaxes(title_text="Time (s)", row=3, col=2)
fig2.update_layout(title=f"Errors ({mode_suffix.upper()}): position error / velocity error / control", height=900, width=1200, legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0))
fig2.update_xaxes(showgrid=True)
fig2.update_yaxes(showgrid=True)

fig1.show()
fig2.show()
robot.close()