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
# Scegli la modalità: 
# "IDENTIFY" -> Usa traiettoria semplice, calcola e salva i parametri.
# "NONE"     -> Usa traiettoria complessa SENZA shaper (crea la baseline).
# "ZV", "ZVD", "ZVDD", "EI" -> Usa traiettoria complessa CON shaper (caricando i parametri).
SHAPER_TYPE = "NONE" 

model_name = "crane"  

# Scelta automatica del programma (traiettoria)
if SHAPER_TYPE == "IDENTIFY":
    program_name = "identification_trj1"
else:
    program_name = "test_trj1"

print(f"[INIT] Modalità selezionata: {SHAPER_TYPE}")
print(f"[INIT] Traiettoria in uso: {program_name}.txt")

# Inserisci qui il nome ESATTO del file salvato senza shaper (con estensione .mat)
MANUAL_BASELINE_FILE = "test_trj1_baseline.mat"


# ==============================================================================
# FUNZIONE DI ESTRAZIONE PARAMETRI DAI DATI (T, zeta)
# ==============================================================================
def extract_tz_from_data(mat_file_path, Tc):
    try:
        print(f"\n--- Estrazione parametri dal file: {mat_file_path} ---")
        data = loadmat(mat_file_path)
        t = data["time"].flatten()
        
        # Calcoliamo l'errore di posizione (asse x)
        error_x = (data["reference_position"] - data["joint_position"])[:, 0]
        
        # Analizziamo la vibrazione libera partendo da t = 20.0 secondi
        start_idx = np.searchsorted(t, 20.0)
        error_x_tail = error_x[start_idx:]
        t_tail = t[start_idx:]
        
        # Cerca i picchi PRINCIPALI (distanti almeno 1.0 sec per evitare falsi positivi)
        min_dist = int(1.0 / Tc)
        
        # 'prominence' forza l'algoritmo a ignorare le microsospensioni e il rumore numerico.
        peaks, _ = find_peaks(error_x_tail, distance=min_dist, prominence=0.01)
        t_peaks = t_tail[peaks]
        
        # --- NOVITÀ: Generazione codice MATLAB da copiare e incollare ---
        if len(peaks) > 0:
            valori_picchi = error_x_tail[peaks]
            
            # Creiamo stringhe formattate in stile MATLAB come VETTORI RIGA [val1, val2, val3, ...]
            str_t_peaks = "[" + ", ".join([f"{val:.4f}" for val in t_peaks]) + "]"
            str_valori = "[" + ", ".join([f"{val:.6f}" for val in valori_picchi]) + "]"
            
            print("\n=== COPIA E INCOLLA IN MATLAB ===")
            print("% Tempi in cui si verificano i picchi (vettore riga)")
            print(f"peaks_times = {str_t_peaks};")
            print("")
            print("% Ampiezza dell'errore di posizione in quei picchi (vettore riga)")
            print(f"valori_picchi = {str_valori};")
            print("=================================\n")
        # ----------------------------------------------------------------

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
                
            print(f"--> Parametri estratti: T = {T_esatto:.4f} s, zeta = {zeta_esatto:.6f}\n")
            return T_esatto, zeta_esatto
        else:
            print("Non ci sono abbastanza picchi principali nella vibrazione residua. Uso default.")
            return 1.5, 0.0
            
    except Exception as e:
        print(f"Errore durante l'analisi dei dati ({e}). Uso default.")
        return 1.5, 0.0

# ==============================================================================
# CLASSE INPUT SHAPER (ZV - Zero Vibration)
# ==============================================================================
class ZVShaper:
    def __init__(self, omega_n, zeta, Tc):
        omega_d = omega_n * np.sqrt(1 - zeta**2) if zeta < 1.0 else omega_n
        K = np.exp(-(zeta * np.pi) / np.sqrt(1 - zeta**2)) if zeta < 1.0 else 1.0
        
        A1 = 1 / (1 + K)
        A2 = K / (1 + K)
        
        self.delays = []
        self.weights = []
        
        # 1. Primo impulso
        self.delays.append(0)
        self.weights.append(A1)
        
        # 2. Secondo impulso
        t2 = np.pi / omega_d
        n2_esatto = t2 / Tc 
        N2 = int(np.floor(n2_esatto)) 
        alpha2 = n2_esatto - N2 
        
        self.delays.append(N2)
        self.weights.append(A2 * (1 - alpha2)) 
        self.delays.append(N2 + 1)
        self.weights.append(A2 * alpha2) 
        
        self.buffer_size = max(self.delays) + 1
        self.history = deque([(0.0, 0.0, 0.0)] * self.buffer_size, maxlen=self.buffer_size)
        self.initialized = False

    def shape(self, q, dq, ddq):
        if not self.initialized:
            self.history = deque([(q, dq, ddq)] * self.buffer_size, maxlen=self.buffer_size)
            self.initialized = True
            
        self.history.appendleft((q, dq, ddq))
        
        shaped_q = sum(w * self.history[d][0] for w, d in zip(self.weights, self.delays))
        shaped_dq = sum(w * self.history[d][1] for w, d in zip(self.weights, self.delays))
        shaped_ddq = sum(w * self.history[d][2] for w, d in zip(self.weights, self.delays))
        
        return shaped_q, shaped_dq, shaped_ddq

# ==============================================================================
# CLASSE INPUT SHAPER (ZVD - Zero Vibration and Derivative)
# ==============================================================================
class ZVDShaper:
    def __init__(self, omega_n, zeta, Tc):
        omega_d = omega_n * np.sqrt(1 - zeta**2) if zeta < 1.0 else omega_n
        K = np.exp(-(zeta * np.pi) / np.sqrt(1 - zeta**2)) if zeta < 1.0 else 1.0
        
        A1 = 1 / ((1 + K)**2)
        A2 = (2 * K) / ((1 + K)**2)
        A3 = (K**2) / ((1 + K)**2)
        
        self.delays = []
        self.weights = []
        
        # 1. Primo impulso
        self.delays.append(0)
        self.weights.append(A1)
        
        # 2. Secondo impulso (t = pi / omega_d)
        t2 = np.pi / omega_d
        n2_esatto = t2 / Tc 
        N2 = int(np.floor(n2_esatto))
        alpha2 = n2_esatto - N2 
        
        self.delays.append(N2)
        self.weights.append(A2 * (1 - alpha2)) 
        self.delays.append(N2 + 1)
        self.weights.append(A2 * alpha2)
        
        # 3. Terzo impulso (t = 2*pi / omega_d)
        t3 = 2 * np.pi / omega_d
        n3_esatto = t3 / Tc
        N3 = int(np.floor(n3_esatto))
        alpha3 = n3_esatto - N3
        
        self.delays.append(N3)
        self.weights.append(A3 * (1 - alpha3))
        self.delays.append(N3 + 1)
        self.weights.append(A3 * alpha3)
        
        self.buffer_size = max(self.delays) + 1
        self.history = deque([(0.0, 0.0, 0.0)] * self.buffer_size, maxlen=self.buffer_size)
        self.initialized = False

    def shape(self, q, dq, ddq):
        if not self.initialized:
            self.history = deque([(q, dq, ddq)] * self.buffer_size, maxlen=self.buffer_size)
            self.initialized = True
            
        self.history.appendleft((q, dq, ddq))
        
        shaped_q = sum(w * self.history[d][0] for w, d in zip(self.weights, self.delays))
        shaped_dq = sum(w * self.history[d][1] for w, d in zip(self.weights, self.delays))
        shaped_ddq = sum(w * self.history[d][2] for w, d in zip(self.weights, self.delays))
        
        return shaped_q, shaped_dq, shaped_ddq

# ==============================================================================
# CLASSE INPUT SHAPER (ZVDD - Zero Vibration, Derivative & Second Derivative)
# ==============================================================================
class ZVDDShaper:
    def __init__(self, omega_n, zeta, Tc):
        omega_d = omega_n * np.sqrt(1 - zeta**2) if zeta < 1.0 else omega_n
        K = np.exp(-(zeta * np.pi) / np.sqrt(1 - zeta**2)) if zeta < 1.0 else 1.0
        
        # Calcolo ampiezze ideali degli impulsi ZVDD (pattern binomiale)
        Denom = (1 + K)**3
        A1 = 1 / Denom
        A2 = (3 * K) / Denom
        A3 = (3 * (K**2)) / Denom
        A4 = (K**3) / Denom
        
        self.delays = []
        self.weights = []
        
        # 1. Primo impulso
        self.delays.append(0)
        self.weights.append(A1)
        
        # 2. Secondo impulso
        t2 = np.pi / omega_d
        n2_esatto = t2 / Tc 
        N2 = int(np.floor(n2_esatto))
        alpha2 = n2_esatto - N2 
        
        self.delays.append(N2)
        self.weights.append(A2 * (1 - alpha2)) 
        self.delays.append(N2 + 1)
        self.weights.append(A2 * alpha2)
        
        # 3. Terzo impulso
        t3 = 2 * np.pi / omega_d
        n3_esatto = t3 / Tc
        N3 = int(np.floor(n3_esatto))
        alpha3 = n3_esatto - N3
        
        self.delays.append(N3)
        self.weights.append(A3 * (1 - alpha3))
        self.delays.append(N3 + 1)
        self.weights.append(A3 * alpha3)

        # 4. Quarto impulso
        t4 = 3 * np.pi / omega_d
        n4_esatto = t4 / Tc
        N4 = int(np.floor(n4_esatto))
        alpha4 = n4_esatto - N4
        
        self.delays.append(N4)
        self.weights.append(A4 * (1 - alpha4))
        self.delays.append(N4 + 1)
        self.weights.append(A4 * alpha4)
        
        self.buffer_size = max(self.delays) + 1
        self.history = deque([(0.0, 0.0, 0.0)] * self.buffer_size, maxlen=self.buffer_size)
        self.initialized = False

    def shape(self, q, dq, ddq):
        if not self.initialized:
            self.history = deque([(q, dq, ddq)] * self.buffer_size, maxlen=self.buffer_size)
            self.initialized = True
            
        self.history.appendleft((q, dq, ddq))
        
        shaped_q = sum(w * self.history[d][0] for w, d in zip(self.weights, self.delays))
        shaped_dq = sum(w * self.history[d][1] for w, d in zip(self.weights, self.delays))
        shaped_ddq = sum(w * self.history[d][2] for w, d in zip(self.weights, self.delays))
        
        return shaped_q, shaped_dq, shaped_ddq

# ==============================================================================
# CLASSE INPUT SHAPER (EI - Extra Insensitive due-gobbe)
# ==============================================================================
class EIShaper:
    def __init__(self, omega_n, zeta, Tc, V_max=0.05):
        # V_max = tolleranza vibrazione residua (default 5%)
        omega_d = omega_n * np.sqrt(1 - zeta**2) if zeta < 1.0 else omega_n
        K = np.exp(-(zeta * np.pi) / np.sqrt(1 - zeta**2)) if zeta < 1.0 else 1.0
        
        # Calcolo ampiezze ideali degli impulsi EI (3 impulsi come ZVD)
        Denom = (1 + K)**2
        A1 = (1 + V_max) / Denom
        A2 = (2 * K * (1 - V_max)) / Denom
        A3 = (K**2 * (1 + V_max)) / Denom
        
        self.delays = []
        self.weights = []
        
        # 1. Primo impulso
        self.delays.append(0)
        self.weights.append(A1)
        
        # 2. Secondo impulso (t = pi / omega_d)
        t2 = np.pi / omega_d
        n2_esatto = t2 / Tc 
        N2 = int(np.floor(n2_esatto))
        alpha2 = n2_esatto - N2 
        
        self.delays.append(N2)
        self.weights.append(A2 * (1 - alpha2)) 
        self.delays.append(N2 + 1)
        self.weights.append(A2 * alpha2)
        
        # 3. Terzo impulso (t = 2*pi / omega_d)
        t3 = 2 * np.pi / omega_d
        n3_esatto = t3 / Tc
        N3 = int(np.floor(n3_esatto))
        alpha3 = n3_esatto - N3
        
        self.delays.append(N3)
        self.weights.append(A3 * (1 - alpha3))
        self.delays.append(N3 + 1)
        self.weights.append(A3 * alpha3)
        
        self.buffer_size = max(self.delays) + 1
        self.history = deque([(0.0, 0.0, 0.0)] * self.buffer_size, maxlen=self.buffer_size)
        self.initialized = False

    def shape(self, q, dq, ddq):
        if not self.initialized:
            self.history = deque([(q, dq, ddq)] * self.buffer_size, maxlen=self.buffer_size)
            self.initialized = True
            
        self.history.appendleft((q, dq, ddq))
        
        shaped_q = sum(w * self.history[d][0] for w, d in zip(self.weights, self.delays))
        shaped_dq = sum(w * self.history[d][1] for w, d in zip(self.weights, self.delays))
        shaped_ddq = sum(w * self.history[d][2] for w, d in zip(self.weights, self.delays))
        
        return shaped_q, shaped_dq, shaped_ddq

# ==============================================================================
# CONFIGURAZIONE E AVVIO SIMULAZIONE
# ==============================================================================
with open(f'{model_name}/control_config.yaml', 'r') as file:
    params_yaml = yaml.safe_load(file)

with open(f'{model_name}/control_config.yaml', 'r') as file:
    params_yaml = yaml.safe_load(file)
    controller_params = params_yaml['controller']
    dynamic_params = np.array(params_yaml['model_parameters'])

xml_path = f"{model_name}/model.xml"
robot = MuJoCoMechanicalSystem(xml_path=xml_path, motor_actuators=["motor_1"], motor_joints=["joint_1"], spring_joints=[], ee_site="payload")
# robot = MuJoCoMechanicalSystem(xml_path=xml_path)
robot.show()
robot.initialize()
dof = robot.get_input_number()

Tc = robot.get_sampling_period()

# ----------------------------------------------------git p--------------------------
# INIZIALIZZAZIONE DELLO SHAPER CARICANDO I DATI SALVATI
# ------------------------------------------------------------------------------
shaper = None
if SHAPER_TYPE in ["ZV", "ZVD", "ZVDD", "EI"]:
    params_file = f"{model_name}/identified_params.yaml"

    if os.path.exists(params_file):
        with open(params_file, 'r') as f:
            saved_params = yaml.safe_load(f)
        
        omega_n = saved_params["omega_n"]
        zeta = saved_params["zeta"]
        
        if SHAPER_TYPE == "ZV":
            shaper = ZVShaper(omega_n, zeta, Tc)
        elif SHAPER_TYPE == "ZVD":
            shaper = ZVDShaper(omega_n, zeta, Tc)
        elif SHAPER_TYPE == "ZVDD":
            shaper = ZVDDShaper(omega_n, zeta, Tc)
        elif SHAPER_TYPE == "EI":
            shaper = EIShaper(omega_n, zeta, Tc) 
            
        print(f"[STATUS] Input Shaper {SHAPER_TYPE} ATTIVO (Dati da file): omega_n={omega_n:.3f} rad/s, zeta={zeta:.5f}")
    else:
        print(f"[ERRORE] File parametri non trovato in {params_file}. Esegui prima lo script con SHAPER_TYPE = 'NONE'!")
        SHAPER_TYPE = "NONE"
else:
    print("[STATUS] Esecuzione in modalità BASELINE (Shaper Disattivato)")
# -----------------------------------------------------------------------------


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
    if SHAPER_TYPE != "NONE" and shaper is not None:
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

# Suffix con nome in MAIUSCOLO
mode_suffix = f"shaped_{SHAPER_TYPE.upper()}" if SHAPER_TYPE != "NONE" else "baseline"

if SHAPER_TYPE == "NONE":
    filename = f"{model_name}/tests/{program_name}_{mode_suffix}.mat"
else:
    filename = f"{model_name}/tests/{program_name}_{mode_suffix}_{timestamp}.mat"

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
# IDENTIFICAZIONE E SALVATAGGIO PARAMETRI (Solo se SHAPER_TYPE == "NONE")
# ==============================================================================
if SHAPER_TYPE == "IDENTIFY":
    print("\n--- AVVIO FASE DI IDENTIFICAZIONE ---")
    T_oscillazione, zeta = extract_tz_from_data(filename, Tc)
    omega_n = (2 * np.pi) / T_oscillazione
    
    # Creiamo un dizionario con i dati identificati
    identified_data = {
        "T": float(T_oscillazione),
        "omega_n": float(omega_n),
        "zeta": float(zeta)
    }
    
    # Salviamo in un file YAML
    params_file = f"{model_name}/identified_params.yaml"
    with open(params_file, 'w') as f:
        yaml.dump(identified_data, f)
        
    print(f"[IDENTIFICAZIONE] Parametri fisici salvati con successo in: {params_file}\n")


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
fig1.update_layout(title=f"Tracking ({mode_suffix}): position / velocity / control", height=900, width=1200, legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0))
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
fig2.update_layout(title=f"Errors ({mode_suffix}): position error / velocity error / control", height=900, width=1200, legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0))
fig2.update_xaxes(showgrid=True)
fig2.update_yaxes(showgrid=True)

fig1.show()
fig2.show()
robot.close()