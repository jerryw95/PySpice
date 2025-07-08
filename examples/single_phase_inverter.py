import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from PySpice.Spice.Netlist import Circuit
from PySpice.Unit import *
matplotlib.use('TkAgg')

# Define Global Parameters
dc_voltage = 200 @ u_V  # DC bus voltage
output_frequency = 50 @ u_Hz # Desired AC output frequency
carrier_frequency = 5 @ u_kHz # PWM switching frequency (much higher than output freq)
modulation_index = 0.8 # Mi = V_sine_peak / V_triangle_peak (should be <= 1 for linear modulation)

# Derived parameters
carrier_period = 1 / carrier_frequency
output_period = 1 / output_frequency

circuit = Circuit('Single Phase SPWM Inverter')

# 2. DC Voltage Source
circuit.V('input', 'dc_plus', 'dc_minus', dc_voltage)

# 3. Full-Bridge Inverter (using Voltage Controlled Switches)
# Switches are named S1, S2, S3, S4 (top-left, bottom-left, top-right, bottom-right)
# S1 and S4 are complementary, S2 and S3 are complementary
# S1 and S2 form one leg, S3 and S4 form the other leg

# Define a generic switch model
switch_model_name = 'IdealSwitch'
circuit.model(switch_model_name, 'SW', Ron=1@u_mOhm, Roff=1@u_GOhm, Vt=2.5@u_V, Vh=0.1@u_V) # Vt and Vh need to be chosen carefully based on control signal

# Leg 1 (Left Leg)
# S1: from dc_plus to bridge_out_a
# S2: from bridge_out_a to gnd
circuit.S('1', 'dc_plus', 'bridge_out_a', 'gate1_control', circuit.gnd, model=switch_model_name)
circuit.S('2', 'bridge_out_a', 'dc_minus', 'gate2_control', circuit.gnd, model=switch_model_name) # Assuming gate2 is complement of gate1

# Leg 2 (Right Leg)
# S3: from dc_plus to bridge_out_b
# S4: from bridge_out_b to gnd
circuit.S('3', 'dc_plus', circuit.gnd, 'gate2_control', circuit.gnd, model=switch_model_name)
circuit.S('4', circuit.gnd, 'dc_minus', 'gate1_control', circuit.gnd, model=switch_model_name) # Assuming gate4 is complement of gate3



# 3.1. Sinusoidal Reference Signal (Modulating Signal)
# V(t) = V_amplitude * sin(2*pi*f_output*t)
# We need to ensure its amplitude is compatible with the triangle wave for comparison.
# Let's say the triangle wave swings from -1V to +1V.
# So, the sine wave peak will be modulation_index * 1V.
sine_amplitude = modulation_index # Peak value of the sine wave for comparison

# circuit.V('modulating_sine', 'sine_ref', circuit.gnd, 'SIN({}, {}, {}Hz)'.format(0, sine_amplitude, output_frequency)) # DC offset, Amplitude, Frequency
circuit.SinusoidalVoltageSource('modulating_sine', 'sine_ref', circuit.gnd, offset=0, amplitude=sine_amplitude, frequency=output_frequency)

# 3.2. Triangular Carrier Signal
# V(t) = triangle_wave(freq, peak_amplitude)
# We need to define this as a PWL (Piece Wise Linear) source for PySpice
# A single period of a triangle wave from -Vpeak to +Vpeak
triangle_peak_voltage = 1 # Let's normalize carrier to -1 to +1

# Create a triangle wave using PWL. For a 5kHz carrier, period = 0.2ms
# Points: (0, -1V), (0.5*period, +1V), (period, -1V)
carrier_time_points = [0, carrier_period/2, carrier_period]
carrier_volt_points = [triangle_peak_voltage, -triangle_peak_voltage, triangle_peak_voltage]
# carrier_pwl_string = 'PWL( {values} ) REPEAT_LAST_POINT'.format(
#     values=' '.join(f'{t:.9f}V {v:.9f}' for t, v in zip(carrier_time_points, carrier_volt_points))
# )
circuit.PieceWiseLinearVoltageSource('carrier', 'triangle_carrier', circuit.gnd, repeat_time=0,
                                     values=[(t, v) for t, v in zip(carrier_time_points, carrier_volt_points)])

# 3.3. Comparators using Behavioral Voltage Sources (B-sources)
# Bipolar SPWM:
# Gate1 (S1) and Gate4 (S4) are ON when Sine_ref > Triangle_carrier
# Gate2 (S2) and Gate3 (S3) are ON when Sine_ref < Triangle_carrier
# We'll need a threshold for "on" and "off" of the gate signal for the switches.
# Let's say 5V for ON, 0V for OFF.

# Behavioral Source for Gate1 (S1) and Gate4 (S4)
# If V(sine_ref) > V(triangle_carrier), output 5V (ON), else 0V (OFF)
# The expression is `V = 5 * (V(sine_ref) > V(triangle_carrier))`
# PySpice allows conditional expressions in B-sources: `V = 'V(N1) > V(N2) ? 5 : 0'`
circuit.B('gate1_signal', 'gate1_control', circuit.gnd, v='(V(sine_ref) > V(triangle_carrier)) ? 5 : 0') # Gate for S1 and S4

# Behavioral Source for Gate2 (S2) and Gate3 (S3)
# If V(sine_ref) < V(triangle_carrier), output 5V (ON), else 0V (OFF)
circuit.B('gate2_signal', 'gate2_control', circuit.gnd, v='(V(sine_ref) < V(triangle_carrier)) ? 5 : 0') # Gate for S2 and S3

# The above assumes ideal switches with instant turn-on/off at 2.5V from the gate.
# Our switch model has Vt=2.5V, Vh=0.1V. So, 0V and 5V are good.

# Note on Dead Time:
# Implementing dead time directly with B-sources in PySpice can be complex
# as it often requires delay functions not directly available in standard SPICE B-sources.
# For a full-fledged simulation, you might need to:
# a) Generate 4 separate gate signals.
# b) Use multiple B-sources and logic to introduce delays.
# c) Consider using more advanced models for switches (e.g., MOSFETs) that might
#    incorporate some delay or have external gate driver circuits.
# For this example, we'll keep it ideal without explicit dead time,
# but be aware of its importance in real inverters.

# 4. LC Output Filter
# Choose L and C values for desired corner frequency and ripple reduction.
# Output frequency is 50Hz. Let's aim for a filter corner frequency much higher than 50Hz,
# but much lower than the switching frequency.
# Fc = 1 / (2 * pi * sqrt(LC))
# Let Fc = 500 Hz (1/10th of switching frequency)
# If L = 10mH
# C = 1 / ( (2*pi*Fc)^2 * L )
filter_frequency = 500 @ u_Hz
filter_L = 10 @ u_mH
# filter_C = 1 / ((2 * np.pi * filter_frequency.value)**2 * filter_L.value) @ u_F
# Round to a more standard capacitor value for simulation clarity
filter_C = 10 @ u_uF # Adjusted for a reasonable value

circuit.L('filter', 'bridge_out_a', 'output_mid', filter_L, initial_condition=0@u_A)
circuit.C('filter', 'output_mid', circuit.gnd, filter_C, initial_condition=0@u_V)

# 5. Load
load_resistance = 10 @ u_Ohm # Resistive load
circuit.R('load', 'output_mid', circuit.gnd, load_resistance)

# print(circuit)

# 6. Simulation
simulator = circuit.simulator(temperature=25, nominal_temperature=25)
# simulator.initial_condition(output_mid=0@u_V)

# Simulate for a few cycles of the output frequency
end_time = 3 * output_period # 3 cycles
step_time = carrier_period / 100 # 100 points per carrier cycle for good resolution

print(f"Simulating for {end_time*1e3:.2f} ms with step {step_time*1e6:.2f} us")
analysis = simulator.transient(step_time=step_time, end_time=end_time, use_initial_condition=True)

print(analysis)

# 7. Plotting Results
fig, ax = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

# Plot PWM control signals and carrier/reference
ax[0].plot(analysis.time * 1e3, analysis.sine_ref, label='Sine Reference (V)')
ax[0].plot(analysis.time * 1e3, analysis.triangle_carrier, label='Triangle Carrier (V)', linestyle='--')
ax[0].plot(analysis.time * 1e3, analysis.gate1_control, label='Gate1 (S1/S4) PWM (V)')
ax[0].plot(analysis.time * 1e3, analysis.gate2_control, label='Gate2 (S2/S3) PWM (V)', linestyle=':')
ax[0].set_ylabel('Voltage (V)')
ax[0].set_title('Sinusoidal PWM Generation')
ax[0].grid(True)
ax[0].legend()

# Plot Inverter Bridge Output (before filter)
ax[1].plot(analysis.time * 1e3, analysis.bridge_out_a, label='Leg Output (V_ao)', alpha=0.7)
# ax[1].plot(analysis.time * 1e3, analysis.bridge_out_b, label='Right Leg Output (V_bo)', alpha=0.7)
ax[1].set_ylabel('Voltage (V)')
ax[1].set_title('Inverter Leg Outputs (Before Filter)')
ax[1].grid(True)
ax[1].legend()

# Plot Filtered Output Voltage across Load
ax[2].plot(analysis.time * 1e3, analysis.output_mid, label='Filtered Output Voltage (V_o)')
ax[2].set_xlabel('Time (ms)')
ax[2].set_ylabel('Voltage (V)')
ax[2].set_title('Filtered AC Output Voltage Across Load')
ax[2].grid(True)
ax[2].legend()

plt.tight_layout()
plt.show()

# # Optional: Perform FFT to check harmonic content (requires scipy)
# try:
#     from scipy.fft import fft, fftfreq
#     load_voltage = np.array(analysis.output_mid)
#     time_data = np.array(analysis.time)
#     sample_rate = 1 / (time_data[1] - time_data[0]) # Approximate sample rate
#
#     N = len(load_voltage)
#     yf = fft(load_voltage)
#     xf = fftfreq(N, 1/sample_rate)
#
#     # Take only positive frequencies and magnitude
#     xf_pos = xf[0:N//2]
#     yf_pos = 2.0/N * np.abs(yf[0:N//2])
#
#     plt.figure(figsize=(10, 6))
#     plt.plot(xf_pos, yf_pos)
#     plt.title('FFT of Filtered Output Voltage')
#     plt.xlabel('Frequency (Hz)')
#     plt.ylabel('Magnitude (V)')
#     plt.xlim(0, carrier_frequency.value * 2) # Show up to twice the carrier frequency
#     plt.grid(True)
#     plt.show()
#
#     # Find fundamental component
#     fundamental_index = np.argmin(np.abs(xf_pos - output_frequency.value))
#     fundamental_voltage = yf_pos[fundamental_index]
#     print(f"Fundamental output voltage (peak): {fundamental_voltage:.2f} V")
#
# except ImportError:
#     print("Scipy not installed. Skipping FFT analysis. Run 'pip install scipy' to enable.")