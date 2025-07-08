# ####################################################################################################
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
####################################################################################################

import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()

####################################################################################################

from PySpice.Doc.ExampleTools import find_libraries
from PySpice.Probe.Plot import plot
from PySpice.Spice.Library import SpiceLibrary
from PySpice.Spice.Netlist import Circuit, SubCircuit
from PySpice.Unit import *

####################################################################################################

libraries_path = find_libraries()
spice_library = SpiceLibrary(libraries_path)

####################################################################################################
class Switch_snubber(SubCircuit):
    NODES = ('n+', 'n-', 'nc+', 'nc-')
    def __init__(self, name, R=1 @u_uOhm, C = 0@u_uF):
        SubCircuit.__init__(self, name, *self.NODES)
        if C.value > 0:
            self.R(1, 'n+', 'n', R)
            self.C(1, 'n', 'n-', C)
        else:
            self.R(1, 'n+', 'n-', R)
        self.S(1, 'n+', 'n-', 'nc+', 'nc-', model='IdealSwitch')

class Diode_snubber(SubCircuit):
    NODES = ('n+', 'n-')
    def __init__(self, name, R=1 @u_uOhm, C = 0@u_uF):
        SubCircuit.__init__(self, name, *self.NODES)
        if C.value > 0:
            self.R(1, 'n+', 'n', R)
            self.C(1, 'n', 'n-', C)
        else:
            self.R(1, 'n+', 'n-', R)
        self.Diode(1, 'n+', 'n-', model='MyDiode')

#?# circuit_macros('buck-converter.m4')

circuit = Circuit('Boost Converter')

# circuit.include(spice_library['1N5822']) # Schottky diode
# circuit.include(spice_library['irf150'])

# From Microchip WebSeminars - Buck Converter Design Example

Vin = 250@u_V
Vout = 500@u_V
ratio = Vin / Vout


frequency = 10@u_kHz
period = frequency.period
duty_cycle = (1 - ratio) * period

print('ratio =', ratio)
print('period =', period.canonise())
print('duty_cycle =', duty_cycle.canonise())


circuit.model('IdealSwitch', 'SW', Ron=1@u_mOhm, Roff=1@u_GOhm, Vt=2.5@u_V, Vh=0.1@u_V) # Vt and Vh need to be chosen carefully based on control signal

circuit.model('MyDiode', 'D', Is=1@u_nA, N=1.984, Rs=0.1@u_Ohm)

circuit.subcircuit(Switch_snubber('Switch_snubber', R=1e5@u_Ohm))
circuit.subcircuit(Diode_snubber('Diode_snubber', R=1e5@u_Ohm))



circuit.V('in', 'in', circuit.gnd, Vin)
circuit.L('in', 'in', 1, 300@u_uH)

# circuit.S('switch', 1, circuit.gnd, 'switch_control', circuit.gnd, model='IdealSwitch')
circuit.X('switch', 'Switch_snubber', 1, circuit.gnd, 'switch_control', circuit.gnd)
circuit.PulseVoltageSource('pulse', 'switch_control', circuit.gnd, 0@u_V, 5@u_V, duty_cycle, period)

circuit.Diode('Diode', 1, 'out', model='MyDiode')
circuit.R('out', 'out', circuit.gnd, 10@u_Ohm)
circuit.C('out', 'out', circuit.gnd, 200@u_uF, initial_condition=50@u_V) # , initial_condition=0@u_V

print(circuit)

simulator = circuit.simulator(temperature=25, nominal_temperature=25)
analysis = simulator.transient(step_time=period/1000, end_time=period*100, use_initial_condition=False)

fig, ax = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

ax[0].plot(analysis.time, analysis.out, label='output (V)')
ax[0].plot(analysis.time, analysis['in'], label='input (V)')
ax[0].grid()
ax[0].set_xlabel('t [s]')
ax[0].set_ylabel('[V]')
ax[0].legend()

ax[1].plot(analysis.time, analysis['1'], label='Switch (V)')
ax[1].grid()
ax[1].set_xlabel('t [s]')
ax[1].set_ylabel('[V]')
ax[1].legend()

ax[2].plot(analysis.time, analysis.switch_control, label='Switch control (V)')
ax[2].grid()
ax[2].set_xlabel('t [s]')
ax[2].set_ylabel('[V]')
ax[2].legend()

plt.tight_layout()
plt.show()
#
# #f# save_figure('figure', 'buck-converter.png')




# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
# matplotlib.use('TkAgg')
# from PySpice.Spice.Netlist import Circuit
# from PySpice.Unit import *
# from PySpice.Doc.ExampleTools import find_libraries
# from PySpice.Spice.Library import SpiceLibrary
#
# # --- 1. Define Circuit Parameters ---
# input_voltage = 12 @ u_V      # Input DC voltage
# output_voltage_desired = 24 @ u_V # Desired output voltage (ideal)
#
# switching_frequency = 10 @ u_kHz # PWM switching frequency
# switching_period = 1 / switching_frequency
#
# # Calculate ideal duty cycle for continuous conduction mode (CCM)
# # Vout = Vin / (1 - D)  => 1 - D = Vin / Vout => D = 1 - (Vin / Vout)
# duty_cycle = 1 - (input_voltage.value / output_voltage_desired.value)
# print(f"Calculated Ideal Duty Cycle: {duty_cycle:.3f}")
#
# # Pulse width for the switch ON time
# pulse_width = duty_cycle * switching_period
#
# # Component values (design based on typical boost converter design equations)
# # These values are chosen for reasonable ripple and stability at 50kHz.
# inductor_value = 300 @ u_uH  # Inductor
# capacitor_value = 100 @ u_uF # Output capacitor
# load_resistance = 10 @ u_Ohm # Load resistor
#
# libraries_path = "D:\\Work_Learn\\projects\\Datacenter_power\\PySpice\\examples\\libraries"
# spice_library = SpiceLibrary(libraries_path)
# # --- 2. Create the Circuit ---
# circuit = Circuit('Basic Boost Converter')
#
# # --- 3. Add Components ---
#
# # 3.1. DC Input Voltage Source
# circuit.V('in', 'vin', circuit.gnd, input_voltage)
#
# # 3.2. Inductor
# # Connected between input voltage source and the switching node
# # circuit.L('L1', 'vin', 'switch_node', inductor_value)
# # Set initial current through inductor to 0A for a clean start
# circuit.L('L1', 'vin', 'drain_node', inductor_value, initial_condition=0@u_A)
#
#
# # 3.3. Ideal Switch (MOSFET equivalent, controlled by PWM)
# # The switch connects 'switch_node' to 'gnd'.
# # When ON, it shorts the inductor to ground, allowing current to build.
# switch_model_name = 'MyBoostSwitch'
# circuit.model(switch_model_name, 'SW',
#               Ron=1@u_mOhm,      # Very low on-resistance
#               Roff=1@u_GOhm,  # Very high off-resistance
#               Vt=2.5@u_V,        # Threshold voltage for control (e.g., 2.5V for a 0-5V gate signal)
#               Vh=0.1@u_V)       # Hysteresis for smooth switching
# circuit.include(spice_library['1N5822']) # Schottky diode
# circuit.include(spice_library['irf150'])
#
# circuit.S('S', 'drain_node', circuit.gnd, 'pwm_control', circuit.gnd, model='MyBoostSwitch')
# # circuit.X('Q', 'irf150', 'drain_node', 'base_node',  circuit.gnd)
# circuit.R('R_base', 'pwm_control', 'base_node', 1@u_Ohm) # PWM output to base through resistor
#
#
# # 3.4. Diode
# # The diode connects 'switch_node' to the output capacitor, allowing current flow
# # only when the inductor voltage goes higher than the output voltage.
# # An ideal diode model often works well for basic simulations.
# diode_model_name = 'IdealDiode'
# # The D device takes (positive_node, negative_node, model_name)
# circuit.model(diode_model_name, 'D',
#               Is=1@u_pA,    # Saturation current (small for ideal)
#               Rs=0.001@u_Ohm, # Series resistance (small for ideal)
#               N=1)          # Emission coefficient (close to 1 for ideal)
#               # If you need reverse breakdown, use BV and IBV (e.g., BV=100V, IBV=10uA)
#
# circuit.X('D', '1N5822','drain_node', 'output_node')
# # circuit.Diode('D', 'drain_node', 'output_node', model = 'IdealDiode')
#
# # 3.5. Output Capacitor
# # Filters the pulsed voltage after the diode
# # circuit.C('C1', 'output_node', circuit.gnd, capacitor_value)
# # Set initial voltage across capacitor to 0V for a clean start
# circuit.C('C1', 'output_node', circuit.gnd, capacitor_value, initial_condition=0@u_V)
#
# # 3.6. Load Resistor
# circuit.R('Rload', 'output_node', circuit.gnd, load_resistance)
#
# # --- 4. PWM Signal Generation ---
# # This pulsed voltage source will drive the gate of the switch
# # It generates a square wave to turn the switch ON and OFF
# circuit.PulseVoltageSource('pwm_gen', 'pwm_control', circuit.gnd,
#           initial_value=0@u_V,        # Low state (switch OFF)
#           pulsed_value=5@u_V,         # High state (switch ON)
#           rise_time= 10@u_ns,
#           fall_time= 10@u_ns,
#           pulse_width=pulse_width,      # ON time
#           period=switching_period)      # Total period
#
# # --- 5. Simulation ---
# simulator = circuit.simulator(temperature=25, nominal_temperature=25)
#
# # Simulate for several switching cycles to see it reach steady-state
# # A few hundred cycles are usually sufficient for boost converters to settle
# num_cycles_to_simulate = 100
# end_time = num_cycles_to_simulate * switching_period
# step_time = switching_period / 1000 # 1000 points per switching period for good resolution
#
# print(f"Simulating for {end_time*1e6:.2f} us ({num_cycles_to_simulate} cycles) with step {step_time*1e9:.2f} ns")
#
# # Crucial: Use initial_state=True to respect the specified initial conditions (L1=0A, C1=0V)
# analysis = simulator.transient(step_time=step_time, end_time=end_time, use_initial_condition=False)
#
# # --- 6. Plotting Results ---
# fig, ax = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
#
# # Plot Inductor Current
# # # PySpice allows you to access inductor current using I(L_name)
# ax[0].plot(analysis.time * 1e6, analysis.drain_node, label='Switch node (V)')
# ax[0].set_ylabel('Voltage (V)')
# ax[0].set_title('Boost Converter Simulation')
# ax[0].grid(True)
# ax[0].legend()
#
# # Plot Switch Control Voltage (PWM)
# ax[1].plot(analysis.time * 1e6, analysis.pwm_control, label='PWM Control (V)')
# ax[1].set_ylabel('Voltage (V)')
# ax[1].grid(True)
# ax[1].legend()
#
# # Plot Output Voltage
# ax[2].plot(analysis.time * 1e6, analysis.output_node, label='Output Voltage (V_out)')
# ax[2].axhline(y=output_voltage_desired.value, color='r', linestyle='--', label='Desired Output Voltage')
# ax[2].set_xlabel('Time (us)')
# ax[2].set_ylabel('Voltage (V)')
# ax[2].grid(True)
# ax[2].legend()
#
# plt.tight_layout()
# plt.show()
