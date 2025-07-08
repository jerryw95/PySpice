import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()

from PySpice.Doc.ExampleTools import find_libraries
from PySpice.Probe.Plot import plot
from PySpice.Spice.Library import SpiceLibrary
from PySpice.Spice.Netlist import Circuit, SubCircuit
from PySpice.Unit import *

figure, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(18, 9))

libraries_path = find_libraries()
spice_library = SpiceLibrary(libraries_path)

class Switch_snubber(SubCircuit):
    NODES = ('n+', 'n-', 'nc+', 'nc-')
    def __init__(self, name, R=1 @u_uOhm, C = 0@u_uF):
        SubCircuit.__init__(self, name, *self.NODES)
        if C.value > 0:
            self.R(1, 'n+', 'n', R)
            self.C(1, 'n', 'n-', C)
        else:
            self.R(1, 'n+', 'n-', R)
        self.S(1, 'n+', 'n-', 'nc+', 'nc-', model='MySwitch')

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

circuit = Circuit('uninterruptible power supply')

circuit.model('MySwitch', 'SW', Ron=1@u_mOhm, Roff=1@u_GOhm, Vt=2.5@u_V, Vh=0.1@u_V)
circuit.model('MyDiode', 'D', Is=1@u_nA, N=1.984, Rs=0.1@u_Ohm)

circuit.subcircuit(Switch_snubber('Switch_snubber', R=1e8@u_Ohm))
circuit.subcircuit(Diode_snubber('Diode_snubber', R=1e8@u_Ohm))


################################################################################################
# Adding the AC source and rectifier components
source = circuit.SinusoidalVoltageSource('input', 'src_pos', 'src_neg', amplitude=230, frequency=50)
circuit.X('D1', 'Diode_snubber', 'src_pos', 'rec_out')
circuit.X('D2', 'Diode_snubber', circuit.gnd, 'src_neg')
circuit.X('D3', 'Diode_snubber', 'src_neg', 'rec_out')
circuit.X('D4', 'Diode_snubber', circuit.gnd, 'src_pos')
circuit.C('1', 'rec_out', circuit.gnd, 10@u_mF,initial_condition=230@u_V)
# circuit.R('load', 'rec_out', circuit.gnd, 100@u_Ω)

################################################################################################
# Adding the Boost converter component
frequency = 10@u_kHz
period = frequency.period
ratio = 0.6
duty_cycle = (1 - ratio) * period
circuit.X('D5', 'Diode_snubber', 'rec_out', 'boost_in')
circuit.L('L_boost', 'boost_in', 'boost_switch', 200@u_uH)
circuit.X('switch', 'Switch_snubber','boost_switch', circuit.gnd, 'switch_control', circuit.gnd)
circuit.X('D6', 'Diode_snubber', 'boost_switch', 'boost_out')
circuit.PulseVoltageSource('pulse', 'switch_control', circuit.gnd, 0@u_V, 5@u_V, duty_cycle, period)
circuit.R('boost', 'boost_out', circuit.gnd, 10@u_Ohm)
circuit.C('boost', 'boost_out', circuit.gnd, 200@u_uF, initial_condition=300@u_V) # , initial_condition=0@u_V


################################################################################################
# Adding the DC-AC inverter
dc_voltage = 300 @ u_V  # DC bus voltage
output_frequency = 50 @ u_Hz # Desired AC output frequency
carrier_frequency = 10 @ u_kHz # PWM switching frequency (much higher than output freq)
modulation_index = 0.8 # Mi = V_sine_peak / V_triangle_peak (should be <= 1 for linear modulation)
carrier_period = 1 / carrier_frequency
output_period = 1 / output_frequency
circuit.X('S1', 'Switch_snubber', 'boost_out', 'bridge_out_a', 'gate1_control', circuit.gnd)
circuit.X('S2', 'Switch_snubber', 'bridge_out_b', circuit.gnd, 'gate1_control', circuit.gnd)
circuit.X('S3', 'Switch_snubber', 'boost_out', 'bridge_out_b', 'gate3_control', circuit.gnd)
circuit.X('S4', 'Switch_snubber', 'bridge_out_a', circuit.gnd, 'gate3_control', circuit.gnd)
# Generate modulation signals
sine_amplitude = modulation_index # Peak value of the sine wave for comparison
circuit.SinusoidalVoltageSource('modulating_sine', 'sine_ref', circuit.gnd, offset=0, amplitude=sine_amplitude, frequency=output_frequency)
triangle_peak_voltage = 1
carrier_time_points = [0, carrier_period/2, carrier_period]
carrier_volt_points = [triangle_peak_voltage, -triangle_peak_voltage, triangle_peak_voltage]
circuit.PieceWiseLinearVoltageSource('carrier', 'triangle_carrier', circuit.gnd, repeat_time=0,
                                     values=[(t, v) for t, v in zip(carrier_time_points, carrier_volt_points)])
circuit.B('gate1_signal', 'gate1_control', circuit.gnd, v='(V(sine_ref) > V(triangle_carrier)) ? 5 : 0') # Gate for S1 and S2
circuit.B('gate3_signal', 'gate3_control', circuit.gnd, v='(V(sine_ref) < V(triangle_carrier)) ? 5 : 0') # Gate for S3 and S4
# LC filtering and load
circuit.L('filter', 'bridge_out_a', 'output', 10 @ u_mH, initial_condition=0@u_A)
circuit.C('filter', 'output', 'bridge_out_b', 10 @ u_uF, initial_condition=0@u_V)
circuit.R('load', 'output', 'bridge_out_b', 10 @ u_Ohm)


################################################################################################
# Simulation settings and visualization
simulator = circuit.simulator(temperature=25, nominal_temperature=25)
analysis = simulator.transient(step_time=period/200, end_time=source.period*2, use_initial_condition=True)

ax1.set_title('Full-Wave Rectification with filtering')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Voltage [V]')
ax1.grid()
ax1.plot(analysis.time, analysis['src_pos'])
ax1.plot(analysis.time, analysis.rec_out)
ax1.legend(('input', 'output'), loc=(.05,.1))


ax2.set_title('DC-DC boost converter with filtering')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Voltage [V]')
ax2.grid()
ax2.plot(analysis.time, analysis.boost_switch)
ax2.plot(analysis.time, analysis.boost_out)
ax2.legend(('switch', 'output'), loc=(.05,.1))

ax3.set_title('DC-AC inverter with filtering')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Voltage [V]')
ax3.grid()
ax3.plot(analysis.time, analysis.bridge_out_a - analysis.bridge_out_b)
ax3.plot(analysis.time, analysis.output - analysis.bridge_out_b)
ax3.legend(('inverter_out', 'output'), loc=(.05,.1))

plt.tight_layout()

plt.show()