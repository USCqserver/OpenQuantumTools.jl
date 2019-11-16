using QuantumAnnealingTools, Test

control = QuantumAnnealingTools.InstPulseControl([0.5], (x)->σx)
@test control(1) == σx
