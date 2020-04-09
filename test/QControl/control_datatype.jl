using QuantumAnnealingTools, Test
using DiffEqBase

mutable struct TestDEArray{T,N} <: DEDataArray{T,N}
    x::Array{T,N}
    test::Int
end

de_wrapper(x) = TestDEArray(x, 1)
# create FluctuatorControl
bath = EnsembleFluctuator([1.0, 2.0], [2.0, 1.0])
f_control = QuantumAnnealingTools.FluctuatorControl(10, 1, bath)
f_de_control = QuantumAnnealingTools.FluctuatorControl(10, 1, bath, :noise)
# create InstPulseControl
p_control = InstPulseControl([0.5], (x)->σx)
p_de_control = InstDEPulseControl([0.5], (x)->σx, :state)
# check_de_data_error test
u0 = PauliVec[1][1]
@test QuantumAnnealingTools.check_de_data_error(u0, p_control, nothing) == nothing
@test_throws ErrorException QuantumAnnealingTools.check_de_data_error(u0, p_de_control, nothing)
@test_throws ErrorException QuantumAnnealingTools.check_de_data_error(u0, p_de_control, de_wrapper)

# check ControlSet
@test_throws ErrorException control_set = ControlSet(p_control, p_control)
control_set = ControlSet(p_control, f_control)
@test control_set.pulse_control == p_control
@test control_set.fluctuator_control == f_control
