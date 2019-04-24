using QTool, Test

A = (s)->(1-s)
B = (s)->s
dθ= (s)->π/2
gap = (s)-> (cos(2*π*s) + 1)/2
h_test = AdiabaticFrameHamiltonian([dθ], [gap], [-real(σx)], [-real(σz)])
@test h_test(0.5) ≈ -π * real(σx) / 2
set_tf!(h_test, 10)
@test h_test(0.0) ≈ -π * real(σx) / 2 - 10 * real(σz)
pausing_test = construct_pausing_hamiltonian(0.5, 1.0, h_test)
@test pausing_test.tf_ext == 20
@test pausing_test(0.6) ≈ zeros(2, 2)
@test pausing_test(1.6) ≈ -π*real(σx)/2-(cos(1.2*π)+1)*10*real(σz)
h_test = AdiabaticFrameHamiltonian([dθ], [gap], [-real(σx)], [-real(σz)], unitless=false)
set_tf!(h_test, 10)
@test h_test(5) ≈ -π * real(σx) / 20
pausing_test = construct_pausing_hamiltonian(0.5, 1.0, h_test)
@test pausing_test(5) ≈ -π * real(σx)/40
@test pausing_test(6) ≈ zeros(2, 2)
@test pausing_test(16) ≈ -π*real(σx)/40-(cos(1.2*π)+1)*real(σz)/2
