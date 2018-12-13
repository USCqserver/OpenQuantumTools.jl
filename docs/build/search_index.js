var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "HOME",
    "title": "HOME",
    "category": "page",
    "text": ""
},

{
    "location": "#QTool.Hamiltonian",
    "page": "HOME",
    "title": "QTool.Hamiltonian",
    "category": "type",
    "text": "Object to hold general time dependent Hamiltonians, whose form is assumed to be a summation of time dependent function times constant matrices. \n\n\n\n\n\n"
},

{
    "location": "#QTool.Hamiltonian-Tuple{Real}",
    "page": "HOME",
    "title": "QTool.Hamiltonian",
    "category": "method",
    "text": "Evaluate the Hamiltonian at time t \n\n\n\n\n\n"
},

{
    "location": "#QTool.comm!-Tuple{Any,Any,Any}",
    "page": "HOME",
    "title": "QTool.comm!",
    "category": "method",
    "text": "comm!(Y, A, B)\n\nCalculate the commutator of matrices A and B. Write the result in Y\n\n\n\n\n\n"
},

{
    "location": "#QTool.comm-Tuple{Any,Any}",
    "page": "HOME",
    "title": "QTool.comm",
    "category": "method",
    "text": "comm(A, B)\n\nCalculate the commutator of matrices A and B\n\n\n\n\n\n"
},

{
    "location": "#QTool.construct_hamming_weight_op-Tuple{Int64,String}",
    "page": "HOME",
    "title": "QTool.construct_hamming_weight_op",
    "category": "method",
    "text": "construct_hamming_weight_op(num_qubit::Int64, op::String)\n\nConstruct the Hamming weight operator for system of size num_qubit. The type of the Hamming weight operator is specified by op: \"x\", \"y\" or \"z\".\n\nExamples\n\njulia> construct_hamming_weight_op(2,\"z\")\n4×4 Array{Complex{Float64},2}:\n 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im\n 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im\n 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im\n 0.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im\n\n\n\n\n\n"
},

{
    "location": "#QTool.eigen_state_eval-Tuple{Any,AbstractArray{Float64,1}}",
    "page": "HOME",
    "title": "QTool.eigen_state_eval",
    "category": "method",
    "text": "eigen_state_eval(hfun, t[, levels=[1,]])\n\nCalculate the eigen states of Hamiltonian hfun at each points of t. The output levels are specified by levels argument.\n\n\n\n\n\n"
},

{
    "location": "#QTool.eigen_value_eval-Tuple{Any,AbstractArray{Float64,1}}",
    "page": "HOME",
    "title": "QTool.eigen_value_eval",
    "category": "method",
    "text": "eigen_value_eval(hfun, t[, levels])\n\nCalculate the eigen values of Hamiltonian hfun at each points of t. The output levels are specified by levels argument.\n\n\n\n\n\n"
},

{
    "location": "#QTool.inst_population-Tuple{Any,Any,Any}",
    "page": "HOME",
    "title": "QTool.inst_population",
    "category": "method",
    "text": "function inst_population(t, states, hamiltonian; level=1)\n\nFor a time series quantum states given by states, whose time points are given by t, calculate the population of instantaneous eigenstates of hamiltonian. The levels of the instantaneous eigenstates are specified by level, which can be any slice index.\n\n\n\n\n\n"
},

{
    "location": "#QTool.matrix_decompose-Union{Tuple{T}, Tuple{Array{T,2},Array{Array{T,2},1}}} where T<:Number",
    "page": "HOME",
    "title": "QTool.matrix_decompose",
    "category": "method",
    "text": "matrix_decompose(mat::Array{T,2}, basis::Array{Array{T,2},1})\n\nDecompse matrix mat onto matrix basis basis\n\nExamples\n\njulia> matrix_decompose(1.0*σx+2.0*σy+3.0*σz, [σx,σy,σz])\n3-element Array{Complex{Float64},1}:\n 1.0 + 0.0im\n 2.0 + 0.0im\n 3.0 + 0.0im\n\n\n\n\n\n"
},

{
    "location": "#QTool.q_translate-Tuple{String}",
    "page": "HOME",
    "title": "QTool.q_translate",
    "category": "method",
    "text": "q_translate(h::String)\n\nConvert a string h representing multi-qubits Pauli matrices summation into its numerical form.\n\nExamples\n\njulia> q_translate(\"X+2.0Z\")\n2×2 Array{Complex{Float64},2}:\n 2.0+0.0im   1.0+0.0im\n 1.0+0.0im  -2.0+0.0im\n\n\n\n\n\n"
},

{
    "location": "#QTool.temperature_2_beta-Tuple{Any}",
    "page": "HOME",
    "title": "QTool.temperature_2_beta",
    "category": "method",
    "text": "temperature_2_beta(T; unit=:h)\n\nConvert temperature from mK to unit*β\n\n\n\n\n\n"
},

{
    "location": "#QTool.temperature_2_freq-Tuple{Any}",
    "page": "HOME",
    "title": "QTool.temperature_2_freq",
    "category": "method",
    "text": "temperature_2_freq(T)\n\nConvert temperature from mK to GHz(physical unit)\n\n\n\n\n\n"
},

{
    "location": "#HOME-1",
    "page": "HOME",
    "title": "HOME",
    "category": "section",
    "text": "Modules = [QTool]\nPrivate = false"
},

]}
