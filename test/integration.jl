using QTool, Test

@test scipy_quad((x)->x,0,1)[1] == 0.5
res = integrate_1d((x)->exp(-2*x),0,1)
@test isapprox(res[1], (1-exp(-2))/2, atol=1e-8, rtol=1e-6)
res = integrate_1d((x)->exp(-2*x),0,Inf)
@test isapprox(res[1], 0.5, atol=1e-8, rtol=1e-6)
res = integrate_1d((x)->exp(-x^2),-Inf,Inf)
@test isapprox(res[1], sqrt(pi), atol=1e-8, rtol=1e-6)
