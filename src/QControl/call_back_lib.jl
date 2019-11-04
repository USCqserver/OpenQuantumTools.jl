function ensemble_fluctuators_affect!(integrator)
    integrator.p.control
    for c in full_cache(integrator)
        c.pause = 1 - c.pause
    end
    u_modified!(integrator, false)
end
