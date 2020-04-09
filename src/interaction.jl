struct Interaction
    coupling::AbstractCouplings
    bath::AbstractBath
end

function create_redfield(inter::Interaction, unitary, tf)
    create_redfield(inter.coupling, unitary, tf, inter.bath)
end

struct InteractionSet{T<:Tuple}
    interactions::T
end

function InteractionSet(inters::Interaction...)
    InteractionSet(inters)
end

function create_redfield(inter::InteractionSet, unitary, tf)
    RedfieldSet([
        create_redfield(i.coupling, unitary, tf, i.bath)
        for i in inter.interactions
    ]...)
end
