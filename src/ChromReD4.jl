module ChromReD4
export Nucleus
export plot
export update!
export init_nucleus


import Plots; const plt = Plots; plt.gr();
using Random

struct Chromatin
    positions::Array{Int, 2}
end


mutable struct Nucleus
    chromatin_pos::Chromatin
    chromatin_state::Array{Int, 2}
    chromatin_dist::Array{Int, 2}
    size::Array{Int, 1}
    rad4_state::Array{Float16, 2}
end


function diffusion(nucleus::Nucleus, diffu_coeff::Float64=6.7)
    x_influence = circshift(nucleus.rad4_state, (0, 1)) .+ circshift(nucleus.rad4_state, (0, -1))
    y_influence = circshift(nucleus.rad4_state, (1, 0)) .+ circshift(nucleus.rad4_state, (-1, 0))
    return diffu_coeff .* (x_influence .+ y_influence .- 4nucleus.rad4_state)
end


function associate(nucleus::Nucleus, prob::Float64=.4)
    interact_idc = findall((nucleus.rad4_state .> 0) .& (nucleus.chromatin_state .> 0))
    asso_candidates = rand(length(interact_idc)) .< prob

    rad4_delta = zeros(nucleus.size...)
    chromatin_delta = zeros(nucleus.size...)
    if length(interact_idc[asso_candidates]) > 0
        rad4_delta[interact_idc[asso_candidates]] .= -1
        chromatin_delta[interact_idc[asso_candidates]] .= -1
    end
    return rad4_delta, chromatin_delta
end


function dissociate(nucleus::Nucleus, prob::Float64=.7)
    interact_idc = findall(nucleus.chromatin_dist .- nucleus.chromatin_state .> 0)  
    disso_candidates = rand(length(interact_idc)) .< prob

    rad4_delta = zeros(nucleus.size...)
    chromatin_delta = zeros(nucleus.size...)
    if length(interact_idc[disso_candidates]) > 0
        rad4_delta[interact_idc[disso_candidates]] .= 1
        chromatin_delta[interact_idc[disso_candidates]] .= 1
    end
    return rad4_delta, chromatin_delta
end


function init_chromatin(num_iter::Int=10, size_n::Array{Int, 1}=[100, 100])::Chromatin
    verify_boundary_conditions = pos -> min.(size_n, max.(1, pos))

    chromatin = Chromatin(Array{Int, 2}(undef, num_iter, 2))

    start_pos = [rand(1:size_n[1]), rand(1:size_n[2])]
    chromatin.positions[1, :] = start_pos
    for i = 2:num_iter
        chromatin.positions[i, :] = verify_boundary_conditions(chromatin.positions[i-1, :] .+ [rand(-1:1), rand(-1:1)])
    end
    
    return chromatin
end


function init_nucleus(len_chromatin::Int=10000, size_n::Array{Int, 1}=[100, 100])::Nucleus
    chromatin = init_chromatin(len_chromatin, size_n)
    rad4_state = rand(0:1, size_n...)
    chromatin_dist = zeros(size_n...)
    chromatin_state = zeros(size_n...)
    for i in 1:size(chromatin.positions, 1)
        x, y = chromatin.positions[i, :]
        chromatin_state[y, x] += 1
        chromatin_dist[y, x] += 1
    end
    nucleus = Nucleus(chromatin, chromatin_state, chromatin_dist, size_n, rad4_state)
    return nucleus
end


function update!(nucleus::Nucleus; asso_prob::Float64=.1, disso_prob::Float64=.7, diffu::Float64=7.2)
    asso_rad4, asso_chrom = associate(nucleus, asso_prob)
    disso_rad4, disso_chrom = dissociate(nucleus, disso_prob)
    diffu_update = diffusion(nucleus, diffu)
    nucleus.rad4_state += asso_rad4 + disso_rad4 + diffu_update
    nucleus.chromatin_state += asso_chrom + disso_chrom
    return nothing
end


function plot(nucleus::Nucleus, title_appendix::String="")
    graph = plt.plot(
        nucleus.chromatin_pos.positions[:, 1],
        nucleus.chromatin_pos.positions[:, 2]
    )

    plt.heatmap!(
        1:size(nucleus.rad4_state,1),
        1:size(nucleus.rad4_state,2),
        nucleus.rad4_state,
        c=plt.cgrad([:blue, :white, :red]), 
        interpolate=true, 
        alpha=.4,
        clim=(-5, 5)
    )

    plt.title!("Rad4 Distribution $title_appendix")
    plt.xlabel!("x Direction")
    plt.ylabel!("y Direction")
    plt.xlims!((0, nucleus.size[1]))
    plt.ylims!((0, nucleus.size[2]))

    return graph
end
end