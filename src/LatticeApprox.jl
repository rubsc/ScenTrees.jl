
using Statistics, LinearAlgebra, PyPlot

mutable struct Lattice
    name::String
    state::Array{Array{Float64,2},1}
    probability::Array{Array{Float64,2},1}
    Lattice(name::String, state::Array{Array{Float64,2},1}, probability::Array{Array{Float64,2},1}) = new(name, state, probability)
end

"""
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
	LatticeApproximation()

Returns a valuated scanario lattice. It takes as inputs a vector of branching structure and the sample size.

"""


function LatticeApproximation(states::Array{Int64,1},nScenarios::Int64)
    WassersteinDistance = 0.0
    rWasserstein = 2
    lns = length(states)
    LatState = [zeros(states[j],1) for j = 1:lns]                                             # States of the lattice at each time t
    LatProb = vcat([zeros(states[1],1)],[zeros(states[j-1],states[j]) for j = 2:lns])          # Probabilities of the lattice at each time t

    #Stochastic approximation
    for n = 1:nScenarios
        #Z = 300 .+ 50 * vcat(0.0,cumsum(randn(lns-1,1),dims = 1))
        Z = vcat(0.0,cumsum(randn(lns-1,1),dims = 1))                                           #draw a new sample Gaussian path
        idtm1 = 1
        dist = 0.0
        for t = 1:length(states)                                                                #walk along the gradient
            sumLat = sum(LatProb[t],dims = 2)
            sqStates = 1.3 * sqrt(n) / states[t]
            tmp = Int64[id for (id,ls) in enumerate(sumLat) if ls < sqStates]
            LatState[t][tmp] .= Z[t]                                                            #corrective action to include lost nodes
            mindist,indx = findmin(vec(abs.(LatState[t] .- Z[t])))                              #find the closest lattice entry
            dist = dist + mindist^2                                                             #Euclidean distance for the paths
            LatProb[t][idtm1,indx] = LatProb[t][idtm1,indx] .+ 1.0                              #increase the probability
            idtm1 = indx
            LatState[t][indx] = LatState[t][indx] - 2 / (3000 + n)^0.75 *rWasserstein * mindist^(rWasserstein-1) * (LatState[t][indx] - Z[t])
        end
        dist = dist^(1/2)
        WassersteinDistance = (WassersteinDistance*(n-1) + dist^rWasserstein)/n
    end
    LatProb = LatProb ./ nScenarios                                               #scale the probabilities to 1.0
    return Lattice("Lattice Approximation of $states, \n distance=$(round(WassersteinDistance^(1/rWasserstein),digits = 4)) at $(nScenarios) scenarios",LatState,LatProb)
end

=======
	LatticeApproximation(states::Array{Int64,1},path::Function,nScenarios::Int64)
=======
	LatticeApproximation(states::Array{Int64,1},path::Function,nIterations::Int64)
>>>>>>> master
=======
	LatticeApproximation(bstructure::Array{Int64,1}, path::Function, nIterations::Int64, r::Int64 = 2)
>>>>>>> master
=======
	lattice_approximation(bstructure::Array{Int64,1}, path::Function, nIterations::Int64, r::Int64 = 2)
>>>>>>> master

Returns an approximated lattice approximating the stochastic process provided.

Args:
bstructure - branching structure of the scenario lattice e.g., bstructure = [1,2,3,4,5] represents a 5-staged lattice
path - function generating samples from known distribution with length equal to the length of bstructure of the lattice.
nIterations - number of iterations for stochastic approximation algorithm.
r - parameter for the transportation distance.
"""
function lattice_approximation(bstructure::Array{Int64,1}, path::Function, nIterations::Int64, r::Int64 = 2)
    t_dist = 0.0                               # multistage distance
    T = length(bstructure)     # number of stages in the scenario lattice
    states = [zeros(bstructure[j], 1) for j = 1 : T]   # States of the lattice at each time t
    init_lat = path()         # initialize the states of the nodes of the lattice
    for t = 1:T
        states[t] .= init_lat[t]
    end
    probabilities = vcat([zeros(bstructure[1], 1)], [zeros(bstructure[j-1], bstructure[j]) for j = 2 : T]) # Probabilities of the lattice at each time t
    Z = Array{Float64}(undef,T,1)  # initialize a 1D vector to hold the state values
    #Stochastic approximation procedure to update the values on the nodes of the scenario lattice.
    for n = 1:nIterations
        Z .= path(); last_index = 1; dist = 0.0
        for t = 1 : T #walk along the gradient
            # corrective action to include lost nodes
            states[t][findall(sum(probabilities[t], dims = 2) .< 1.3 * sqrt(n) / bstructure[t])] .= Z[t]
            min_dist, new_index = findmin(abs.(vec(states[t] .- Z[t])))  # find the closest lattice entry
            dist += min_dist^2                             # Euclidean distance for the paths
            probabilities[t][last_index,new_index] += 1.0        # increase the probability
            last_index = new_index
            states[t][new_index] = states[t][new_index] - (r * min_dist^(r-1) * (states[t][new_index] - Z[t]) / (3000 + n)^0.75)
        end
        dist = dist^(1/2)
        t_dist = (t_dist*(n-1) + dist^r)/n
    end
    probabilities = probabilities / nIterations						                # scale the probabilities to 1.0
    return Lattice("Lattice $bstructure with distance=$(t_dist^(1/r)) after $(nIterations) iterations ", states, probabilities)
end

"""
	PlotLattice(lt::Lattice,fig=1)

<<<<<<< HEAD
Returns a plot of a lattice. The arguments is only a lattice.
<<<<<<< HEAD
"""

<<<<<<< HEAD
>>>>>>> e9b1bc9cdc5c989ee6e99a1505eeecf47d22e288
function PlotLattice(lt::Lattice)
    pt = subplot2grid((1,length(lt.state)),(0,0),colspan = length(lt.state))
    pt.spines["top"].set_visible(false)                                           #remove the box at the top
    pt.spines["right"].set_visible(false)                                         #remove the box at the right
=======
# function PlotLattice(lt::Lattice)
#     pt = subplot2grid((1,length(lt.state)),(0,0),colspan = length(lt.state))
#     pt.spines["top"].set_visible(false)                                           #remove the box at the top
#     pt.spines["right"].set_visible(false)                                         #remove the box at the right
#     for t = 2:length(lt.state)
#         for i=1:length(lt.state[t-1])
#             for j=1:length(lt.state[t])
#                 pt.plot([t-2,t-1],[lt.state[t-1][i],lt.state[t][j]])
#             end
#         end
#     end
#     xlabel("stage")
#     ylabel("states")
#     xticks(0:length(lt.state)-1)
# end
					
=======
"""					
>>>>>>> master
=======
Returns a plot of a lattice.
"""
<<<<<<< HEAD
>>>>>>> master
function PlotLattice(lt::Lattice,fig = 1)
=======
function plot_lattice(lt::Lattice,fig = 1)
>>>>>>> master
    if !isempty(fig)
        figure(figsize=(6,4))
    end
    lts = subplot2grid((1,4),(0,0),colspan = 3)
    title("states")
    xlabel("stage,time",fontsize=11)
    #ylabel("states")
    xticks(1:length(lt.state))
    lts.spines["top"].set_visible(false)                                                         # remove the box at the top
    lts.spines["right"].set_visible(false)                                                       # remove the box at the right
>>>>>>> master
    for t = 2:length(lt.state)
        for i=1:length(lt.state[t-1])
            for j=1:length(lt.state[t])
                lts.plot([t-1,t],[lt.state[t-1][i],lt.state[t][j]])
            end
        end
    end

    prs = subplot2grid((1,4), (0,3))
    title("probabilities")
    prs.spines["top"].set_visible(false)
    prs.spines["left"].set_visible(false)
    prs.spines["right"].set_visible(false)
# Use the states and probabilites at the last stage to plot the marginal distribution
    stts = lt.state[end]
    n = length(stts)                                    # length of terminal nodes of the lattice.
    h = 1.05 * std(stts) / (n^0.2) + eps()                  #Silverman rule of thumb
    #lts.set_ylim(minimum(stts)-1.5*h, maximum(stts)+1.5*h)
    #prs.set_ylim(minimum(stts)-3.0*h, maximum(stts)+3.0*h)
    proba = sum(lt.probability[end], dims=1)
    yticks(())                                          #remove the ticks on probability plot
    t = LinRange(minimum(stts) - h, maximum(stts) + h, 100) #100 points on probability plot
    density = zero(t)
    for (i, ti) in enumerate(t)
        for (j, xj) in enumerate(stts)
            tmp = (xj - ti) / h
            density[i] += proba[j]* 35/32 * max(1.0 -tmp^2, 0.0)^3 / h #triweight kernel
        end
    end
    plot(density, t)
    prs.fill_betweenx(t, 0 , density)
end
