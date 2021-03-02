using LinearAlgebra
using Clp
using MathProgBase
include("ontoSimplex.jl")


function distFunction(states1::Vector{Float64}, states2::Vector{Float64})::Array{Float64,2}
    n1= length(states1); n2= length(states2)
    dMatrix= Array{Float64}(undef, (n1,n2))
    for i= 1:n1, j= 1:n2
		dMatrix[i,j]= abs(states2[j]- states1[i])
    end
    return dMatrix
end


#	computes the Wasserstein distance
function Wasserstein(p1::Vector{Float64}, p2::Vector{Float64}, distMatrix::Array{Float64,2}, rWasserstein::Float64=1.)
	ontoSimplex!(p1); ontoSimplex!(p2)
    n1= length(p1);   n2= length(p2)

    A= kron(ones(n2)', Matrix{Float64}(I, n1, n1))
    B= kron(Matrix{Float64}(I, n2, n2), ones(n1)')

    x= linprog(vec(distMatrix.^rWasserstein), [A;B], '=', [p1;p2], ClpSolver())
    return (distance= (x.objval)^(1/ rWasserstein), π= reshape(x.sol, (n1, n2)))
end


#	Sinkhorn-Knopp iteration algorithm
function Sinkhorn(p1::Vector{Float64}, p2::Vector{Float64}, distMatrix::Array{Float64,2}, rWasserstein::Float64= 1., λ::Float64= 1.)
	ontoSimplex!(p1); ontoSimplex!(p2)
	r= p1; c= p2
	K= exp.(-λ * (distMatrix.^ rWasserstein))
	for i= 1:1000		# Sinkhorn iteration
		c= c./ (p2'*c)		# rescale
		r= p1./ (K* c)		# vector operation
		c= p2./ (K'* r)		# vector operation
	end
#	println("r=", r, p1'*r); println("c=", c, p2'*c)
	π= Diagonal(r)*K*Diagonal(c)
	return (distance= (sum(π.* (distMatrix.^rWasserstein))) ^(1/rWasserstein), π= π)
end


# #	Sinkhorn-Knopp algorithm
# function Sinkhorn(p1::Vector{Float64}, p2::Vector{Float64}, distMatrix::Array{Float64,2}, rWasserstein::Float64= 1., λ::Float64= 1.)
# 	ontoSimplex!(p1); ontoSimplex!(p2)
# 	π= exp.(-λ * (distMatrix.^ rWasserstein))
# 	for i= 1:100
# 		π.*= (p2'./ sum(π, dims=1)) # scale to sum of column= 1
# 		π.*= (p1 ./ sum(π, dims=2)) # scale to sum of lines= 1
# 	end
# 	return (distance= sum(π.* distMatrix)^(1/ rWasserstein), π= π)
# end
