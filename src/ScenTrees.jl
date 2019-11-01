module ScenTrees
include("TreeStructure.jl")
include("TreeApprox.jl")
include("StochPaths.jl")
include("LatticeApprox.jl")
<<<<<<< HEAD
export TreeApproximation!,LatticeApproximation,Tree,Lattice,stage,height,leaves,nodes,root,
        partTree,buildProb!,treeplot,plotD,PlotLattice,GaussianSamplePath,RunningMaximum


"""
	GaussianSamplePath()

Returns samples from the Gaussian Random Walk
"""

"""
	RunningMaximum()

Returns the maximum of consecutive numbers from the Gaussian distribution
"""

=======
include("KernelDensityEstimation.jl")
include("bushinessNesDistance.jl")
export TreeApproximation!,LatticeApproximation,Tree,Lattice,stage,height,leaves,nodes,root,
        partTree,buildProb!,treeplot,plotD,PlotLattice,GaussianSamplePath1D,GaussianSamplePath2D,
        RunningMaximum1D,RunningMaximum2D,path,bushinessNesDistance,LogisticKernel,KernelScenarios
>>>>>>> e9b1bc9cdc5c989ee6e99a1505eeecf47d22e288
end
