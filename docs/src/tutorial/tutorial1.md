<<<<<<< HEAD
# Introduction

Multistage stochastic optimization involves approximating a stochastic process in by a finite structure.
We call these structures a scenario tree and a scanrio lattice.
Scenario lattices plays an important role in approximating Markovian data processes.

We will look at the basic functionality of `ScenTrees.jl` just to highlight on how it works.

## Goal

The goal is to generate a valuated probability tree or a lattice which represents the stochastic process in the best way possible.

To measure the quality of the approximation, we use a process distance between the stochastic process and the scanario tree or lattice.

## Description of a scenario tree

A scenario tree is characterized by the following:
=======
```@meta
CurrentModule = ScenTrees
```

# Introduction

A stochastic program is a mathematical program that involves some uncertain data. These parameters may be mostly accurately described by random variables. In most cases, it is difficult to optimize directly in terms of the distributions of these random variables. Hence, in most cases, these distributions are approximated by discrete distributions with a finite number of scenarios for the random variables. This discretization procedure is what is often called `scenario generation`. Uncertainty in long-term capacity planning is inescapable. The random parameters can be conceived to follow a multistage stochastic process over some time space so that the discrete scenarios represent sample paths. The approach we take is to form an approximation of the original stochastic process by discretization.

In multistage stochastic optimization, we are interested in approximations of stochastic processes by finite structures. These processes are random and they have uncertain scenarios and a decision maker needs to make decisions at different stages of the process. It is useful to depict the possible sequences of data for this processes in form of a `scenario tree` in the case of a discrete time stochastic process and a `scenario lattice` for Markovian data processes.

A scenario tree is a set of nodes and branches used in models of decision making under uncertainty. Every node in the tree represents a possible state of the world at a particular point in time and a position where a decision can be made. Each tree node has a single predecessor and multiple successors whereas a lattice can have many predecessor.

A scenario tree/lattice is organized in levels which corresponds to stages ``1,\ldots,T``. Each node in a stage has a specified number of predecessors as defined by the branching structure. A node represents a possible state of the stochastic process and the vertices represents the possibility of transition between the two connected nodes. A scenario tree differs from a scenario lattice by the condition that each node in stage ``t`` must have one predecessor in stage ``t-1``. For a lattice, that is not the case; all the nodes in stage ``t-1`` share the same children in stage ``t``.

## Goal of `ScenTrees.jl`

We model stochastic processes by scenario trees and scenario lattices. The distributions of these processes may be continuous and involves parameters that are uncertain.

The goal  of `ScenTrees.jl` is to approximate the distributions of these stochastic processes by discrete distributions with finite number of scenarios of the random variables. We generate a valuated probability scenario tree or a scenario lattice which represents the stochastic process in the best way possible using the stochastic approximation algorithm. These processes are random and represent uncertainty at a particular state and at a certain point in time.  

These approximations should be tractable, which is small enough to allow for reasonable calculation times, but is large enough to capture the important features of the problem. We use the concept of multistage distance to determine the quality of the approximations

### Introductory example

Consider a simple Gaussian random walk in 5 stages. The starting value of this process is known and fixed, say at ``0`` and the other values are random. The following plot shows 100 sample paths of this process:

![100 sample paths from Gaussian random walk](../assets/100GaussianPaths.png)

We generate and improve a scenario tree or a scenario lattice using this stochastic process. The number of iterations for the algorithm equals the number of sample paths that we want to generate from the stochastic process. Also, the number of stages in the stochastic process equals the number of stages in the scenario tree or the scenario lattice.

The user is free to choose any branching structure for the scenario tree/lattice. The branching structure shows how many branches each node in the tree has at each stage of the tree. For example, we can use a branching structure of ``1x2x2x2x2`` for the scenario tree. This means that each node in the tree has two children. Basically, this is a `binary tree`. It has been shown that the elements in the branching structure have a direct relationship with the quality of the resulting scenario tree/lattice. A scenario tree/lattice with many branches has a better approximation quality than a scenario tree with less branches.

Using the binary branching structure stated above, we obtain the following valuated probability tree that represents the above stochastic process:

![Scenario Tree 1x2x2x2x2](../assets/TreeExample.png)

*Figure 1: Scenario Tree 1x2x2x2x2*

The above tree is optimal and therefore can be used by a decision maker for a decision making process depending on the type of problem he/she is handling. To measure the quality of this approximation, we use the concept of multistage distance between the stochastic process and the scenario tree or lattice, which we introduce in the following subsection.

### Multistage distance

To measure the distance of stochastic processes, it is not sufficient to only consider the distance between their laws. It is also important to consider the information accumulated over time i.e., what the filtration has to tell us over time. The Wasserstein distance do not correctly separate stochastic processes having different filtration. It ignores filtration and hence does not distinguish stochastic processes. Multistage distance comes in handy in the situations for measuring distances for stochastic processes. Multistage distance is also called the `process distance` or `nested distance`.

Multistage distance was introduced by [Georg Ch. Pflug (2009)](https://doi.org/10.1137/080718401). It turns out that this distance is very important to measure the distance between multistage stochastic processes as it incorporates filtration introduced by the processes. We use this distance in our algorithm to measure the quality of approximation of the scenario tree and scenario lattice. Generally, a scenario tree/lattice with a minimal distance to the stochastic process is consider to have a better quality approximation.

The distance between the above scenario tree and the original process is `0.0894`. This shows that the scenario tree above approximates the stochastic process well. This tree can therefore be used for decision making under uncertainty.

## Description of a scenario tree

A scenario tree is described by the following:
>>>>>>> e9b1bc9cdc5c989ee6e99a1505eeecf47d22e288

1. Name of the tree
2. Parents of the nodes in the tree
3. Children of the parents in the tree
4. States of the nodes in the tree
5. Probabilities of transition from one node to another.

<<<<<<< HEAD
<<<<<<< HEAD
The tree is also characterized by its _nodes_, _stages_, _height_,_leaves_ and the _root_ of the tree or the nodes.

Each tree has stages starting from `0` where the root node is. 
=======
A scenario tree is a mutable struct of type `Tree()`. To create a non-optimal scenario tree, we need to fix the branching structure and the dimension of the states of nodes you are wroking on. This typs `Tree()` has different methods:
=======
A scenario tree is a mutable struct of type `Tree()`. To create a non-optimal scenario tree, we need to fix the branching structure and the dimension of the states of nodes you are working on. The type `Tree()` has different methods:
>>>>>>> master
```julia
julia> using Pkg
julia> Pkg.add("ScenTrees")
julia> using ScenTrees
julia> methods(Tree)
# 4 methods for generic function "(::Type)":
[1] Tree(name::String, parent::Array{Int64,1},
children::Array{Array{Int64,1},1}, state::Array{Float64,2}, probability::Array{Float64,2})
[2] Tree(identifier::Int64)
[3] Tree(spec::Array{Int64,1})
[4] Tree(spec::Array{Int64,1}, dimension)
```
<<<<<<< HEAD
<<<<<<< HEAD
All the methods correspond to the way you can create a scenario tree. For the first method, the length of states must be equal to the length of the probabilities. In the 2nd method, you can call any of our predefined trees by just calling on the identifier (these identifiers are `0,301,302,303,304,305,306,307,401,402,4022,404,405`). And finaly the most important methods are the 3rd and 4th method. If you know the branching structure of your scenario tree, then you can create an non-optimal starting tree using it. If you don't state the dimension you ae working on, then it is defaulted into `1`. For example, `Tree([1,2,2,2,2])` creates a binary tree with states of dimension one as in Figure 1 above
>>>>>>> e9b1bc9cdc5c989ee6e99a1505eeecf47d22e288
=======
All the methods correspond to the way you can create a scenario tree. For the first method, the length of states must be equal to the length of the probabilities. In the 2nd method, you can call any of our predefined trees by just calling on the identifier (these identifiers are `0,301,302,303,304,305,306,307,401,402,4022,404,405`). And finally the most important methods are the 3rd and 4th method. If you know the branching structure of your scenario tree, then you can create an non-optimal starting tree using it. If you don't state the dimension you are working on, then it is defaulted into `1`. For example, `Tree([1,2,2,2,2])` creates a binary tree with states of dimension one as in Figure 1 above
>>>>>>> master
=======
All the methods correspond to the way you can create a scenario tree. For the first method, the length of states must be equal to the length of the probabilities. In the 2nd method, you can call any of our predefined trees by just calling on the identifier (these identifiers are `0, 301, 302, 303, 304, 305, 306, 307, 401, 402, 4022, 404, 405`). And finally the most important methods are the 3rd and 4th method. If you know the branching structure of your scenario tree, then you can create an non-optimal starting tree using it. If you don't state the dimension you are working on, then it is defaulted into `1`. For example, `Tree([1,2,2,2,2])` creates a binary tree with states of dimension one as in Figure 1 above
>>>>>>> master

## Description of a scenario lattice

A scenario lattice differs from a scenario tree in that every node in stage `t` is a child for each node in stage `t-1`. So the nodes in stage `t-1` share the same children.

Due to the above, we only describe a scenario lattice by:

<<<<<<< HEAD
1. Name of the lattice 
2. States of the nodes of the lattice
3. Probabilities of transition from one node to another in the lattice

# Usage

Since we have the basics of the scenario tree and the scenario lattice and since we created `ScenTrees.jl` 
with an intention of being user-friendly, we will give an example of its usage and explain each part of it.

```julia
using ScenTrees

ex1 = Tree([1,3,3,3,3]);
sol1 = TreeApproximation!(ex1, GaussianSamplePath, 100000,2,2);

treeplot(sol1)
```
In the above, we are creating a scenario tree with the branching structure `1x3x3x3x3` as `ex1`. 
We want to approximate the Gaussian random walk with this tree.
So, the tree approximation process takes 4 inputs:

1. A tree
2. The process to be approximated
3. Sample size
4. Which norm (default to Euclidean norm = 2)
5. The value of `r` for process distance ( default to `r=2`)

Those are basically the inputs that we are giving to the function _TreeApproximation!_. 
What we obtain as a result of this is a valuated tree which can be visualized by the function _treeplot_ or _plotD_ if incase you were dealing with a 2-dimensional state space.

![Stochastic approximation](../assets/exampleTree1.png)

Just as easy as that! 

It is more even simple for lattice approximation as in this we only take as inputs the branching structure of the lattice and the sample size as follows:

```julia
sol2 = LatticeApproximation([1,2,3,4,5], 500000)

PlotLattice(sol2)
```

# Plotting

As shown in the above, there are three different functions for plotting. We have _treeplot_, _plotD_ and _PlotLattice_ functions. 
Each of these functions is special in its own way. Both _treeplot_ and _plotD_ are for plotting scenario trees while _PlotLattice_ is only for plotting lattice.
In this package, we deal also with trees of 2D state space. So you can visualize them after approximation using _plotD_ function. 
So the main difference between _treeplot_ and _plotD_ is that  _treeplot_ is for trees only in 1D state space while _plotD_ can be used for trees in 1D and 2D but specifically created for trees in 2D state space.

!!! info
	You need to install the [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl) package for this plots.

You can save the plots using the the `Plots.jl` function `savefig`:
```julia
Plots.savefig("example1.pdf")
```

This ends the tutorial for this package.
=======
1. Name of the lattice
2. States of the nodes of the lattice
3. Probabilities of transition from one node to another in the lattice

A scenario lattice has only one method.
```julia
julia> methods(Lattice)
 1 method for generic function "(::Type)":
[1] Lattice(name::String, state::Array{Array{Float64,2},1},
probability::Array{Array{Float64,2},1})
```
This method is not very important because we only need it to produce the results of the lattice approximation process. We will see later that for lattice approximation, we need the branching structure and so the structure of the lattice is not very important as in the case of a scenario tree.

## Exported functions

Since we have the basics of the scenario tree and the scenario lattice and since we created `ScenTrees.jl` with an intention of being user-friendly, we present the exported functions that are visible to the user i.e., that are public, and the user can call these functions depending on what he/she wants to achieve with this package:

1. Tree (associated are: nodes, stage, height, leaves, root, part_tree, build_probabilities!),
2. tree_approximation!
3. Lattice,
4. lattice_approximation,
5. kernel_scenarios (for conditional density estimation method)
6. Plotting utilities (these functions include: tree_plot, plot_hd and plot_lattice),
7. Examples of process functions (gaussian_path1D, gaussian_path2D, running_maximum1D, running_maximum2D, path) and,
8. bushiness_nesdistance (returns a graph showing how different factors affects the multistage distance.)

- The most important functions in this module are `tree_approximation!()` and `lattice_approximation()` since these are the two functions which are used to approximate scenario trees and scenario lattices respectively.

- The other important function is the `Tree(bstructure, dimension)` function which gives the basic starting structure of a scenario tree.

### Querying the documentation of each function

All of the above functions have been documented in their respective scripts and the user can find out what each function does by putting a `?` before the function. For example, `?leaves` will give an explanation of what the function `leaves` does.

<<<<<<< HEAD
<<<<<<< HEAD
In the upcoming tutorials, we will have a look in detail on the functionalities of the main functions of this library.
>>>>>>> e9b1bc9cdc5c989ee6e99a1505eeecf47d22e288
=======
In the upcoming tutorials, we will have a look in detail on the functionalities of the main functions of this package.
>>>>>>> master
=======
In the upcoming tutorials, we will have a look in detail on what each function of this package does.
>>>>>>> master
