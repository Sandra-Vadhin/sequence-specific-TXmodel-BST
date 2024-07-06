# it's training time -

function learn_optim(index::Int, model::BSTModel, training_data::Array{Float64,2}) #    pₒ::Union{Nothing,Array{Float64,1}} = nothing)

 # setup static -
 sfa = model.static_factors_array
 sfa = sfav
 # find DNP leak rxn and set to zero
 E_DNP_leak_idx = findfirst(x->x=="E_DNP_leak",model.list_of_static_species)
 sfa[E_DNP_leak_idx] = 0
 # return sfa
 model.static_factors_array = sfa
 
 # set up initial conditions array -
 # ica = ica from func then update with lit vals
 ℳ = model.number_of_dynamic_states
 xₒ = zeros(ℳ)
 xₒ = ica         
 model.initial_condition_array = xₒ
 
    # what is the output array?
    Y = training_data

G = model.G

# get nonzero gammas for dynamic species in non-transport rxns
# gidx = findall(x->x!=0, model.G[1:157,1:200])
gidx = findall(x->x!=0&&x!=100, model.G[1:158,1:26])
g_length = length(Tuple.(gidx))  
model.G[152,3] = 1.0

g = G[gidx]
a = model.α
# a = ones(223)
#= g = vec(readdlm("ensemble/gamma3.txt"))      # look at sample_ensemble.jl for specific G values
a = vec(readdlm("ensemble/alpha3.txt")) =#
#= a[203] = 0.0
a[213] = 0.0
a[214] = 0.0
 =#

 g = vec(readdlm("EnsembleTXTL_3_ITP_trna/Best/saturation_constant.dat"))
 a = vec(readdlm("EnsembleTXTL_3_ITP_trna/Best/rate_constant.dat"))

# fusion -
parameters = vcat(a,g)
# parameters = vcat(parameters[1:2],parameters[4:end])

# setup initial parameter values and bounds array -
nparam = length(parameters)

κ = Array{Vector,1}(undef,3)

κ[1] = parameters
κ[2] = zeros(nparam)
k2 = κ[2]
k2[3] = 0.5*2.38
κ[2] = k2
κ[3] = zeros(nparam)
k3 = κ[3]
k3[1:2] .= 100000.0
k3[3] = 2.0*2.38
k3[1:26] .= 500000.0
k3[3] = 2.0*model.α[3]
k3[27:90] .= 1000.0
#=  k3[28] = 100.0
k3[30] = 100.0
k3[31] = 100.0
k3[33:90] .= 100.0 =#
k3[91:end] .= 10.0
κ[3] = k3
 
    # set default set as the start -
 #=    if (isnothing(pₒ) == true)
        P = length(κ[1])
        σ = 0.1 # move up to 10%
        pₒ = κ[1].*(1 .+ σ*rand(-1:1,P))
    end =#

    # setup the objective function -
    inner_optimizer = NelderMead()
    obj_function(p) =  loss_scalar(p, Y, model,index,gidx) # sfa[1],## data) FII,index)
    results = optimize(obj_function, κ[2], κ[3], κ[1], Fminbox(inner_optimizer),
    Optim.Options(show_trace = true, show_every = 10, iterations=1))

#=      results = optimize(obj_function, κ[2], κ[3], κ[1], SAMIN(nt=1,ns=1,rt=0.9),
        Optim.Options(show_trace = true, show_every = 100, iterations=1000)) =#
        #= obj = (p,d) -> obj_function(p)

        f = OptimizationFunction(obj)
        prob = Optimization.OptimizationProblem(f, κ[1], lb = κ[2], ub = κ[3])
        sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 2) =#
#=         obj = (p,d) -> obj_function(p)

        f = OptimizationFunction(obj) #,Optimization.AutoZygote())

prob = Optimization.OptimizationProblem(f, κ[1], lb = κ[2], ub = κ[3])
sol = solve(prob, MultistartOptimization.TikTak(100), Optim.NelderMead(), Optim.Options(show_trace = true, show_every= 1, iterations = 5))
     =#

      
    
    # grab the best parameters -
    p_best = Optim.minimizer(results)

   # p_best = sol.u
    
    # run the sim w/the best parameters -
    model.α[1:90] = p_best[1:90]
   
    #= model.α[1:2] = p_best[1:2]
    model.α[4:28] = p_best[3:27] =#
    for i in 1:g_length-1
    model.G[gidx[i]] = p_best[i+90]
    end
  
    # run the model -
    (T,U) = evaluate_opt(model, tspan = (0.0, 16.0))
    #= CF = Array{Float64,1}
    CF = amplitude(T,U[:,2],sfa[1],FIIa,FII) =#
    U_data = U[:,txtl_idx]
    data = [T U]
    Xₘ = data
    Yₘ = (U_data) # model_output_vector(T, CF) # properties of the CF curve 
    Ymd = zeros(size(Y))
    Ymd[1,:] = Yₘ[1,:]
    Ymd[2,:] = Yₘ[201,:]
    Ymd[3,:] = Yₘ[401,:]
    Ymd[4,:] = Yₘ[801,:]
    Ymd[5,:] = Yₘ[1601,:]

    
    return (p_best, T, Xₘ, Yₘ, Ymd, Y)

#=     file = open("summ.txt", "w")

write(file, Optim.trace(results))
 =#
# close(f)
end