include("Include.jl")
include("DataProcessing.jl")
include("OGrxns.jl")
include("getalphas.jl")

model_file = joinpath(pwd(),"txonly.bst")
model = build(model_file)

# load the training data -

# which case? "control", "tta", or "dnp"
case = "control"

t_array = collect(0:0.01:16.0)

# interpolate data, grab indices, and get initial concentrations
processed_data = ProcessData(model,case,t_array)

dynamic_idx = processed_data[:Species_index_metabolite]
data_species_idx = processed_data[:data_Species_index_metabolite]
aa_idx = processed_data[:Species_index_aa]
itp_met = processed_data[:Conc_met]
itp_aa = processed_data[:Conc_aa]
itp_mrna = processed_data[:Conc_mrna]
itp_protein = processed_data[:Conc_protein]
energyidx = processed_data[:index_energy]
glycolysis_idx = processed_data[:index_gly]
reducing_idx = processed_data[:index_reducing]
TCA_idx = processed_data[:index_tca]
enzA_idx = processed_data[:index_enzyme_activity]
Vmax2 = processed_data[:enzyme_vmax_2h]
Vmax8 = processed_data[:enzyme_vmax_8h]
met_aa_data_idx = processed_data[:data_idx]
ica = processed_data[:init_data]
enzIDx = processed_data[:enz_idx]
sfav = processed_data[:sfa_vals]
sats = processed_data[:sat]
rxn = processed_data[:rxn]
Km = processed_data[:Km]
data2control = processed_data[:data2c]
all_raw_data = processed_data[:raw_data]
data_idx_all = processed_data[:all_data_idx]

txtl_idx = zeros(Int64,31)
txtl_idx[1:9] = [24,31,38,41,60,71,72,141,142]
txtl_idx[10:end-2] = Int.(vec(aa_idx))
txtl_idx[end-1] = 152
txtl_idx[31] = 148

tx_idx = [24,31,38,41,71,72,141,142,152]
tx_data_idx = zeros(Int64,9)
for i in 1:9
    tx_data_idx[i] = findfirst(x->x==tx_idx[i],data_idx_all)
end


txtl_data_idx = zeros(Int,31)
for i in 1:31
    txtl_data_idx[i] =findfirst(x->x==txtl_idx[i],data_idx_all)
end


all_raw_data = Matrix(all_raw_data)
all_raw_data[:,end-1] = all_raw_data[:,end-1]./1e6
all_raw_data[:,end] = all_raw_data[:,end]./1e3
txtl_raw = all_raw_data[:,txtl_data_idx]
tx_data = all_raw_data[:,tx_data_idx]


all_itp = hcat(itp_met,itp_aa,itp_mrna,itp_protein)

tx_itp = all_itp[:,tx_data_idx]

ica[1] = 5e-6
ica[152] = 0
 ica[147] = 0
 model.initial_condition_array = ica
 model.static_factors_array = sfav

 # estimate orders and rates based on enzymatic rate laws
y = zeros(200)
b = zeros(157)

result = DiffResults.JacobianResult(y,b)
result = ForwardDiff.jacobian!(result,ogrxnsetup,ica[1:157])
jac = DiffResults.jacobian(result)
vel = DiffResults.value(result)
J = Matrix{Float64}(transpose(jac))
nz = findall(x->x!=0&&!isnan(x),J)
gg1 = zeros(157,200)

for i in 1:200
   for j in 1:157
      gg1[j,i] = J[j,i]*ica[j]/vel[i]
   end
end

ag1 = getalphas(ica[1:157],vel,gg1)

 G = model.G
 # set some G's close to zero - very saturated
h2o_G_idx = findall(x->x!=0,model.G[74,1:4])
G[74,h2o_G_idx] .= 0.001
# h__idx = findall(x->x!=0,model.G[76,1:25])
# G[76,h__idx] .= 0.0
# o2_G_idx = findall(x->x!=0,model.G[109,1:25])
# G[109,o2_G_idx] .= 0.001
# pi_G_idx = findall(x->x!=0,model.G[116,1:25])
# G[116,pi_G_idx] .= 0.0
# tRNA_G_idx = findall(x->x!=0,model.G[153,:])
# G[153,tRNA_G_idx] .= 1.0
 G[152,3] = 100.0
# G[findall(x->x>=10.0,sats)] .= 0.0
model.G = G
model.α .= 1.0
model.α[1] = ag1[176]
model.α[2] = ag1[177]
model.α[3] = 2.38
model.α[4] = ag1[70]
model.α[5] = 100.0

 gidx = findall(x->x!=0&&x!=100.0, model.G[1:158,1:4])
#=  gidx2 = findall(x->x!=0, model.G[1:158,4:26])
 gidx = vcat(gidx1,gidx2) =#
 
  # Create Ensemble_TX_MRNA_3raw directory if it doesn't already exist
 if ~isdir("Ensemble_TX_MRNA_3raw")
     mkdir("Ensemble_TX_MRNA_3raw")
 end
 
 # Create Best directory if it doesn't already exist
 if ~isdir("Ensemble_TX_MRNA_3raw/Best")
     mkdir("Ensemble_TX_MRNA_3raw/Best")
 end
   
 # Create or read the relevant information
   if ~(isfile("Ensemble_TX_MRNA_3raw/Best/rate_constant.dat")&&isfile("Ensemble_TX_MRNA_3raw/Best/saturation_constant.dat"))
     # Take parameters to be varied from data_dictionary_best
 
 
     model.initial_condition_array = ica
     model.static_factors_array = sfav
     model.G[152,3] = 1.0
     g_length = length(Tuple.(gidx))  

     g = model.G[gidx] 
     a = model.α[1:24]
     global parameters = vcat(a,g)
     
         global rate_best = a
         global rate = copy(rate_best)
         global sat_best = g
         global sat = copy(sat_best)
        #=  global cont_best = data_dictionary_best[:CONTROL_PARAMETER_ARRAY]
         global cont = copy(cont_best) =#
         # Solve the balance equations and calculate cost
     
         Tsim_best,X_best = evaluate_opt(model,tspan = (0.0,16.0))
         Error_best = loss_scalar_N(parameters, all_raw_data[:,end-1], Tsim_best,X_best,gidx)
         global cost_best = Error_best
         global cost = copy(cost_best)
         # Save the relevant information
         writedlm("Ensemble_TX_MRNA_3raw/Best/rate_constant.dat",rate)
         writedlm("Ensemble_TX_MRNA_3raw/Best/saturation_constant.dat",sat)
         # writedlm("Ensemble_TX_MRNA_3raw/Best/control_constant.dat",cont)
         writedlm("Ensemble_TX_MRNA_3raw/Best/Tsim",Tsim_best)
         writedlm("Ensemble_TX_MRNA_3raw/Best/X",X_best)
         writedlm("Ensemble_TX_MRNA_3raw/Best/cost",cost_best)
     else
         # Read the parameters from disk
         global rate_best = readdlm("Ensemble_TX_MRNA_3raw/Best/rate_constant.dat")
         global rate = copy(rate_best)
         global sat_best = readdlm("Ensemble_TX_MRNA_3raw/Best/saturation_constant.dat")
         global sat = copy(sat_best)
         # global cost_best = readdlm("Ensemble_TX_MRNA_3raw/Best/cost")
         # global cost_best = cost_best[1]
         # global cost = copy(cost_best)
         # global cont_best = readdlm("Ensemble_TX_MRNA_3raw/Best/control_constant.dat")
         # global cont = copy(cont_best)
         # Check that params is of correct length
      #=    if ~((num_rate==length(rate))&&(num_sat==length(findall(x->x>0,sat)))&&(num_cont==length(cont)))
             throw("Wrong number of parameters")
         end =#
         # Add parameters to data_dictionary_best
         
         model.α = vec(rate)
         model.G[gidx] = vec(sat)
         model.G[152,3] = 1.0
         
         g = model.G[gidx] 
         a = model.α[1:24]
         global parameters = vcat(a,g)
            # data_dictionary_best[:CONTROL_PARAMETER_ARRAY] = cont
         # Solve the balance equations and calculate cost
         Tsim_best,X_best = evaluate_opt(model,tspan = (0.0,16.0))
         Error_best = loss_scalar_N(parameters, all_raw_data[:,end-1], Tsim_best,X_best,gidx)
         global cost_best = Error_best
         global cost = copy(cost_best)
         # Save the relevant information
         writedlm("Ensemble_TX_MRNA_3raw/Best/Tsim",Tsim_best)
         writedlm("Ensemble_TX_MRNA_3raw/Best/X",X_best)
         writedlm("Ensemble_TX_MRNA_3raw/Best/cost",cost_best)
         end
     
     nparam = length(parameters)
     κ = Array{Vector,1}(undef,3)
     κ[1] = parameters
     κ[2] = zeros(nparam)
     k2 = κ[2]
     # k2[1] = 70000.0
     # k2[2] = 30000.0
     # k2[1:28] = parameters[1:28]
     k2[3] = 0.8*model.α[3]
     κ[2] = k2
     κ[3] = zeros(nparam)
     k3 = κ[3]
     k3[1:2] .= 200000.0
     k3[3] = 1.2*model.α[3]
     k3[4] = 600000.0
     k3[5] = 1000.0
     k3[6:24] .= 100.0
     # k3[1:28] = parameters[1:28]
     k3[25:end] .= 3.0
     # k3[35] = 0.0
     # k3[38:39] .= 0.0
     κ[3] = k3
     
     if ~isfile("Ensemble_TX_MRNA_3raw/num_dir")
         writedlm("Ensemble_TX_MRNA_3raw/num_dir",0)
     end
     
     num_dir = convert(Int64,readdlm("Ensemble_TX_MRNA_3raw/num_dir")[1]) # Number of pre-existing directories
     
     num_iter = 1000000
     num_to_reset = 100 # Number of iterations before resetting to best-fit set
     
     alpha = 10 # Increasing alpha decreases acceptance probability (for cost_new > cost)
     
     variance_lower_bound = .01
     variance_upper_bound = 0.2
      
      for i = num_dir+1:num_dir+num_iter
     
     variance = exp(log(variance_lower_bound)+rand()*(log(variance_upper_bound)-log(variance_lower_bound)))
         
         if rem(i,num_to_reset) == 0 # Reset to best
             global rate = copy(rate_best)
             global cost = copy(cost_best)
         end
         
         # Copy data_dictionary
         global rate = deepcopy(rate_best)
         global sat = deepcopy(sat_best)
         
         param = vec(vcat(rate,sat))
         
         # Create perturbation vectors
         perturb = exp.(variance*randn(length(param)))
         param_new = param.*perturb
         param_new = max.(κ[2],param_new)
         param_new = min.(κ[3],param_new)
        # param_new[2] = min(0.8*param_new[1],param_new[2])
         rate_new = param_new[1:24]
         sat_new = param_new[25:end]
     
          dd = deepcopy(model)
              
         dd.α = rate_new
         dd.G[gidx] = sat_new
         dd.initial_condition_array = ica
             
         # Solve the balance equations
          dodo = @timed        Tsim,X = evaluate_opt(dd,tspan = (0.0,16.0))
          
     
         # println("CAT: ",X[end,98])
         if dodo.time < 40 # Should prevent sets with CVODE errors from joining ensemble
                              # May also exclude some sets without CVODE errors
             Error_new = loss_scalar_N(param_new, all_raw_data[:,end-1], Tsim,X,gidx)
             cost_new = sum(Error_new)
             
             acc_prob = exp(alpha*(cost-cost_new)/cost)
             
             cost_round = round(cost,digits=2)
             cost_new_round = round(cost_new,digits=2)
             cost_best_round = round(cost_best,digits=2)
             acc_prob_round = round(acc_prob,digits=2)
             var_round = round(variance,digits=2)
             
             # If a new best is achieved, overwrite parameters and cost and save the relevant information to Best directory
             if cost_new < cost_best
                 # Save to Best directory
                 global rate_best = copy(rate_new)
                 global sat_best = copy(sat_new)
                 # global cont_best = copy(cont)
                 global cost_best = copy(cost_new)
                 writedlm("Ensemble_TX_MRNA_3raw/Best/rate_constant.dat",rate_best)
                 writedlm("Ensemble_TX_MRNA_3raw/Best/saturation_constant.dat",sat_best)
                 # writedlm("Ensemble_TX_MRNA_3raw/Best/control_constant.dat",cont_best)
                 writedlm("Ensemble_TX_MRNA_3raw/Best/Tsim",Tsim)
                 writedlm("Ensemble_TX_MRNA_3raw/Best/X",X)
                 writedlm("Ensemble_TX_MRNA_3raw/Best/cost",cost_best)
                 println("$i: cost_new = $cost_new_round, cost = $cost_round, best = $cost_best_round, acc_prob = $acc_prob_round, var = $var_round, NEW BEST")
             else
                 println("$i: cost_new = $cost_new_round, cost = $cost_round, best = $cost_best_round, acc_prob = $acc_prob_round, var = $var_round")
             end
             # If new cost is lower than previous cost, choose new params as reference point for perturbation
             if rand() < acc_prob || isnan(acc_prob) # Accept all better sets and some worse sets
                 # Create directory for next sample and save the relevant information
                 mkdir("Ensemble_TX_MRNA_3raw/$i")
                 writedlm("Ensemble_TX_MRNA_3raw/$i/rate_constant.dat",rate_new)
                 writedlm("Ensemble_TX_MRNA_3raw/$i/saturation_constant.dat",sat_new)
                 # writedlm("Ensemble_TX_MRNA_3raw/$i/control_constant.dat",cont)
                 writedlm("Ensemble_TX_MRNA_3raw/$i/Tsim",Tsim)
                 writedlm("Ensemble_TX_MRNA_3raw/$i/X",X)
                 writedlm("Ensemble_TX_MRNA_3raw/$i/cost",cost_new)
                 # Overwrite variables
                 global cost = copy(cost_new)
             end
             writedlm("Ensemble_TX_MRNA_3raw/num_dir",i) # Record new number of directories
         end
     end