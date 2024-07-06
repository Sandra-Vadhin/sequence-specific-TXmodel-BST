using DataFrames
using CSV
using DelimitedFiles
using SmoothingSplines
using PyCall
using BSTModelKit
#np = pyimport("numpy")
@pyimport numpy as np

function ProcessData(model,case,t_array)

    data_met = DataFrame(CSV.File("config/SpeciesDict/$(case)_metabolites.dat",drop=[:FAD]))[1:5,:]
    data_aa = DataFrame(CSV.File("config/SpeciesDict/$(case)_aa.dat"))[1:5,:]
    data_mrna_protein = DataFrame(CSV.File("config/SpeciesDict/$(case)_mRNA_protein.dat"))[1:5,:]
    # deletecols!(data_met,:FAD)

    no_mets = length(data_met[1,:])
    no_aa = length(data_aa[1,:])
    data_raw = hcat(data_met,data_aa,data_mrna_protein)

    #===============Interpolate data for Smoothing function==========#
    data_t = [0.0;2.0;4.0;8.0;16.0;]; # Time points of data
    t_vec = collect(0.0:0.01:16.0)

    # np.interp(time points where you want interpolation, given time vector, given data)
    sim_met = zeros(length(t_vec),no_mets)
    for i = 1:no_mets
        sim_met[:,i] = np.interp(t_vec,data_t,data_met[:,i])
    end

    sim_aa = zeros(length(t_vec),no_aa)
    for i = 1:no_aa
        sim_aa[:,i] = np.interp(t_vec,data_t,data_aa[:,i])
    end

    #sim_mrna = zeros(length(t_vec),1)
    sim_mrna = np.interp(t_vec,data_t,data_raw[:,end-1]./1e6)

    #sim_protein = zeros(length(t_vec),1)
    sim_protein = np.interp(t_vec,data_t,data_raw[:,end]./1e3)
    
    #===============Smooth concentration values==========#
    smooth_met = zeros(length(t_vec),no_mets)
    for i = 1:no_mets
        spl_met = fit(SmoothingSpline,t_vec,sim_met[:,i],100.0)
        y_met = SmoothingSplines.predict(spl_met)
        diff = y_met[1] - sim_met[1,i]
        smooth_met[:,i] = y_met .- diff
    end

    smooth_aa = zeros(length(t_vec),no_aa)
    for i = 1:no_aa
        spl_aa = fit(SmoothingSpline,t_vec,sim_aa[:,i],100.0)
        y_aa = SmoothingSplines.predict(spl_aa)
        diff = y_aa[1] - sim_aa[1,i]
        smooth_aa[:,i] = y_aa .- diff
    end

    #smooth_mrna = zeros(length(t_vec),1)
    
        spl_mrna = fit(SmoothingSpline,t_vec,sim_mrna,100.0)
        y_mrna = SmoothingSplines.predict(spl_mrna)
        diff = y_mrna[1] - sim_mrna[1]
        smooth_mrna = y_mrna .- diff

        spl_protein = fit(SmoothingSpline,t_vec,sim_protein,100.0)
        y_protein = SmoothingSplines.predict(spl_protein)
        diff = y_protein[1] - sim_protein[1]
        smooth_protein = y_protein .- diff
    

    #===============Interpolate Smooth Curves for Flux==========#
    #np.interp(time points where you want interpolation, given time vector, given data)
    
    sim_met_smooth = zeros(length(t_array),no_mets)
    for i = 1:no_mets
        sim_met_smooth[:,i] = np.interp(t_array,t_vec,smooth_met[:,i])
    end

    sim_aa_smooth = zeros(length(t_array),no_aa)
    for i = 1:no_aa
        sim_aa_smooth[:,i] = np.interp(t_array,t_vec,smooth_aa[:,i])
    end

     
        sim_mrna_smooth = np.interp(t_array,t_vec,smooth_mrna)
        sim_protein_smooth = np.interp(t_array,t_vec,smooth_protein)
    

    no_data_points = length(sim_met_smooth[:,1])
#=     #=============Extract Flux Values=========================#
    flux_met = zeros(no_data_points-1,no_mets)
    flux_aa = zeros(no_data_points-1,no_aa)
    for flux_idx = 1:(no_data_points-1)
        flux_met[flux_idx,:] .= (sim_met_smooth[flux_idx+1,:] - sim_met_smooth[flux_idx,:])/(t_array[flux_idx+1]-t_array[flux_idx])
        flux_aa[flux_idx,:] .= (sim_aa_smooth[flux_idx+1,:] - sim_aa_smooth[flux_idx,:])/(t_array[flux_idx+1]-t_array[flux_idx])
    end =#
    #================Find Species Index=======================#
    network_idx = DataFrame(CSV.File("config/Reactions.txt",header=[Symbol("idx"),Symbol("rxn_name"),Symbol("substrate"),Symbol("arrow"),Symbol("product")]))

    AA_name = readdlm("config/SpeciesDict/AA_name_internal.txt")
    Met_name  = readdlm("config/SpeciesDict/Met_name_internal.txt")
    dynSpecies_name = model.list_of_dynamic_species # readdlm("dynspecies.txt")
    Species_idx = model.total_species_list # readdlm("totalspecies.txt")
    enz_name = model.list_of_static_species

    names_met = names(data_met)
    names_aa = names(data_aa)
    names_data = vcat(names_met,names_aa)

    data_idx = zeros(Int64,length(names_data))

    for i in 1:length(data_idx)
        speciesname = names_data[i]
        data_idx[i] = findfirst(x->x==speciesname,model.list_of_dynamic_species)
    end
    #init_data = zeros(model.number_of_dynamic_states)

    all_data_idx = vcat(data_idx,152,148)


    init_data = readdlm("config/SpeciesDict/met_mM_ss.dat")
    for i in 1:length(data_idx)
        ica_idx = data_idx[i]
        if i < 42
        init_data[ica_idx] = data_met[1,i]
        elseif i > 41
        init_data[ica_idx] = data_aa[1,i-41]
        end
    end

    init_data[findfirst(x->x=="M_maltodextrin6_c",model.total_species_list)] = 3e-5
    init_data[findfirst(x->x=="RIBOSOME",model.total_species_list)] = 2.15e-3
    init_data[findfirst(x->x=="RNAP",model.total_species_list)] = 0.0675e-3
    init_data = vec(init_data)
    init_data = push!(init_data,35*1e-6)


    data2c = zeros(model.number_of_dynamic_states)
    for i in 1:length(data_idx)
        ica_idx = data_idx[i]
        if i < 42
        data2c[ica_idx] = data_met[2,i]
        elseif i > 41
        data2c[ica_idx] = data_aa[2,i-41]
        end
    end 
    data2c[152] = 575.2173316e-6
    data2c[148] = 6.284153543e-3


    dataMet_Idx = zeros(no_mets,1)
    dataMetIdx = convert(Array{Int64,2},dataMet_Idx)
    for i = 1:no_mets
        dataMetIdx[i,1] = findall(Met_name[i].==Species_idx[:,1])[1]
    end

    AA_Idx = zeros(no_aa,1)
    AAIdx = convert(Array{Int64,2},AA_Idx)
    for i = 1:no_aa
        AAIdx[i,1] = findall(AA_name[i].==Species_idx[:,1])[1]
    end

    Met_Idx = zeros(length(dynSpecies_name),1)
    MetIdx = convert(Array{Int64,2},Met_Idx)
    for i = 1:length(Met_Idx)
        MetIdx[i,1] = findall(dynSpecies_name[i].==Species_idx[:,1])[1]
    end

    enz_idx = zeros(Int64,length(enz_name))
    for i = 1:length(enz_name)
        enz_idx[i] = findall(enz_name[i].==Species_idx[:,1])[1]
    end

    enz_data = readdlm("enzyme_name_conc_nM_kcat_1_over_s_rxn_no.dat")

    enzdata_name = enz_data[1,:]
    E = enz_data[2,:]./1e6
    E_idx = findall(E.<2.4e-5) #Find all enzymes that were not reported
    E[E_idx] .= 5e-5 #Set all unknown enzymes to 50nM

    sfa_vals = zeros(length(model.static_factors_array))
    
    for i in 1:length(sfa_vals)
        enzname_i = enz_name[i]
        enz_model_idx = findfirst(x->x==enzname_i,enzdata_name)
        sfa_vals[i] = E[enz_model_idx]
    end

    #=====================Separate Metabolite Index into subsections==============================#
    #Energy
    energy = readdlm("config/SpeciesDict/energy_idx.dat")
    no_energy = length(energy)
    Energy_Idx = convert(Array{Int64,2},zeros(no_energy,1))
    for i = 1:no_energy
        Energy_Idx[i] = getindex(findall(energy[i].==Met_name))[1]
    end
    e_idx = sort(Energy_Idx, dims = 1)
    EnergyIdx = MetIdx[e_idx]
    # flux_energy = flux_met[:,e_idx]

    #Glycolysis
    glycolysis = readdlm("config/SpeciesDict/gly_index.dat")
    no_gly = length(glycolysis)
    Gly_Idx = convert(Array{Int64,2},zeros(no_gly,1))
    for i = 1:no_gly
        Gly_Idx[i] = getindex(findall(glycolysis[i].==Met_name))[1]
    end
    g_idx = sort(Gly_Idx, dims = 1)
    GlyIdx = MetIdx[g_idx]
    # flux_gly = flux_met[:,g_idx]

    #Reducing
    reducing = readdlm("config/SpeciesDict/reducing_index.dat")
    no_reducing = length(reducing)
    Reducing_Idx = convert(Array{Int64,2},zeros(no_reducing,1))
    for i = 1:no_reducing
        Reducing_Idx[i] = getindex(findall(reducing[i].==Met_name))[1]
    end
    r_idx = sort(Reducing_Idx, dims = 1)
    ReducingIdx = MetIdx[r_idx]
    # flux_reducing = flux_met[:,r_idx]

    #TCA
    tca = readdlm("config/SpeciesDict/tca_index.dat")
    no_tca = length(tca)
    TCA_Idx = convert(Array{Int64,2},zeros(no_tca,1))
    for i = 1:no_tca
        TCA_Idx[i] = getindex(findall(tca[i].==Met_name))[1]
    end
    t_idx = sort(TCA_Idx, dims = 1)
    TCAIdx = MetIdx[t_idx]
    # flux_tca = flux_met[:,t_idx]

    #===========================Load in Enzyme Activity Assays =======================================#
    rxn_enzyme_data = convert(Array{Int64,2},readdlm("config/EnzymeDict/rxn_index_enzyme_data.dat"))
    vmax_2h = readdlm("config/EnzymeDict/$(case)_enzyme_vmax.dat")[:,1]
    vmax_8h = readdlm("config/EnzymeDict/$(case)_enzyme_vmax.dat")[:,2]
    

    Km = readdlm("config/KineticDict/Km_mM_matrix.dat");

sat = zeros(Float64,157,200)
for i in 1:length(model.list_of_dynamic_species)-1
  for j in 1:200
      s = init_data[i]
      km = Km[i,j]
      if km!=0
      sat[i,j] = s/(km)
      else 
        sat[i,j] = 0.0
  end
end
end

kcat = enz_data[3,:] .* 3660 #.*3600;
rxn = zeros(200)
for i = 1:200
    rxn[i] = kcat[i]
end

    #==================Store In Dictionary========================#
    ConcData_dictionary  = Dict{Symbol,Any}()
    ConcData_dictionary[:Species_index_metabolite] = MetIdx
    ConcData_dictionary[:data_Species_index_metabolite] = dataMetIdx
    ConcData_dictionary[:Species_index_aa] = AAIdx
    ConcData_dictionary[:Conc_met] = sim_met_smooth
    ConcData_dictionary[:Conc_aa] = sim_aa_smooth
    ConcData_dictionary[:Conc_mrna] = sim_mrna_smooth
    ConcData_dictionary[:Conc_protein] = sim_protein_smooth
    ConcData_dictionary[:index_energy] = EnergyIdx
    ConcData_dictionary[:index_gly] = GlyIdx
    ConcData_dictionary[:index_reducing] = ReducingIdx
    ConcData_dictionary[:index_tca] = TCAIdx
    ConcData_dictionary[:index_enzyme_activity] = rxn_enzyme_data
    ConcData_dictionary[:enzyme_vmax_2h] = vmax_2h
    ConcData_dictionary[:enzyme_vmax_8h] = vmax_8h
    ConcData_dictionary[:data_idx] = data_idx
    ConcData_dictionary[:init_data] = init_data
    ConcData_dictionary[:enz_idx] = enz_idx
    ConcData_dictionary[:sfa_vals] = sfa_vals
    ConcData_dictionary[:sat] = sat
    ConcData_dictionary[:rxn] = rxn
    ConcData_dictionary[:Km] = Km
    ConcData_dictionary[:data2c] = data2c
    ConcData_dictionary[:raw_data] = data_raw
    ConcData_dictionary[:all_data_idx] = all_data_idx

    return ConcData_dictionary
    #TXTL[:ConcDictionary] = ConcData_dictionary
  #  return TXTL
end
