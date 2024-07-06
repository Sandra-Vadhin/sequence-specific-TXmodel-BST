function loss_scalar_N(κ::Array{Float64,1}, Y::Array{Float64,1}, T, Xsim::Array{Float64,2}, gidx)

    # κ map
    tmp_alpha = κ[1:24]
    tmp_g = κ[25:end]
   
    # set new parameters -
    α = model.α
    α[1:24] = tmp_alpha
    model.α = α

    # set G values -
    G = model.G;

   #Gs -
      for i in 1:length(tmp_g)
         model.G[gidx[i]] = tmp_g[i]
     end
 
 # get the G matrix -
 
 # run the model -
  U_data = Xsim[:,152]
 Yₘ = (U_data) # model_output_vector(T, CF) # properties of the CF curve 

 # check to see if solution length matches data length
 if length(Yₘ[:,1]) != 1601
     ϵ = Inf64
        
 else
     # data_tps = [0.0, 2.0, 4.0, 8.0, 16.0]
     data_vals = Y
     global sim_itp = zeros(size(Y))
     datapoints = [1, 201, 401, 801, 1601]
     # Interpolate simulation array to fit dataset
     for i in 1:length(1)
     sim = Yₘ[:,i]
     itp = interpolate((T,),sim,Gridded(Linear()))
     global sim_itp[1,i] = sim[1]
     global sim_itp[2,i] = sim[201]
     global sim_itp[3,i] = sim[401]
     global sim_itp[4,i] = sim[801]
     global sim_itp[5,i] = sim[1601]
     global sim_itp = sim_itp
     end
     # Calculate error
     norm_diff = (sim_itp.-data_vals)./maximum(data_vals,dims=1)
    weight = ones(1,1)
    # weight[9] = 100
    # weight[2] = 4
    # weight[30] = 8
    # weight[31] = 2
    # weight[26] = 10
    #= weight[1,62] = 10
    weight[1,63] = 10   =# 
    global ϵ = sum((norm_diff.^2).*weight)
    # global ϵ = norm((norm_diff.^2).*weight)
    # global ϵ = norm((Yₘ.-data_vals)./maximum(data_vals,dims=1))
    # global ϵ = norm((Yₘ.-data_vals)./maximum(data_vals,dims=1))
         # ϵ = norm(((Y .- Yₘ).^2)/Yₘ[argmax(Yₘ)])
     
 end

 # @info ϵ

 # return -
 return ϵ
end


function loss_scalar(κ::Array{Float64,1}, Y::Array{Float64,2},  model::BSTModel,  index::Int, gidx)

       # κ map
       tmp_alpha = κ[1:24]
       tmp_g = κ[25:end]
      
       # set new parameters -
       α = model.α
       α = tmp_alpha
     #=   α[1:2] = tmp_alpha[1:2]
       α[4:28] = tmp_alpha[3:27] =#
       model.α = α
   
       # set G values -
       G = model.G;
   
      #Gs -
         for i in 1:length(tmp_g)
           G[gidx[i]] = tmp_g[i]
        end
    model.G = G
    # get the G matrix -
    
    # run the model -
    (T,U) = evaluate_opt(model, tspan = (0.0, 16.0))
    U_data = U[:,txtl_idx]
    Yₘ = (U_data) # model_output_vector(T, CF) # properties of the CF curve 

    # check to see if solution length matches data length
    if length(Yₘ[:,1]) != 1601
        ϵ = Inf64
           
    else
        data_tps = [0.0, 2.0, 4.0, 8.0, 16.0]
        data_vals = Y
        global sim_itp = zeros(size(Y))
        datapoints = [1, 201, 401, 801, 1601]
        # Interpolate simulation array to fit dataset
        for i in 1:length(txtl_idx)
        sim = Yₘ[:,i]
        # itp = interpolate((T,),sim,Gridded(Linear()))
        global sim_itp[1,i] = sim[1]
        global sim_itp[2,i] = sim[201]
        global sim_itp[3,i] = sim[401]
        global sim_itp[4,i] = sim[801]
        global sim_itp[5,i] = sim[1601]
        global sim_itp = sim_itp 
        end
        # Calculate error
        norm_diff = (sim_itp.-data_vals)./maximum(data_vals,dims=1)
        global ϵ = sum(norm_diff.^2)
        # ϵ = sum(norm_diff.^2)
        # ϵ = norm(((Y .- Yₘ).^2)/Yₘ[argmax(Yₘ)])
        
    end

    # @info ϵ

    # return -
    return ϵ
end

function loss(κ::Array{Float64,1}, Y::Array{Float64,1},  model::BSTModel)

    # κ map
    # 1 - 9 : α vector
    model.α = κ[1:223]

    #Gs -

    # get the G matrix -
    G = model.G

    ## get nonzero G idxs for dynamic
    # G of nonzero idxs set to K mapping

    # run the model -
    (T,U) = evaluate(model, tspan = (0.0, 16.0))
    sol_w_data = U[:,data_idx_all]
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,2],sfa[1],FIIa,FII[index])
    data = [T U]
    Xₘ = hcat(data,CF)
    Yₘ = model_output_vector(T, CF) # properties of the CF curve 
    sum(abs2, X .- all_itp)
    ϵ = norm((Y .- Yₘ).^2)

    # return -
    return ϵ
end