using DelimitedFiles
using Plots
using Statistics
using DataFrames
using CSV

num = readdlm("Ensemble_TX_MRNA_3raw/num_dir")
num = Int(num[1])
sim_array = zeros(1601,num)
simidx = zeros(num)
for i in 1:num
    if isfile("Ensemble_TX_MRNA_3raw/$i/X")
        sim_i = readdlm("Ensemble_TX_MRNA_3raw/$i/X")
        sim_i = sim_i[:,152]
        sim_array[:,i] = sim_i
        simidx[i] = i
    end
end

realsim = findall(x->x!=0,sim_array[1601,:])
sim_array = sim_array[:,realsim]
data = readdlm("config/SpeciesDict/_mRNA_protein.dat")
mean = data[2:6,1]
sd = data[7:11,1]


t_sim = collect(0.0:0.01:16)
case = "control"

t = [0;2;4;8;15.9]
Species_mean = Statistics.mean(sim_array, dims = 2)
Species_err = Statistics.std(sim_array, dims = 2)
Species_pos = Species_mean + 1.96*Species_err
Species_neg = Species_mean - 1.96*Species_err
    markercolor = "#1C66A8"
    shade = "powderblue"
    lcolor = markercolor
    linestyle = "-"

#=     plot(t, mean, yerr=sd,fmt="o",markersize = 3,capsize=2,elinewidth=1,color=markercolor)
    plot!(t_sim,Species_mean[:,1]*1e6,color=lcolor,linestyle=linestyle) =#
    #fill_between(t_sim,Species_pos[:,1]*1e6,Species_neg[:,1]*1e6,color=shade,alpha=0.5,linewidth=0)
    

fig1 = scatter(t, mean,yerr=sd, grid=false,xticks=0.0:1:16, yticks=0.0:100:1000, xlim = (0,16), ylim = (0,2),label = "",fmt="o",markersize = 3,capsize=2,color=markercolor, bg="aliceblue", background_color_outside="white", framestyle = :box, xlabel="Time (h)", ylabel="mRNA (nM)")
       # plot!(t_sim,Species_mean[:,1]*1e6, xticks=0.0:2:16, yticks=0.0:100:1000, xlim = (0,16), ylim = (0,1000),label="",lw = 1.25,c=colorant"#89CCE2", bg="aliceblue", background_color_outside="white", framestyle = :box, xlabel="Time (h)", ylabel="mRNA (nM)")
       plot!(t_sim,Species_mean[:,1],ribbon = (Species_pos[:,1]*1e6 .- Species_neg[:,1]*1e6)./2, color=:navy,fillalpha=0.35,linecolor=:navy,label="")
 
     




   savefig(fig1,"mrna.pdf")

