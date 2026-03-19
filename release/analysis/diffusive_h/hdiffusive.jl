using Pkg; Pkg.activate(".")

using NCDatasets
using CSV, DataFrames
using Dates
using Statistics
using PyPlot

ds = NCDataset("data/EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc")
# ds[:ustar]
psltime = DateTime(2020,1,1) .+ Second.(round.(ds[:time][:]))

tb = CSV.read("table2.txt", DataFrame, delim=' ', header=false, ignorerepeated=true)

kappa = 0.4
cptime = tb[!,1]
lambdai = 1e-3 .* tb[!,8] # 1/s

# psl indices of cold pools
idx = [searchsortedlast(psltime, cpt) for cpt in cptime]
ustar = [mean(ds[:ustar][i:i+5]) for i in idx]

hD = kappa * ustar ./ lambdai
hi = tb[!,end]

plot(hi, hD, marker="o", linestyle="none")
plot([0,150],[0,150],"k-")
plot([0,150],[0,5*150],"k-")
xlabel(L"h_i"); ylabel(L"h_D")

