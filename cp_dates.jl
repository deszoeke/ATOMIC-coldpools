using Dates
using NCDatasets
using MAT
using PyPlot
using PyCall

"convert matlab datenumber to DateTime"
mat2dt(mdn) = DateTime(0,1,1) - Day(1) + Second(round(Int64,mdn*86_400))

f = matread("./data/conserved_variables_10minLIMITED.mat")
g = matread("data/mixing_fractions_vars.mat")

ii = [2:11; 13:16]
map(mat2dt, [g["tcold"][ii,1]  g["tcoldw"][ii,1] g["mf_t_cp"][ii,1]])

