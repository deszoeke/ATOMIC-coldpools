using PyPlot

d = CSV.read("table1.csv", DataFrame, header=0)

function spxy(x,y; ax=gca(), xla="x", yla="y")
        r = cor(x,y)

        ax.plot(x,y, marker="o", linestyle="none")
        ax.set_xlabel(xla)
        ax.set_ylabel(yla)
        return r, ax
end

clf()
ax = subplot(2,2,1)
r,_ = spxy(d[!,2], d[!,3], ax=ax, xla=L"T_0-T_{min}"*" (°C)", yla=L"\delta"*"D (10\$^{-3}\$")
xlm = ax.set_xlim(0, 4)
ylm = ax.set_ylim(0, 9)
xt = [(xlm)...]'*[0.75, 0.25]
yt = [(ylm)...]'*[0.25, 0.75]
ax.text(xt, yt, "R\$^2\$=$(@sprintf("%0.2f", r^2))" )

# label each dot with a shorthand date
for row in eachrow(d[[(2:6)..., 8],:]) # not 10
    text(row[2]+0.1, row[3], row[1][[1, 5,6, 8,9]])
end
for row in eachrow(d[[1,7],:])
    text(row[2]-0.8, row[3]-0.7, row[1][[1, 5,6, 8,9]])
end
let row = d[[1,7],:]
    text(row[2]-0.8, row[3]-0.7, row[1][[1, 5,6, 8,9]])
end

[s[[1, 5,6, 8,9]] for s in d[!,1]]
#ax = subplot(2,2,2)
#

