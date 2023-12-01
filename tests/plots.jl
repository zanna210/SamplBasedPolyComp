using GLMakie, DelimitedFiles, Polynomials

dms = [10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,1000]
colors = ["red", "blue", "green"]
scatter_colors = ["orange", "cyan", "light green"]
funs = ["ConvexHull"]
f = Figure(resolution=(1200,800))
Axis(f[1,1], title="ConvexHull time", xlabel="Dimensionality of the instances", ylabel="elapsed time (s)")

for j in 1:size(funs)[1]
        
    #n_samples=[100,200,300,400,500]
    #brb = "C:/Users/david/University/TESI/tests/time results/$(fun)/$(fun)_200x$(dm).txt"
    #times = readdlm(brb)
    
    avgs = zeros(size(dms)[1])
    i = 0
    for k in dms
        i +=1
        T = vec(readdlm("C:/Users/david/University/TESI/tests/time results/$(funs[j])/$(funs[j])_200x$k.txt"))
        avgs[i] = sum(T) / size(T)[1]
        scatter!(f[1,1], [k for _ in 1:size(T)[1]], T, color=scatter_colors[j], markersize=7, marker=:xcross)
    end
    
    scatter!(f[1,1], dms, avgs, color=colors[j], markersize=10)
    lines!(f[1,1], dms, avgs, color=colors[j], label="$(funs[j]) average elapsed time")

    #=
    xs = 0:10:dms[size(dms)[1]]
    p = fit(dms, T ,1)
    #ys = coeffs(p)[size(coeffs(p))[1]] * xs.^2
    ys = p.(xs)
    #lines!(f[1,1], xs, ys, color="dark green", label="x^2 reference")
    =#
end
axislegend(position=:lt)

display(f)
save("C:/Users/david/University/TESI/tests/time results/ConvexHull_plot.png", f)