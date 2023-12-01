using DelimitedFiles
print("\n\n\n")
for n in [60,70,80,90,150,250,350,450]
    for k in 1:10
        
        name = "C:/Users/david/University/TESI/test samples/V_Rep/V_200x$(n)_$k.csv"
        local V = round.(rand(200,n), digits=2)
        writedlm(name, V, ";")
        
        #=
        name = "C:/Users/david/University/TESI/test samples/H_Rep/Ab_200x$(n)_$k.csv"
        local Ab = round.(rand(200,n+1), digits=2)
        Ab[:,n+1] = round.(Ab[:,n+1]*10, digits=2)
        writedlm(name, Ab, ";")
        =#
        #=
        name = "C:/Users/david/University/TESI/test samples/M/M_$(n)x$(n)_$k.csv"
        local M = round.(rand(n,n), digits=2)
        writedlm(name, M, ";")
        =#
    end
end
#=
V = round.(10*rand(5,5), digits=2)
writedlm("test.csv", V, ";")
display(V)
=#