#module Plotting
    #export poly_plot, proj_plot

    using GLMakie
    using Polyhedra

    #---------------------------------------------------------------------

    function poly_plot(f, p, empty, c="yellow", alfa=0.5)
        mp = Polyhedra.Mesh(p)
        if empty
            mesh(f[1,1], mp, color=c, alpha=alfa)
        else
            mesh!(f[1,1], mp, color=c, alpha=alfa)
        end
        return f
    end

    #---------------------------------------------------------------------

    function proj_plot(f, projs, empty, c="red")
        if empty
            scatter(f[1, 1], projs, color=c, markersize=7)
        else
            scatter!(f[1, 1], projs, color=c, markersize=7)
        end
        return f
    end

    print("\nDone loading Plotting module\n")
#end #module