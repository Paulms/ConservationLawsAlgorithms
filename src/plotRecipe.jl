@recipe function f(sol::AbstractFVSolution, t = size(sol.t,1))
    seriestype  :=  :path
    xguide --> "x"
    yguide --> "u"
    labels = String[]
    for i in 1:size(sol.u,2)
      push!(labels,"u$i")
    end
    label --> reshape(labels,1,length(labels))
    sol.prob.mesh.x, sol.u[t]
end
