using Test
using Polyhedra
using CDDLib

function test_incident_functions(p::CDDLib.Polyhedron)
    T = Polyhedra.coefficient_type(p)
    tol = Polyhedra._default_tol(T)

    for (hsidx, hs) in zip(eachindex(halfspaces(p)), halfspaces(p))
        inc_pts = @inferred incidentpoints(p, hsidx)
        inc_ptinds = @inferred incidentpointindices(p, hsidx)
        for pts in (inc_pts, get.(p, inc_ptinds))
            for pt in pts
                @test Polyhedra.isincident(pt, hs, tol=tol)
            end
        end
        # Test against the generic versions in Polyhedra.jl
        @test inc_pts == incidentpoints(p, hsidx, tol=tol)
        @test inc_ptinds == incidentpointindices(p, hsidx, tol=tol)
    end

    for (ptidx, pt) in zip(eachindex(points(p)), points(p))
        inc_hss = @inferred incidenthalfspaces(p, ptidx)
        inc_hsinds = @inferred incidenthalfspaceindices(p, ptidx)
        for hss in (inc_hss, get.(p, inc_hsinds))
            for hs in hss
                @test Polyhedra.isincident(pt, hs, tol=tol)
            end
        end
        # Test against the generic versions in Polyhedra.jl
        @test inc_hss == incidenthalfspaces(p, ptidx, tol=tol)
        @test inc_hsinds == incidenthalfspaceindices(p, ptidx, tol=tol)
    end
end

@testset "Incidence tests" begin
    @testset "Incidence $precision" for precision in [:float, :exact]
        A = [1 1; 1 -1; -1 0]; b = [1, 0, 0]
        incidence_extected = Set([
            BitSet([1, 2]),
            BitSet([1, 3]),
            BitSet([2, 3])
        ])
        p = polyhedron(hrep(A, b), CDDLib.Library(precision))
        @test p isa CDDLib.Polyhedron{precision == :float ? Float64 : Rational{BigInt}}
        vrep(p)
        @test Set(copyincidence(p.poly)) == incidence_extected
    end

    @testset "Input incidence $precision" for precision in [:float, :exact]
        V = [[1//2, 1//2], [0, 1], [0, 0]]
        p = polyhedron(vrep(V), CDDLib.Library(precision))
        hrep(p)
        hs = collect(halfspaces(p))

        incidence_computed = copyinputincidence(p.poly)
        for (vidx, v) in enumerate(points(p))
            for i in incidence_computed[vidx]
                @test Polyhedra.isincident(v, hs[i], tol=0)
            end
        end
    end

    @testset "get[hv]incidence $precision" for precision in [:float, :exact]
        A = [1 1; 1 -1; -1 0]; b = [1, 0, 0]
        p_H = polyhedron(hrep(A, b), CDDLib.Library(precision))
        vrep(p_H)

        V = [[1//2, 1//2], [0, 1], [0, 0]]
        p_V = polyhedron(vrep(V), CDDLib.Library(precision))
        hrep(p_V)

        # Homogeneous cone
        A = [-1 0; 0 -1]; b0 = [0, 0]
        p_hc = polyhedron(hrep(A, b0), CDDLib.Library(precision))

        # Non-homogeneous cone
        A = [-1 0; 0 -1]; b1 = [-1, -1]
        p_nhc = polyhedron(hrep(A, b1), CDDLib.Library(precision))

        for p in [p_H, p_V, p_hc, p_nhc]
            @inferred CDDLib.gethincidence(p)
            @inferred CDDLib.getvincidence(p)

            hs = collect(halfspaces(p))
            vs = [collect(rays(p))..., collect(points(p))...]

            T = Polyhedra.coefficient_type(p)
            tol = Polyhedra._default_tol(T)

            for (vidx, v) in enumerate(vs)
                for i in p.hincidence[vidx]
                    @test Polyhedra.isincident(v, hs[i], tol=tol)
                end
            end

            for (hidx, h) in enumerate(hs)
                for i in p.vincidence[hidx]
                    @test Polyhedra.isincident(vs[i], h, tol=tol)
                end
            end
        end
    end

    @testset "incident* functions $precision" for precision in [:float, :exact]
        @testset "Points + lines" begin
            A = [-1 0; 1 0]; b = [0, 1]
            p = polyhedron(hrep(A, b), CDDLib.Library(precision))
            T = Polyhedra.coefficient_type(p)
            tol = Polyhedra._default_tol(T)

            for (hsidx, hs) in zip(eachindex(halfspaces(p)), halfspaces(p))
                inc_pts = @inferred incidentpoints(p, hsidx)
                inc_ptinds = @inferred incidentpointindices(p, hsidx)
                for pts in (inc_pts, get.(p, inc_ptinds))
                    for pt in pts
                        @test Polyhedra.isincident(pt, hs, tol=tol)
                    end
                end
                # Test against the generic versions in Polyhedra.jl
                @test inc_pts == incidentpoints(p, hsidx, tol=tol)
                @test inc_ptinds == incidentpointindices(p, hsidx, tol=tol)
            end

            for (ptidx, pt) in zip(eachindex(points(p)), points(p))
                inc_hss = @inferred incidenthalfspaces(p, ptidx)
                inc_hsinds = @inferred incidenthalfspaceindices(p, ptidx)
                for hss in (inc_hss, get.(p, inc_hsinds))
                    for hs in hss
                        @test Polyhedra.isincident(pt, hs, tol=tol)
                    end
                end
                # Test against the generic versions in Polyhedra.jl
                @test inc_hss == incidenthalfspaces(p, ptidx, tol=tol)
                @test inc_hsinds == incidenthalfspaceindices(p, ptidx, tol=tol)
            end
        end
    end
end
