_incel(p::Polyhedron{T}, inc::Polyhedra.IncidentElements{T}, idx) where {T} = get(p, idx)
_incel(p::Polyhedron{T}, inc::Polyhedra.IncidentIndices{T}, idx) where {T} = idx

function Base.get(p::CDDLib.Polyhedron{T}, inc::Polyhedra.Incident{T, ElemT}; tol) where {T, ElemT<:HRepElement{T}}
    ine = CDDLib.getine(p)
    incidence = CDDLib.getincidence(p)
    incT = Polyhedra._inctype(inc)
    incs = incT[]
    for i in incidence[inc.idx.value]
        idx = Polyhedra.Index{T, ElemT}(i)
        isvalid(ine, idx) && push!(incs, _incel(p, inc, idx))
    end
    return incs
end
