facts("Low-level board tests") do
    A1 = -eye(Int, 9) # x >= 0
    b1 = zeros(Int, 9)
    A2 = eye(Int, 9) # x <= 1
    b2 = ones(Int, 9)
    A3 = zeros(Int, 9, 9)
    b3 = 3 * ones(Int, 9)
    i = 1
    for a = 1:3
      for b = (a+1):3
        for c = 1:3
          for d = (c+1):3
            ac = a + (c-1) * 3
            ad = a + (d-1) * 3
            bc = b + (c-1) * 3
            bd = b + (d-1) * 3
            A3[i, ac] = 1
            A3[i, ad] = 1
            A3[i, bc] = 1
            A3[i, bd] = 1
            i += 1
          end
        end
      end
    end
    A = [A1; A2; A3]
    b = [b1; b2; b3]
    ine = Polyhedra.SimpleHRepresentation(A, b)
    poly = CDDPolyhedra(ine)
    ext  = Polyhedra.SimpleVRepresentation{9, Rational{Int}}(copygenerators(poly))
    target = ones(Int, 9) * (3 // 4)
    ok = false
    for v in vrep(ext)
      if v == target
        ok = true
      end
    end
    @fact ok --> true

    cutA = ones(Int, 1, 9)
    cutb = 6
    Acut = [cutA; A]
    bcut = [cutb; b]
    inecut = Polyhedra.SimpleHRepresentation(Acut, bcut)
    (isredundant, certificate) = redundant(inecut, 1)
    @fact isredundant --> false
    @fact Array{Rational{Int}}(certificate) --> target
    redundantrows(inecut)
    @fact IntSet() --> redundantrows(inecut)
    (issredundant, scertificate) = sredundant(inecut, 1)
    @fact issredundant --> false
    @fact Array{Rational{Int}}(scertificate) --> target
end
