## Example from Bradlyn, Science, Eq. (S66)
# SG199
g3 = S"y+1,z,x+1"       # {3₁₁₁⁻|1,0,1}  
g1 = S"x-1/2,-y+1/2,-z" # 2₁₀₀|-½,½,0
g2 = S"-x,y+1/2,-z-1/2" # {2₀₁₀|0,½,-½}


G3 = [0.0 0 1; 1 0 0; 0 1 0]
G1 = [-1.0 0 0; 0 -1 0; 0 0 1]
G2 = inv(G3)*G1*(G3)

lg = LittleGroup(199, KVec("0.5,0.5,0.5"), "P", [g1,g2,g3,])#g1, g2, g3]) #
lgir = LGIrrep(
    "P?",
    lg,
    complex.([G1,G2,G3,]),#, G2, G3]), # 
    [zeros(3) for _ in 1:3],
    REAL,
    false
    )

H = kdotp(lgir)

## -----
# SG220

g3 = S"z,x,y"
g4 = S"-x+1/2,z+1,-y+1"
m = S"y+1/2,x+1/2,z+1/2"
G3 = -1.0im*[0 1 0; 0 0 -1; 1 0 0]
G4 = [0 0 1; 0 1 0; -1 0 0]*cispi(1/4)
M = [1 0 0; 0 0 1; 0 1 0]*cispi(3/4)

lg = LittleGroup(220, KVec("0.5,0.5,0.5"), "P", [g3, g4, m])
lgir = LGIrrep(
    "P?",
    lg,
    complex.([G3, G4, M]),
    [zeros(3) for _ in 1:3],
    REAL,
    false
    )
kdotp(lgir)