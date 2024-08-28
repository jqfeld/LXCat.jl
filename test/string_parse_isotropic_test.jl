using LXCat, Dates, Test

string = """
SPECIES: Ar^+ / Ar
PROCESS: Ar+ + Ar -> , Isotropic
PARAM.:  Mi = 39.948, Mi/M = 1, complete set
COMMENT: Phelps database retrieved on 02.09.2022. Some datapoints omitted.
COMMENT: Ar+ in Ar: Phelps, J. Appl. Phys. 76, 747 (1994).
COMMENT: 2E-19/(2.*col(A))^0.5/(1.+2.*col(A))+3E-19*2.*col(A)/(1.+2.*col(A)/3.)^2.3.
UPDATED: 2016-03-25 18:22:30
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 0.000000e+0	1.413940e-17
 1.000000e-4	1.413940e-17
 1.200000e-4	1.290690e-17
 1.400000e-4	1.194900e-17
 1.600000e-4	1.117690e-17
 1.900000e-4	1.025600e-17
1.200000e-3	4.073430e-18
 1.400000e-3	3.769930e-18
 1.600000e-3	3.525210e-18
 1.900000e-3	3.233280e-18
1.000000e-2	1.392390e-18
 1.200000e-2	1.267810e-18
 1.400000e-2	1.170900e-18
 1.600000e-2	1.092730e-18
 1.900000e-2	9.994930e-19
5.200000e-2	5.905990e-19
 6.200000e-2	5.391950e-19
 7.400000e-2	4.925980e-19
 8.800000e-2	4.516950e-19
3.100000e+0	1.524680e-19
 3.700000e+0	1.359720e-19
 4.400000e+0	1.200290e-19
 7.400000e+0	7.721610e-20
 1.600000e+1	3.482290e-20
 2.600000e+1	1.991740e-20
 4.400000e+1	1.054790e-20
 6.200000e+1	6.891660e-21
 8.800000e+1	4.434790e-21
 1.400000e+2	2.455490e-21
 1.600000e+2	2.069450e-21
 1.000000e+4	9.687080e-24
-----------------------------
"""

cs = parse_string(string)

@test cs.type.projectile == "Ar^+" 
@test cs.type.target == "Ar" 

#TODO: maybe add newlines to the comments?
@test cs.comment == "Phelps database retrieved on 02.09.2022. Some datapoints omitted.Ar+ in Ar: Phelps, J. Appl. Phys. 76, 747 (1994).2E-19/(2.*col(A))^0.5/(1.+2.*col(A))+3E-19*2.*col(A)/(1.+2.*col(A)/3.)^2.3."

@test cs.updated == DateTime("2016-03-25T18:22:30")
@test cs.cross_section(1.6e2) ≈ 2.069450e-21
