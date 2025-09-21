using LXCat, Interpolations, Dates, Test

string = """EXCITATION
Ar -> Ar*(11.5eV)
  1.150e+1   / threshold energy
COMMENT: All excitation is grouped into this one level.
UPDATED: 2010-06-23 11:41:34
------------------------------------------------------------
 1.0000e+3	1.7700e-21
 1.5000e+3	1.3600e-21
 2.0000e+3	1.1000e-21
 3.0000e+3	8.3000e-22
 5.0000e+3	5.8000e-22
 7.0000e+3	4.5000e-22
 1.0000e+4	3.5000e-22
------------------------------------------------------------
"""

cs = parse_string(string)

@test cs.type.projectile == "e"
@test cs.type.target == "Ar"
@test cs.type.excited_state == "Ar*(11.5eV)"
@test cs.type.threshold_energy == 1.15e1
@test cs.comment == "All excitation is grouped into this one level."
@test cs.updated == DateTime("2010-06-23T11:41:34")
@test cs.cross_section(1e3) == 1.77e-21
@test collect(Interpolations.knots(cs.cross_section)) == [1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 7000.0, 10000.0]

bidirectional_string = """EXCITATION
N2(v=0) <-> N2(v=1)
  2.500000e+0  3.500000e-1 / threshold energy and statistical weight ratio
UPDATED: 2021-01-01 12:34:56
------------------------------------------------------------
 0.0000e+0      0.0000e+0
 2.5000e+0      2.5000e-21
 5.0000e+0      5.0000e-21
------------------------------------------------------------
"""

bidirectional_cs = parse_string(bidirectional_string)

@test bidirectional_cs.type.projectile == "e"
@test bidirectional_cs.type.target == "N2(v=0)"
@test bidirectional_cs.type.excited_state == "N2(v=1)"
@test bidirectional_cs.type.threshold_energy == 2.5
@test bidirectional_cs.type.stat_weight_ratio == 3.5e-1
@test bidirectional_cs.comment == ""
@test bidirectional_cs.updated == DateTime("2021-01-01T12:34:56")
@test bidirectional_cs.cross_section(3.75) ≈ 3.75e-21 atol=eps()
@test collect(Interpolations.knots(bidirectional_cs.cross_section)) == [0.0, 2.5, 5.0]
