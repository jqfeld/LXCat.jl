using LXCat, Dates, Test

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
@test cs.type.threshold_energy ≈ 1.15e1 
@test cs.comment == "All excitation is grouped into this one level."
@test cs.updated == DateTime("2010-06-23T11:41:34")
@test cs.cross_section(1e3) ≈ 1.77e-21
