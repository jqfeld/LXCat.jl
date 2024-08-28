using LXCat, Dates, Test

string = """
IONIZATION
Ar -> Ar^+
  1.580e+1 / threshold energy
COMMENT: RAPP-SCHRAM
UPDATED: 2010-03-02 16:19:07
------------------------------------------------------------
 7.0000e+2	1.1500e-20
 1.0000e+3	8.6000e-21
------------------------------------------------------------
"""

cs = parse_string(string)

@test cs.type.projectile == "e" 
@test cs.type.target == "Ar" 
@test cs.type.excited_state == "Ar^+" 
@test cs.type.threshold_energy ≈ 1.58e1 
@test cs.comment == "RAPP-SCHRAM"
@test cs.updated == DateTime("2010-03-02T16:19:07")
@test cs.cross_section(1e3) ≈ 8.6e-21
