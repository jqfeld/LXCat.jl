using LXCat, Dates, Test

string = """ELASTIC
CH4
  3.424e-5 / mass ratio
UPDATED: 2022-09-02 11:05:07
------------------------------------------------------------
 0.0000e+0	1.1100e-19
 1.0    121.2e-18
 1.0000e+3	3.14000e-21
------------------------------------------------------------
"""

cs = parse_string(string)

@test cs.type.projectile == "e" 
@test cs.type.target == "CH4" 
@test cs.type.mass_ratio ≈ 3.424e-5
@test cs.comment == ""
@test cs.updated == DateTime("2022-09-02T11:05:07")
@test cs.cross_section(1e3) ≈ 3.14e-21
