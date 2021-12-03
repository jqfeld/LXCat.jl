using LXCat, Interpolations, Dates, Test

string = """ELASTIC
CH4
  3.424e-5 / mass ratio
UPDATED: 2011-02-08 11:05:07
------------------------------------------------------------
 0.0000e+0	4.0000e-19
 1.0000e+3	3.2000e-21
------------------------------------------------------------
"""

cs = parse_string(string)

@test cs.ground_state == "CH4" 
@test cs.mass_ratio == 3.424e-5
@test cs.comment == ""
@test cs.updated == DateTime("2011-02-08T11:05:07")
@test cs.cross_section(1e3) == 3.2e-21
