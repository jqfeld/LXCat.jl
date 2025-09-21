using LXCat, Interpolations, Dates, Test

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
@test cs.type.mass_ratio == 3.424e-5
@test cs.comment == ""
@test cs.updated == DateTime("2022-09-02T11:05:07")
@test cs.cross_section(1e3) == 3.14e-21
@test collect(Interpolations.knots(cs.cross_section)) == [0.0, 1.0, 1000.0]

unsorted_string = """ELASTIC
He
  2.500000e-5 / mass ratio
UPDATED: 2020-01-01 00:00:00
------------------------------------------------------------
 2.0000e+0      2.0000e-20
 0.0000e+0      0.0000e+0
 1.0000e+0      1.0000e-20
------------------------------------------------------------
"""

unsorted_cs = parse_string(unsorted_string)

@test unsorted_cs.type.projectile == "e"
@test unsorted_cs.type.target == "He"
@test unsorted_cs.type.mass_ratio == 2.5e-5
@test unsorted_cs.updated == DateTime("2020-01-01T00:00:00")
@test unsorted_cs.cross_section(0.5) ≈ 5.0e-21 atol=eps()
@test unsorted_cs.cross_section(1.5) ≈ 1.5e-20 atol=eps()
@test collect(Interpolations.knots(unsorted_cs.cross_section)) == [0.0, 1.0, 2.0]
