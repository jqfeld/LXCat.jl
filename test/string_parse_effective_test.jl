using LXCat, Dates, Test

string = """EFFECTIVE
Ar
  1.360e-5 / mass ratio
COMMENT: EFFECTIVE MOMENTUM-TRANSFER CROSS SECTION
UPDATED: 2011-06-06 18:21:14
------------------------------------------------------------
 7.0000e-1	8.6000e-21
 1.0000e+4	1.7500e-21
------------------------------------------------------------
"""

cs = parse_string(string)


@test cs.type.projectile == "e" 
@test cs.type.target == "Ar" 
@test cs.type.mass_ratio == 1.360e-5
@test cs.comment == "EFFECTIVE MOMENTUM-TRANSFER CROSS SECTION"
@test cs.updated == DateTime("2011-06-06T18:21:14")
@test cs.cross_section(1e4) ≈ 1.75e-21
