using LXCatHelpers, Interpolations, Dates, Test

string = """EXCITATION
Ar -> Ar*(11.5eV)
  1.150e+1   / threshold energy
COMMENT: All excitation is grouped into this one level.
UPDATED: 2010-06-23 11:41:34
------------------------------------------------------------
 1.1500e+1	 0.0000e+0
 1.2700e+1	7.0000e-22
 1.3700e+1	1.4100e-21
 1.4700e+1	2.2800e-21
 1.5900e+1	3.8000e-21
 1.6500e+1	4.8000e-21
 1.7500e+1	6.1000e-21
 1.8500e+1	7.5000e-21
 1.9900e+1	9.2000e-21
 2.2200e+1	1.1700e-20
 2.4700e+1	1.3300e-20
 2.7000e+1	1.4200e-20
 3.0000e+1	1.4400e-20
 3.3000e+1	1.4100e-20
 3.5300e+1	1.3400e-20
 4.2000e+1	1.2500e-20
 4.8000e+1	1.1600e-20
 5.2000e+1	1.1100e-20
 7.0000e+1	9.4000e-21
 1.0000e+2	7.6000e-21
 1.5000e+2	6.0000e-21
 2.0000e+2	5.0500e-21
 3.0000e+2	3.9500e-21
 5.0000e+2	2.8000e-21
 7.0000e+2	2.2500e-21
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

@test cs.ground_state == "Ar" 
@test cs.excited_state == "Ar*(11.5eV)" 
@test cs.threshold_energy == 1.15e1 
@test cs.comment == "All excitation is grouped into this one level."
@test cs.updated == DateTime("2010-06-23T11:41:34")
@test cs.cross_section(1e3) == 1.77e-21
