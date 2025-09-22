using LXCat, Dates, Test

db = load_database("test_data.txt")

@test first(filter(x -> x isa CrossSection{Effective} && x.type.target == "Ar", db))(2e4) ≈ 2.9500e-21
@test first(filter(x -> x isa CrossSection{Isotropic} && x.type.projectile == "Ar^+", db))(8.8e-2) ≈ 4.516950e-19
@test first(filter(x -> x isa CrossSection{Excitation} && x.type.target == "He", db))(1.984e+1) ≈ 3.497e-23
@test first(filter(x -> x isa CrossSection{Attachment} && x.type.target == "H2O", db))(8.190000e+0) ≈ 7.800000e-24

