using LXCat, Interpolations, Dates, Test

db = load_database("test_data.txt")


@test first(filter(x -> x isa CrossSection{Effective} && x.type.target == "Ar",db))(2e4) == 2.9500e-21
@test first(filter(x -> x isa CrossSection{Isotropic} && x.type.projectile == "Ar^+",db))(8.8e-2) == 4.516950e-19


