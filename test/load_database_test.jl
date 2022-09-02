using LXCat, Interpolations, Dates, Test

db = load_database("test_data.txt")


@test first(filter(x -> x isa CrossSection{Effective} && x.type.target == "Ar",db))(2e4) == 2.9500e-21


