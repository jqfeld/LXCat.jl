using LXCat, Interpolations, Dates, Test

db = load_database("../data/SigloDataBase-LXCat-04Jun2013.txt")


@test first(filter(x -> x isa CrossSection{Effective} && x.type.target == "Ar",db))(1e1) == 1.5e-19

