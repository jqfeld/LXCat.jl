using LXCat, Interpolations, Dates, Test

db = load_database("../data/SigloDataBase-LXCat-04Jun2013.txt")


@test haskey(db, ("Ar", "Ar^+"))
@test db["Ar", "effective"](1e1) == 1.5e-19

@test haskey(db, ("N2", "N2^+"))

