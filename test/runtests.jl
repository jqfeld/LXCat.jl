using LXCatHelpers, Test, SafeTestsets

@safetestset "String Parsing Effective Cross Section" begin include("string_parse_effective_test.jl") end
@safetestset "String Parsing Elastic Cross Section" begin include("string_parse_elastic_test.jl") end
@safetestset "String Parsing Excitation Cross Section" begin include("string_parse_excitation_test.jl") end
@safetestset "String Parsing Ionization Cross Section" begin include("string_parse_ionization_test.jl") end

@safetestset "Database Parsing" begin include("load_database_test.jl") end
