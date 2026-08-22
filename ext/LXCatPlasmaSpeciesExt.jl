module LXCatPlasmaSpeciesExt

using LXCat
using PlasmaSpecies: Species, StringGas

function LXCat.target_species(c; labels=Dict{String,String}())
  resolve_species(LXCat.target_label(c), labels)
end

function LXCat.product_species(c; labels=Dict{String,String}())
  label = LXCat.product_label(c)
  label === nothing ? nothing : resolve_species(label, labels)
end

# LoKI-notation parse with a StringGas fallback for labels PlasmaSpecies
# cannot parse (equality still works, mass does not).
function resolve_species(label::AbstractString, labels)
  s = String(get(labels, label, label))
  try
    Species(s)
  catch
    Species(; gas=StringGas(s))
  end
end

end
