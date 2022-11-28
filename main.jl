using Revise
using LXCat
using Plots
using Interpolations

path = "e+He__e+e+He(+).txt"

cross_sections = load_database(path)

bounds = Interpolations.bounds(cross_sections[1].cross_section.itp)
start = 0.0
stop = bounds[1][2]
steps = 0.1
length = 100_000
x_vals = range(start=start, stop=stop, length=length)
y_vals = cross_sections[1].cross_section(x_vals)
plot(x_vals, y_vals)

cross_sections_log = load_database(path; interpolation_type=:logarithm)

bounds = Interpolations.bounds(cross_sections_log[1].cross_section.itp)
# start = bounds[1][1]
# stop = bounds[1][2] 
# length = size(cross_sections_log[1].cross_section, 1)
x_vals = range(start=start, stop=stop, length=length)
interpolated = cross_sections_log[1].cross_section(x_vals)
y_vals = 10 .^ (interpolated / 10)
@. y_vals[(isnan(y_vals))] = 0
y_vals
plot!(x_vals, y_vals)


