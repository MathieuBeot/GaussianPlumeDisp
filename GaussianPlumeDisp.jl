module GaussianPlumeDisp

using Plots
using DataFrames
using CSV

data = CSV.read("Param_disp.csv", DataFrame)
data_rural = data[1:end, 1:13]
data_urban = hcat(data[1:end, 1], data[1:end, 14:end])
rename!(data_urban, 1 => "x (m)")

p1 = plot(data_rural[!, 1], data_rural[!, 2], label="A", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 3], label="B", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 4], label="C", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 5], label="D", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 6], label="E", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 7], label="F", leg=:topleft, xscale=:log10, yscale=:log10)
xlabel!("x (m)")
ylabel!("sigma_y")
title!("Rural")

p2 = plot(data_rural[!, 1], data_rural[!, 8], label="A", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 9], label="B", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 10], label="C", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 11], label="D", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 12], label="E", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_rural[!, 1], data_rural[!, 13], label="F", leg=:topleft, xscale=:log10, yscale=:log10)
xlabel!("x (m)")
ylabel!("sigma_z")
title!("Rural")

p3 = plot(data_urban[!, 1], data_urban[!, 2], label="A-B", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_urban[!, 1], data_urban[!, 3], label="C", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_urban[!, 1], data_urban[!, 4], label="D", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_urban[!, 1], data_urban[!, 5], label="E-F", leg=:topleft, xscale=:log10, yscale=:log10)
xlabel!("x (m)")
ylabel!("sigma_y")
title!("Urban")

p4 = plot(data_urban[!, 1], data_urban[!, 6], label="A-B", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_urban[!, 1], data_urban[!, 7], label="C", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_urban[!, 1], data_urban[!, 8], label="D", leg=:topleft, xscale=:log10, yscale=:log10)
plot!(data_urban[!, 1], data_urban[!, 9], label="E-F", leg=:topleft, xscale=:log10, yscale=:log10)
xlabel!("x (m)")
ylabel!("sigma_z")
title!("Urban")

DispersionParam() = plot(p1, p3, p2, p4, layout(2, 2))

Concentration(Q, u, sigma_y, sigma_z, y, z, h) = (Q/(2*pi*u*sigma_y*sigma_z))*exp(-0.5*y^2/sigma_y^2)*(exp(-0.5*(z-h)^2/sigma_z^2)+exp(-0.5*(z+h)^2/sigma_z^2))

concA = Concentration.(5, 2, data_rural[1:11, 2], data_rural[1:11, 8], 0, 0, 30)*10^6
concB = Concentration.(5, 2, data_rural[1:11, 3], data_rural[1:11, 9], 0, 0, 30)*10^6
concC = Concentration.(5, 2, data_rural[1:11, 4], data_rural[1:11, 10], 0, 0, 30)*10^6
concD = Concentration.(5, 2, data_rural[1:11, 5], data_rural[1:11, 11], 0, 0, 30)*10^6
concE = Concentration.(5, 2, data_rural[1:11, 6], data_rural[1:11, 12], 0, 0, 30)*10^6
concF = Concentration.(5, 2, data_rural[1:11, 7], data_rural[1:11, 13], 0, 0, 30)*10^6

p5 = plot(data_rural[1:11, 1], concA, label="A", leg=:topright)
plot!(data_rural[1:11, 1], concB, label="B", leg=:topright)
plot!(data_rural[1:11, 1], concC, label="C", leg=:topright)
plot!(data_rural[1:11, 1], concD, label="D", leg=:topright)
plot!(data_rural[1:11, 1], concE, label="E", leg=:topright)
plot!(data_rural[1:11, 1], concF, label="F", leg=:topright)
xlabel!("x (m)")
ylabel!("c (µg.m^-3)")
title!("Rural")

conc_AB = Concentration.(5, 2, data_urban[1:9, 2], data_urban[1:9, 6], 0, 0, 30)*10^6
conc_C = Concentration.(5, 2, data_urban[1:9, 3], data_urban[1:9, 7], 0, 0, 30)*10^6
conc_D = Concentration.(5, 2, data_urban[1:9, 4], data_urban[1:9, 8], 0, 0, 30)*10^6
conc_EF = Concentration.(5, 2, data_urban[1:9, 5], data_urban[1:9, 9], 0, 0, 30)*10^6

p6 = plot(data_urban[1:9, 1], conc_AB, label="A-B", leg=:topright)
plot!(data_urban[1:9, 1], conc_C, label="C", leg=:topright)
plot!(data_urban[1:9, 1], conc_D, label="D", leg=:topright)
plot!(data_urban[1:9, 1], conc_EF, label="E-F", leg=:topright)
xlabel!("x (m)")
ylabel!("c (µg.m^-3)")
title!("Urban")

PlumeConcentration() = plot(p5, p6, layout = (2, 1))

end
