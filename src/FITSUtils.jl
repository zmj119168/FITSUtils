VERSION < v"0.1.0" && __precompile__()

module FITSUtils
export get_axis
export get_galprop_axis
export get_name_list
export get_spectra
export get_scatter
export count_R
export count_Ekin
export Particle
export den_Particle
export get_density

using Interpolations
using FITSIO
using Printf

mutable struct Particle
  dNdE::Array{T,1} where {T<:Real}
  R::Array{T,1} where {T<:Real}
  dNdR::Array{T,1} where {T<:Real}
  Ekin::Array{T,1} where {T<:Real}
  A::Int
  Z::Int
end

mutable struct den_Particle
  density::Array{T,1} where {T<:Real}
  radius::Array{T,1} where {T<:Real}
end

function Base.:copy(p::Particle)
  Particle(copy(p.dNdE), copy(p.R), copy(p.dNdR), copy(p.Ekin), p.A, p.Z)
end

"""
    get_axis(header::FITSHeader, index::Int)

Get the grid array of a specified axis from an hdu header

# Arguments
* `header`: the header.
* `index`:  the index of specfied axis.
"""
function get_axis(header::FITSHeader, ind::Int)
  return map(i->header["CRVAL$ind"] + header["CDELT$ind"] * (i - 1), 1:header["NAXIS$ind"])
end

"""
    get_galprop_axis(header::FITSHeader)

Get the grid array of a axises from an GALPROP result hdu header

# Arguments
* `header`: the header.
"""
function get_galprop_axis(header::FITSHeader)
  result = Dict{String,Array{Real,1}}()
  pairs = header["NAXIS"] == 5 ? zip(["x", "y", "z", "E"], 1:4) : zip(["x", "E"], [1, 3])

  for pair in pairs
    result[pair[1]] = get_axis(header, pair[2])
  end
  result["E"] = map(e->10^e / 1e3, result["E"])
  return result
end

"""
    get_name_list(header::FITSHeader, index::Int)

Get the list of the exist particle (for GALPROP output)

# Arguments
* `header`: the header.
* `index`:  the index of the particle axis.
"""
function get_name_list(header::FITSHeader, index::Int)
  return map(i->get_name(header, i), 1:header["NAXIS$index"])
end
function get_name(header, ind::Int)
  index = @sprintf "%03d" ind
  return header["NAME$index"]
end

"""
    get_spectra(hdu::ImageHDU)

    Get the spectra of the exist particles (for GALPROP output)

# Arguments
* `hdu`: the hdu for GALPROP output.
* `r`: the place of spectra,default to sun=8.3.
"""
function get_spectra(hdu::ImageHDU;r::Real = 8.3)
  header = read_header(hdu)
  data = read(hdu)
  nlist = get_name_list(header, 4)

  rsun = r
  xaxis = get_axis(header, 1)
  ilow = findlast(x->x<rsun, xaxis)
  iup = ilow + 1
  wlow = (xaxis[iup] - rsun) / (xaxis[iup] - xaxis[ilow])
  wup = (rsun - xaxis[ilow]) / (xaxis[iup] - xaxis[ilow])

  spectra = map((flow, fup)->flow * wlow + fup * wup, data[ilow,1,:,:], data[iup,1,:,:])
  result = Dict{String,Particle}()

  eaxis = map(x->10^x / 1e3, get_axis(header, 3)) # [GeV]
  for i in 1:length(nlist)
    index = @sprintf "%03d" i
    iZ=header["NUCZ$index"]
    iA=header["NUCA$index"]
    idnde = map((e,f)->f / e^2 / 1e3, eaxis, spectra[:,i]) # MeV^2 cm^-2 sr^-1 s^-1 MeV^-1 -> cm^-2 sr^-1 s^-1 GeV^-1
    tmp = Particle(idnde, Array{Real,1}(), Array{Real,1}(), eaxis, iA, iZ)
    count_R(tmp)
    result[nlist[i]] = tmp
  end

  result
end


function get_density(hdu::ImageHDU;low::Real = 10)
  c=2.99792458e10
  header = read_header(hdu)
  data = read(hdu)
  nlist = get_name_list(header, 4)
  xaxis = get_axis(header, 1)
  result = Dict{String,den_Particle}()
  eaxis = map(x->10^x / 1e3, get_axis(header, 3)) # [GeV]
  int=log(low):1e-4:7
  for n in 1:length(nlist)
    index = @sprintf "%03d" n
    iZ=abs(header["NUCZ$index"])
    iA=header["NUCA$index"]
    m0 =iA == 0 ? 0.511e-3 : 0.9382   
    de= Array{Real,1}(undef, 0)
    for r in 1:length(xaxis)
      flux = @. data[r,1,:,n] / eaxis^2 / 1e3               # MeV^2 cm^-2 sr^-1 s^-1 MeV^-1 -> cm^-2 sr^-1 s^-1 GeV^-1     
      logene = log.(eaxis)
      spec = extrapolate(interpolate((range(logene[1],last(logene),length=length(logene)),), log.(flux), Gridded(Linear())), Line())
      fun12(e)=4e12pi*exp(spec(log(e)))/(c*sqrt(1-(m0/(m0+e))^2))
      de=vcat(de,sum([(fun12(exp(i))+fun12(exp(i+di)))*(exp(i+di)-exp(i))/2 for i in int]))   
    end
    result[nlist[n]] = den_Particle(de, xaxis)
  end
  result
end

function count_R(particle::Particle)
  m0 = particle.A == 0 ? 0.511e-3 : 0.9382
  E = @. particle.Ekin + m0
  p = @. sqrt(E^2 - m0^2)
  particle.R = p * particle.A / abs(particle.Z) # [GV]
  particle.dNdR = @. particle.dNdE * p/E * abs(particle.Z)/particle.A # cm^-2 sr^-1 s^-1 GV^-1
  particle
end

function count_Ekin(particle::Particle)
  m0 = particle.A == 0 ? 0.511e-3 : 0.9382
  p = @. particle.R * abs(particle.Z) / particle.A # [GeV]
  E = @. sqrt(p^2 + m0^2)
  particle.Ekin = @. E - m0 # [GV]
  particle.dNdE = @. particle.dNdR * E/p * particle.A/abs(particle.Z) # cm^-2 sr^-1 s^-1 GeV^-1
  particle
end

function get_scatter(density::Array{T,4} where {T<:Real}, axis::Dict{String,Array{T,1}} where {T<:Real}, p::Real)
  iup = findfirst(v->v>p, axis["E"])
  iup = iup == nothing ? length(axis["E"]) : iup
  iup = iup == 1 ? 2 : iup

  ilow = iup - 1
  plow, p, pup = log(axis["E"][ilow]), log(p), log(axis["E"][iup])
  weights = [[ilow, (pup-p)/(pup-plow)], [iup, (p-plow)/(pup-plow)]]

  cube = mapreduce(v->density[:,:,:,floor(Int, v[1])]*v[2], +, weights, init=zeros(typeof(density[1,1,1,1]),size(density)[1:3]))

  scatters = Array{Array{Real,1},1}([[],[],[],[]])
  resize!(scatters, 4)

  for i in CartesianIndices(cube)
    (cube[i] == 0) && continue
    push!(scatters[1], axis["x"][i[1]])
    push!(scatters[2], axis["y"][i[2]])
    push!(scatters[3], axis["z"][i[3]])
    push!(scatters[4], cube[i])
  end

  return scatters
end
end  # FITSUtils
