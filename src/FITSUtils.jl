VERSION < v"0.1.0" && __precompile__()

module FITSUtils
export get_axis
export get_name_list
export get_spectra

using FITSIO
using Printf

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
"""
function get_spectra(hdu::ImageHDU)
  header = read_header(hdu)
  data = read(hdu)

  nlist = get_name_list(header, 4)

  rsun = 8.3
  xaxis = get_axis(header, 1)
  ilow = findlast(x->x<rsun, xaxis)
  iup = ilow + 1
  wlow = (xaxis[iup] - rsun) / (xaxis[iup] - xaxis[ilow])
  wup = (rsun - xaxis[ilow]) / (xaxis[iup] - xaxis[ilow])

  spectra = map((flow, fup)->flow * wlow + fup * wup, data[ilow,1,:,:], data[iup,1,:,:])
  result = Dict{String,Array{Float64,1}}()

  eaxis = map(x->10^x / 1e3, get_axis(header, 3)) # [GeV]
  for i in 1:length(nlist)
    result[nlist[i]] = map((e,f)->f / e^2 / 1e3, eaxis, spectra[:,i]) # MeV^2 cm^-2 sr^-1 s^-1 MeV^-1 -> cm^-2 sr^-1 s^-1 GeV^-1
  end
  result["eaxis"] = eaxis

  return result
end

end  # FITSUtils
