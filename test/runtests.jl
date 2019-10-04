using Test

using FITSUtils
using FITSIO

import FITSUtils:
get_axis

@testset "FITSUtils" begin
  hdu = FITS("test.fits")[1]
  header = read_header(hdu)

  @testset "get_axis" begin
    axis = get_axis(header, 1)
    @test axis[1] == header["CRVAL1"] && (axis[2] - axis[1] - header["CDELT1"]) / header["CDELT1"] < 1e-8
  end
end
