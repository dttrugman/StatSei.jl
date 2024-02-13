###### Map Projection Tools for Seismicity
#   Daniel Trugman, 2024

### Function to setup map projection based on input parameters
#   - Adapted from GrowClust3D.jl by DTT
function setup_projection(lons,lats,mapproj="tmerc",rellipse="WGS84")

    # extract center
    plon0, plat0 = median(lons), median(lats)

    ## setup projection
    if mapproj in ["aeqd", "tmerc"]
        proj = Transformation("+proj=longlat +ellps=$rellipse",
            "+proj=$mapproj +ellps=$rellipse +lat_0=$plat0 +lon_0=$plon0 +units=km")
    elseif mapproj == "merc"
        proj = Transformation("+proj=longlat +ellps=$rellipse",
            "+proj=$mapproj +ellps=$rellipse +lat_ts=$plat0 +lon_0=$plon0 +units=km")
    elseif mapproj == "lcc"
        plat1, plat2 = quantile(lats,[0.10,0.90])
        proj = Transformation("+proj=longlat +ellps=$rellipse",
        "+proj=$mapproj +ellps=$rellipse +lat_0=$plat0 +lon_0=$plon0 +lat_1= $plat1 +lat_2=$plat2 +units=km")
    else
        println("ERROR, map projection not defined! ", mapproj)
        exit()
    end

    # return projection object
    return proj

end

### Forward transform using transform object fproj ###
#   - Adapted from GrowClust3D.jl by DTT
function lonlat2xypos(lons::Vector{Float64},lats::Vector{Float64},fproj) # forward projection object

    # initialize
    n = length(lons)
    xx, yy = zeros(n), zeros(n)

    # project each lon lat point
    for ii in eachindex(xx)
        xx[ii], yy[ii] = fproj(lons[ii],lats[ii])
    end

    # return
    return xx, yy

end
function lonlat2xypos(lons::Vector{Float32},lats::Vector{Float32},fproj) # forward projection object

    # initialize
    n = length(lons)
    xx, yy = zeros(Float32,n), zeros(Float32,n)

    # project each lon lat point
    for ii in eachindex(xx)
        xx[ii], yy[ii] = fproj(lons[ii],lats[ii])
    end

    # return
    return xx, yy

end

### Inverse transform using transform object iproj
#   - Adapted from GrowClust3D.jl by DTT
function xypos2latlon(xx::Vector{Float64},yy::Vector{Float64},iproj) # inverse projection object

    # initialize
    n = length(xx)
    lons, lats = zeros(n), zeros(n)

    # inverse project each lon lat point
    for ii in eachindex(xx)
        lons[ii], lats[ii] = iproj(xx[ii], yy[ii])
    end

    # return
    return lons, lats

end
function xypos2latlon(xx::Vector{Float32},yy::Vector{Float32},iproj) # inverse projection object

    # initialize
    n = length(xx)
    lons, lats = zeros(Float32,n), zeros(Float32,n)

    # inverse project each lon lat point
    for ii in eachindex(xx)
        lons[ii], lats[ii] = iproj(xx[ii], yy[ii])
    end

    # return
    return lons, lats

end