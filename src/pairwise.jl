#### Functions for Nearest Neighbors Algorithm of Zaliapin and Ben-Zion, 2013
#    Daniel Trugman, Nevada Seismological Laboratory (2023)

################################## LLT FUNCTIONS ################################

# ### Calculate Earth Radius at a Reference Latitude
# function calculate_erad(lat::Float64)

#     # geographic constants, WGS84
#     a = 6378.1370000000000 # major axis
#     b = 6356.7523141999999 # minor axis

#     # calculate
#     r = sqrt( ( (a^2 *cosd(lat))^2 + (b^2 *sind(lat))^2 ) / (
#         (a*cosd(lat))^2 + (b*sind(lat))^2 ) )

#     # return in km
#     return r

# end

### Finds nearest neighbors from llt data - f64, using pairwise
function find_neighbors_pllt(lons::Vector{Float64}, lats::Vector{Float64}, times::Vector{Float64}, 
    mags::Vector{Float64}, bval::Float64, dfrac::Float64)

    # setup magnitude scaling and other arrays
    scaLM = -bval*mags
    X = transpose(hcat(lons,lats)) # meters
    T = times/(86400.0) # days

    # add noise for numerical stability
    X .+= 0.00001*randn(size(X))

    # get earth radius
    akm = calculate_erad(median(lats))

    # allocate arrays - output as Float64
    npts = length(lons)
    LNNvec = zeros(Float64,npts)
    LNRvec = zeros(Float64,npts)
    LNTvec = zeros(Float64,npts)
    Pvec = zeros(Int64,npts)

    # Pairwise Distances in float 64    
    LDX = dfrac*log10.(pairwise(Haversine(akm),X)) # km

    # Pairwise Times in float 64
    LDT = log10.(pairwise(Euclidean(),T) .- log10(365.0)) # years

    # Calculate Nearest Neighbors
    @inbounds for jj = 2:npts
        @views idx = argmin(LDX[1:jj-1,jj].+LDT[1:jj-1,jj].+scaLM[1:jj-1]) # find parent
        Pvec[jj] = idx # parent index
        LNNvec[jj] = LDX[idx,jj]+LDT[idx,jj]+scaLM[idx] # nearest neighbor distance
        LNRvec[jj] = LDX[idx,jj]+0.5*scaLM[idx] # rescaled distance
        LNTvec[jj] = LDT[idx,jj]+0.5*scaLM[idx] # rescaled time
        LRkm[jj] = LDX[idx,jj]/dfrac # distance in km
        LTday[jj] = LDT[idx,jj] + log10(365.0) # time offset in days
    end

    # return as DataFrame
    return DataFrame("enum"=>1:npts,"pid"=>Pvec,"distR"=>LNRvec,"distT"=>LNTvec,"distN"=>LNNvec)

end

################################## XYT FUNCTIONS ################################


### Finds nearest neighbors from xyt data - f64
function find_neighbors_pxyt(x::Vector{Float64}, y::Vector{Float64}, times::Vector{Float64}, 
    mags::Vector{Float64}, bval::Float64, dfrac::Float64)

    # setup magnitude scaling and other arrays
    scaLM = -bval*mags
    X = transpose(hcat(x,y)) # kilometers
    T = times/(86400.0) # days

    # add noise (+/-10m) for numerical stability
    X .+= 0.01*randn(size(X))

    # allocate arrays - output as Float64
    npts = length(x)
    LNNvec = zeros(Float64,npts)
    LNRvec = zeros(Float64,npts)
    LNTvec = zeros(Float64,npts)
    LRkm = zeros(Float64,npts)
    LTday = zeros(Float64,npts)
    Pvec = zeros(Int64,npts)

    # Pairwise Distances in float 64    
    LDX = dfrac*log10.(pairwise(Euclidean(),X)) # km

    # Pairwise Times in float 64
    LDT = log10.(pairwise(Euclidean(),T) .- log10(365.0)) # years

    # Calculate Nearest Neighbors
    @inbounds for jj = 2:npts
        @views idx = argmin(LDX[1:jj-1,jj].+LDT[1:jj-1,jj].+scaLM[1:jj-1]) # find parent
        Pvec[jj] = idx # parent index
        LNNvec[jj] = LDX[idx,jj]+LDT[idx,jj]+scaLM[idx] # nearest neighbor distance
        LNRvec[jj] = LDX[idx,jj]+0.5*scaLM[idx] # rescaled distance
        LNTvec[jj] = LDT[idx,jj]+0.5*scaLM[idx] # rescaled time
        LRkm[jj] = LDX[idx,jj]/dfrac # distance in km
        LTday[jj] = LDT[idx,jj] + log10(365.0) # time offset in days
    end

    # return as DataFrame
    return DataFrame("enum"=>1:npts,"pid"=>Pvec,"l10KM"=>LRkm,"l10DAY"=>LTday,
        "distR"=>LNRvec,"distT"=>LNTvec,"distN"=>LNNvec)

end


################################## XYZT FUNCTIONS ################################

### Finds nearest neighbors from xyzt data - f64
function find_neighbors_pxyzt(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64},
    times::Vector{Float64}, mags::Vector{Float64}, bval::Float64, dfrac::Float64)

    # setup magnitude scaling and other arrays
    scaLM = -bval*mags
    X = transpose(hcat(x,y,z)) # kilometers
    T = times/(86400.0) # days

    # add noise (+/-10m) for numerical stability
    X .+= 0.01*randn(size(X))

    # allocate arrays - output as Float64
    npts = length(x)
    LNNvec = zeros(Float64,npts)
    LNRvec = zeros(Float64,npts)
    LNTvec = zeros(Float64,npts)
    LRkm = zeros(Float64,npts)
    LTday = zeros(Float64,npts)
    Pvec = zeros(Int64,npts)

    # Pairwise Distances in float 64    
    LDX = dfrac*log10.(pairwise(Euclidean(),X)) # km

    # Pairwise Times in float 64
    LDT = log10.(pairwise(Euclidean(),T) .- log10(365.0)) # years

    # Calculate Nearest Neighbors
    @inbounds for jj = 2:npts
        @views idx = argmin(LDX[1:jj-1,jj].+LDT[1:jj-1,jj].+scaLM[1:jj-1]) # find parent
        Pvec[jj] = idx # parent index
        LNNvec[jj] = LDX[idx,jj]+LDT[idx,jj]+scaLM[idx] # nearest neighbor distance
        LNRvec[jj] = LDX[idx,jj]+0.5*scaLM[idx] # rescaled distance
        LNTvec[jj] = LDT[idx,jj]+0.5*scaLM[idx] # rescaled time
        LRkm[jj] = LDX[idx,jj]/dfrac # distance in km
        LTday[jj] = LDT[idx,jj] + log10(365.0) # time offset in days
    end

    # return as DataFrame
    return DataFrame("enum"=>1:npts,"pid"=>Pvec,"l10KM"=>LRkm,"l10DAY"=>LTday,
        "distR"=>LNRvec,"distT"=>LNTvec,"distN"=>LNNvec)

end
