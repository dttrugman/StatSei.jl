#### Functions for Nearest Neighbors Algorithm of Zaliapin and Ben-Zion, 2013
#    Daniel Trugman, Nevada Seismological Laboratory (2023)


################################## LLT FUNCTIONS ################################

### Calculate Earth Radius at a Reference Latitude
function calculate_erad(lat::Float64)

    # geographic constants, WGS84
    a = 6378.1370000000000 # major axis
    b = 6356.7523141999999 # minor axis

    # calculate
    r = sqrt( ( (a^2 *cosd(lat))^2 + (b^2 *sind(lat))^2 ) / (
        (a*cosd(lat))^2 + (b*sind(lat))^2 ) )

    # return in km
    return r

end

### Calculate expected rupture radius in km - circular crack model
function calculate_l10Rad(mag::Float64,stress::Float64=5.0) # stress in MPa
    return 0.5*mag + (log10(7.0/16.0) + 9.1 - log10(stress) - 6.0)/3.0 - 3.0 # log10 km 
end
function calculate_l10Rad(mag::Vector{Float64},stress::Float64=5.0) # stress in MPa
    return 0.5*mag .+ (log10(7.0/16.0) + 9.1 - log10(stress) - 6.0)/3.0 .- 3.0 # log10 km
end

### Calculate expected rupture length in km - Noda (2013)
function calculate_l10Len(mag::Float64,stress::Float64=5.0,maxW::Float64=20.0)
    M0 = 10.0^(1.5*mag+9.1) # calculate moments
    itp = LinearInterpolation([1.0,4.0,16.0],[2.53,3.02,5.21]) # Noda, 2013 C-factor vs Aspect Ratio
    ARvals = collect(1.0:0.1:16.0) # fine grid
    Cvals = itp(ARvals) # cvals for fine grid
    Avals = (Cvals*M0/(stress*1.0e6)).^(2.0/3.0) ./ 1.0e6 # area in km
    Wvals = sqrt.(Avals./ARvals) # widths for each aspect ratio, area
    Lvals = Avals./Wvals # lengths for these widths
    idx = argmin(abs.(Wvals.-maxW)) # select aspect ratio = 1 or the closest to W=maxW
    return log10(Lvals[idx]) # return in log10
end
# vectorized version of the above
function calculate_l10Len(mags::Vector{Float64},stress::Float64=5.0,maxW::Float64=20.0) 
    l10L = zeros(Float64,length(mags))
    for (ii,mag) in enumerate(mags)
        l10L[ii] = calculate_l10Len(mag,stress,maxW)
    end
    return l10L
end

### Wells-Coppersmith Rupture Lengths (RLD in Table 2A)
function calculate_l10LWC(mag::Float64,sof::Char='A')
    if sof=='S'
        l10L = -2.57 + 0.62*mag
    elseif sof=='R'
        l10L = -2.42 + 0.58*mag
    elseif sof=='N'
        l10L = -1.88 + 0.50*mag
    else
        l10L = -2.44 + 0.59*mag
    end
    return l10L
end
function calculate_l10LWC(mag::Vector{Float64},sof::Char='A')
    if sof=='S'
        l10L = -2.57 .+ 0.62*mag
    elseif sof=='R'
        l10L = -2.42 .+ 0.58*mag
    elseif sof=='N'
        l10L = -1.88 .+ 0.50*mag
    else
        l10L = -2.44 .+ 0.59*mag
    end
    return l10L
end



### Finds nearest neighbors from llt data - f64
function find_neighbors_llt(lons::Vector{Float64}, lats::Vector{Float64}, times::Vector{Float64}, 
    mags::Vector{Float64}, bval::Float64, dfrac::Float64, Tmax::Float64=100.0)

    # setup magnitude scaling and other arrays
    scaLM = -bval*mags
    X = hcat(lons,lats) # meters
    T = times/(86400.0) # days
    LTmax = log10(Tmax*365.0) # maximum time shift

    # add noise for numerical stability
    X .+= 0.00001*randn(size(X))

    # get earth radius
    akm = calculate_erad(median(lats))

    # allocation
    npts = length(lons)
    LRvec = zeros(Float64,npts)
    LTvec = zeros(Float64,npts)
    LRkm = zeros(Float64,npts)
    LTday = zeros(Float64,npts)
    LNRvec = zeros(Float64,npts)
    LNTvec = zeros(Float64,npts)
    Pvec = zeros(Int64,npts)
            
    # loop over child events
    @inbounds for jj = 2:npts
        
        # Distance and time calculations
        jx = 1 # start index for neighbors calculations
        @inbounds for ii in range(jj-1,1,step=-1)
            LTvec[ii] = log10(T[jj]-T[ii]) 
            LRvec[ii] = dfrac*log10(haversine((X[jj,1],X[jj,2]),(X[ii,1],X[ii,2]),akm))
            if LTvec[ii] >= LTmax; jx = ii; break; end # stop when parents are too early
        end
        
        # find nearest neighbor
        @views jdx = argmin(LTvec[jx:jj-1].+LRvec[jx:jj-1].+scaLM[jx:jj-1])
        idx = jdx + jx - 1

        # finalize results
        Pvec[jj] = idx # parent id
        LRkm[jj] = LRvec[idx]/dfrac # distance in km
        LNRvec[jj] = LRvec[idx] + 0.5*scaLM[idx] # corrected for magnitude
        LTday[jj] = LTvec[idx] # time offset in days
        LNTvec[jj] = LTvec[idx] - log10(365.0) + 0.5*scaLM[idx] # corrected for magnitude, converted to years
        
        # Progress
        if mod(jj,10000) == 0
            @printf("Completed %8d %8d\n",jj,npts)
        end

    end

    # return as DataFrame, log distances
    LNNvec = LNRvec .+ LNTvec
    return DataFrame("enum"=>1:npts,"pid"=>Pvec,"l10KM"=>LRkm,"l10DAY"=>LTday,
        "distR"=>LNRvec,"distT"=>LNTvec,"distN"=>LNNvec)

end

################################## XYT FUNCTIONS ################################


### Finds nearest neighbors from xyt data - f64
function find_neighbors_xyt(x::Vector{Float64}, y::Vector{Float64}, times::Vector{Float64}, 
    mags::Vector{Float64}, bval::Float64, dfrac::Float64, Tmax::Float64=100.0)

    # setup magnitude scaling and other arrays
    scaLM = -bval*mags
    X = hcat(x,y) # kilometers
    T = times/(86400.0) # days
    LTmax = log10(Tmax*365.0) # maximum time shift

    # add noise (+/-10m) for numerical stability
    X .+= 0.01*randn(size(X))

    # allocation
    npts = length(lons)
    LRvec = zeros(Float64,npts)
    LTvec = zeros(Float64,npts)
    LNRvec = zeros(Float64,npts)
    LNTvec = zeros(Float64,npts)
    LRkm = zeros(Float64,npts)
    LTday = zeros(Float64,npts)
    Pvec = zeros(Int64,npts)
            
    # loop over child events
    @inbounds for jj = 2:npts
        
        # Distance and time calculations
        jx = 1 # start index for neighbors calculations
        @inbounds for ii in range(jj-1,1,step=-1)
            LTvec[ii] = log10(T[jj]-T[ii]) 
            LRvec[ii] = dfrac*log10(euclidean((X[jj,1],X[jj,2]),(X[ii,1],X[ii,2])))
            if LTvec[ii] >= LTmax; jx = ii; break; end # stop when parents are too early
        end

        # find nearest neighbor
        @views jdx = argmin(LTvec[jx:jj-1].+LRvec[jx:jj-1].+scaLM[jx:jj-1])
        idx = jdx + jx - 1

        # finalize results
        Pvec[jj] = idx # parent id
        LRkm[jj] = LRvec[idx]/dfrac # distance in km
        LNRvec[jj] = LRvec[idx] + 0.5*scaLM[idx] # corrected for magnitude
        LTday[jj] = LTvec[idx] # time offset in days
        LNTvec[jj] = LTvec[idx] - log10(365.0) + 0.5*scaLM[idx] # corrected for magnitude, converted to years
        
        # Progress
        if mod(jj,10000) == 0
            @printf("Completed %8d %8d\n",jj,npts)
        end

    end

    # return as DataFrame, log distances
    LNNvec = LNRvec .+ LNTvec
    return DataFrame("enum"=>1:npts,"pid"=>Pvec,"l10KM"=>LRkm,"l10DAY"=>LTday,
        "distR"=>LNRvec,"distT"=>LNTvec,"distN"=>LNNvec)

end

################################## XYZT FUNCTIONS ################################

### Finds nearest neighbors from xyzt data - f64
function find_neighbors_xyzt(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, 
    times::Vector{Float64}, mags::Vector{Float64}, bval::Float64, dfrac::Float64, Tmax::Float64=100.0)

    # setup magnitude scaling and other arrays
    scaLM = -bval*mags
    X = hcat(x,y,z) # meters
    T = times/(86400.0) # days
    LTmax = log10(Tmax*365.0) # maximum time shift

    # add noise (+/-10m) for numerical stability
    X .+= 0.01*randn(size(X))

    # allocation
    npts = length(lons)
    LRvec = zeros(Float64,npts)
    LTvec = zeros(Float64,npts)
    LNRvec = zeros(Float64,npts)
    LNTvec = zeros(Float64,npts)
    LRkm = zeros(Float64,npts)
    LTday = zeros(Float64,npts)
    Pvec = zeros(Int64,npts)
            
    # loop over child events
    @inbounds for jj = 2:npts
        
        # Distance and time calculations
        jx = 1 # start index for neighbors calculations
        @inbounds for ii in range(jj-1,1,step=-1)
            LTvec[ii] = log10(T[jj]-T[ii]) 
            LRvec[ii] = dfrac*log10(euclidean((X[jj,1],X[jj,2],X[jj,3]),(X[ii,1],X[ii,2],X[[ii,3]])))
            if LTvec[ii] >= LTmax; jx = ii; break; end # stop when parents are too early
        end

        # find nearest neighbor
        @views jdx = argmin(LTvec[jx:jj-1].+LRvec[jx:jj-1].+scaLM[jx:jj-1])
        idx = jdx + jx - 1

        # finalize results
        Pvec[jj] = idx # parent id
        LRkm[jj] = LRvec[idx]/dfrac # distance in km
        LNRvec[jj] = LRvec[idx] + 0.5*scaLM[idx] # corrected for magnitude
        LTday[jj] = LTvec[idx] # time offset in days
        LNTvec[jj] = LTvec[idx] - log10(365.0) + 0.5*scaLM[idx] # corrected for magnitude, converted to years
        
        # Progress
        if mod(jj,10000) == 0
            @printf("Completed %8d %8d\n",jj,npts)
        end

    end

    # return as DataFrame, log distances
    LNNvec = LNRvec .+ LNTvec
    return DataFrame("enum"=>1:npts,"pid"=>Pvec,"l10KM"=>LRkm,"l10DAY"=>LTday,
        "distR"=>LNRvec,"distT"=>LNTvec,"distN"=>LNNvec)

end


######################## AUXILIARY FUNCTIONS #############################

### Shuffle catalog to obtain threshold: lat, lon, time
function get_thresh_shuffle_llt(lons::Vector{Float64}, lats::Vector{Float64}, 
    times::Vector{Float64}, mags::Vector{Float64}, bval::Float64, dfrac::Float64, 
    n_shuffle::Int64=10, pct_shuffle::Float64=5.0, xy_shuffle::String="permute", 
    Tmax::Float64=5.0)

    # sampling vector
    ndat = length(lons)
    idx = Vector{Int32}(1:ndat)

    # start and end time
    tmin, tmax = times[1], times[end]

    # region bounds
    lonmin, lonmax = minimum(lons), maximum(lons)
    latmin, latmax = minimum(lats), maximum(lats)

    # store data here
    logR = zeros(ndat*n_shuffle)
    logT = zeros(ndat*n_shuffle)

    # loop over bootstraps resamplings
    Threads.@threads for ib in 1:n_shuffle

        # start
        @printf("Starting with iteration %d/%d\n",ib,n_shuffle)

        # shuffle data
        isamp = sample(idx,ndat,replace=false)
        bmags = mags[isamp]
        if xy_shuffle == "permute"
            isamp = sample(idx,ndat,replace=false)
            blons = lons[isamp]
            blats = lats[isamp]
        else
            blons = lonmin .+ (lonmax-lonmin)*rand(ndat)
            blats = latmin .+ (latmax-latmin)*rand(ndat)
        end
        btimes = sort(tmin .+ (tmax-tmin)*rand(ndat))

        # input to neighbor algorithm
        ndf = find_neighbors_llt(
            blons, blats, btimes, bmags, bval, dfrac, Tmax
        )

        # save results
        idx1 = 1 + (ib-1)*ndat
        idx2 = idx1 + ndat - 1
        logT[idx1:idx2] .= ndf.distT
        logR[idx1:idx2] .= ndf.distR

        # done with this iteration
        @printf("Done with iteration %d/%d\n",ib,n_shuffle)

    end

    # finalize
    logN = logT .+ logR
    ndf = DataFrame("distR"=>logR,"distT"=>logT,"distN"=>logN)

    # return
    thresh = percentile(logN,pct_shuffle)
    return thresh, ndf

end

### Shuffle catalog to obtain threshold: x, y, time
function get_thresh_shuffle_xyt(x::Vector{Float64}, y::Vector{Float64}, 
    times::Vector{Float64}, mags::Vector{Float64}, bval::Float64, dfrac::Float64, 
    n_shuffle::Int64=10, pct_shuffle::Float64=5.0, xy_shuffle::String="permute", 
    Tmax::Float64=5.0)

    # sampling vector
    ndat = length(lons)
    idx = Vector{Int32}(1:ndat)

    # start and end time
    tmin, tmax = times[1], times[end]

    # region bounds
    xmin, xmax = minimum(x), maximum(x)
    ymin, ymax = minimum(y), maximum(y)

    # store data here
    logR = zeros(ndat*n_shuffle)
    logT = zeros(ndat*n_shuffle)

    # loop over bootstraps resamplings
    Threads.@threads for ib in 1:n_shuffle

        # start
        @printf("Starting with iteration %d/%d\n",ib,n_shuffle)

        # shuffle data
        isamp = sample(idx,ndat,replace=false)
        bmags = mags[isamp]
        if xy_shuffle == "permute"
            isamp = sample(idx,ndat,replace=false)
            bx = x[isamp]
            by = y[isamp]
        else
            bx = xmin .+ (xmax-xmin)*rand(ndat)
            by = ymin .+ (ymax-ymin)*rand(ndat)
        end
        btimes = sort(tmin .+ (tmax-tmin)*rand(ndat))

        # input to neighbor algorithm
        ndf = find_neighbors_xyt(
            bx, by, btimes, bmags, bval, dfrac, Tmax
        )

        # save results
        idx1 = 1 + (ib-1)*ndat
        idx2 = idx1 + ndat - 1
        logT[idx1:idx2] .= ndf.distT
        logR[idx1:idx2] .= ndf.distR

        # done with this iteration
        @printf("Done with iteration %d/%d\n",ib,n_shuffle)

    end

    # finalize
    logN = logT .+ logR
    ndf = DataFrame("distR"=>logR,"distT"=>logT,"distN"=>logN)

    # return
    thresh = percentile(logN,pct_shuffle)
    return thresh, ndf

end

### Shuffle catalog to obtain threshold: x, y, time
function get_thresh_shuffle_xyzt(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64},
    times::Vector{Float64}, mags::Vector{Float64}, bval::Float64, dfrac::Float64, 
    n_shuffle::Int64=10, pct_shuffle::Float64=5.0, xy_shuffle::String="permute", 
    Tmax::Float64=5.0)

    # sampling vector
    ndat = length(lons)
    idx = Vector{Int32}(1:ndat)

    # start and end time
    tmin, tmax = times[1], times[end]

    # region bounds
    xmin, xmax = minimum(x), maximum(x)
    ymin, ymax = minimum(y), maximum(y)
    zmin, zmax = minimum(z), maximum(z)

    # store data here
    logR = zeros(ndat*n_shuffle)
    logT = zeros(ndat*n_shuffle)

    # loop over bootstraps resamplings
    Threads.@threads for ib in 1:n_shuffle

        # start
        @printf("Starting with iteration %d/%d\n",ib,n_shuffle)

        # shuffle data
        isamp = sample(idx,ndat,replace=false)
        bmags = mags[isamp]
        if xy_shuffle == "permute"
            isamp = sample(idx,ndat,replace=false)
            bx = x[isamp]
            by = y[isamp]
            bz = z[isamp]
        else
            bx = xmin .+ (xmax-xmin)*rand(ndat)
            by = ymin .+ (ymax-ymin)*rand(ndat)
            bz = zmin .+ (zmax-zmin)*rand(ndat)
        end
        btimes = sort(tmin .+ (tmax-tmin)*rand(ndat))

        # input to neighbor algorithm
        ndf = find_neighbors_xyzt(
            bx, by, bz, btimes, bmags, bval, dfrac, Tmax
        )

        # save results
        idx1 = 1 + (ib-1)*ndat
        idx2 = idx1 + ndat - 1
        logT[idx1:idx2] .= ndf.distT
        logR[idx1:idx2] .= ndf.distR

        # done with this iteration
        @printf("Done with iteration %d/%d\n",ib,n_shuffle)

    end

    # finalize
    logN = logT .+ logR
    ndf = DataFrame("distR"=>logR,"distT"=>logT,"distN"=>logN)

    # return
    thresh = percentile(logN,pct_shuffle)
    return thresh, ndf

end


### Gaussian mixture model to obtain threshold seperating modes in distance distribution
function fit_gmm(LNvec::Vector{Float64},ratio=1.0)

    # data to fit
    logNN = LNvec[isfinite.(LNvec)]

    # fit with default options
    gmix = GMM(2,logNN)
    g1 = Normal(gmix.μ[1],sqrt(gmix.Σ[1]))
    g2 = Normal(gmix.μ[2],sqrt(gmix.Σ[2]))
    
    # setup interpolate to find crossover point
    if gmix.μ[1] < gmix.μ[2]
        x1, x2 = round(gmix.μ[1],digits=2), round(gmix.μ[2],digits=2)
        xpts = collect(range(x1,x2,step=0.01))
        ypts = (gmix.w[2]*pdf.(g2,xpts))./(gmix.w[1]*pdf.(g1,xpts))
    else
        x1, x2 = round(gmix.μ[2],digits=2), round(gmix.μ[1],digits=2)
        xpts = collect(range(x1,x2,step=0.01))
        ypts = (gmix.w[1]*pdf.(g1,xpts))./(gmix.w[2]*pdf.(g2,xpts))
    end

    # point where probability ratio hits target value
    itp = linear_interpolation(ypts,xpts)
    thresh = itp(ratio)

    # return
    return thresh, gmix

end

### Given neighbor distances, assign families
function assign_families(LNvec::Vector{Float64}, Pvec::Vector{Int64},
    distLR::Vector{Float64},distLT::Vector{Float64}, # log10 km, log10 days
    mags::Vector{Float64},thresh::Float64, # threshold for nearest neighbor distance
    Rscale::Float64,Tscale::Float64) # parameters for remote triggering checks   

    # initialize families
    npts = length(LNvec)
    fnum = collect(1:npts)
    qidx = collect(1:npts)

    # calculate expected rupture dimensions, log10 km
    #l10R = calculate_l10Rad(mags) # radius
    #l10R = calculate_l10Len(mags) # length
    l10R = calculate_l10LWC(mags) # length, wells-coppersmith
    l10R = ifelse.(l10R.>=1.0,l10R,1.0) # set a floor for this

    # distance and time thresholds
    maxLR = log10(Rscale) .+ l10R # km
    maxLT = log10.(mags*Tscale) # days

    # loop over all events and reassign
    count = 0
    for ii in 2:npts
        
        # parent index
        jj = Pvec[ii]

        # assign to parent cluster under following conditions
        #    - nearest-neighbor distance < threshold AND
        #    - geographic distance in km  < Rscale * parent rupture dimension (or uncertainty floor) OR
        #    - temporal difference in days < Tscale * parent mag
        if (LNvec[ii] < thresh) & ((distLR[ii]<maxLR[jj]) | (distLT[ii]<maxLT[jj]))
            fnum[ii] = fnum[jj]
        elseif (LNvec[ii] < thresh) & ((distLR[ii]>maxLR[jj]) & (distLT[ii]>maxLT[jj]))
            count += 1
        end
    end 
    @printf("Family assignments prevented by RT criteria: %d\n",count)

    # create temporary dataframe with event and cluster number
    tdf = DataFrame("enum"=>qidx,"fnum"=>fnum)

    # compute nbranch = number of events in a cluster
    transform!(groupby(tdf, :fnum), nrow => :nb)

    # unique clusters, sort by nbranch
    select!(tdf,[:fnum,:nb])
    unique!(tdf)
    sort!(tdf,[:nb,:fnum],rev=[true,false])

    # assign new cluster ids, starting with largest
    tdf[!,:fid] =range(1,nrow(tdf),step=1)

    # add back in event ids, reorder by event index
    fdf = innerjoin(DataFrame("enum"=>qidx,"fnum"=>fnum),tdf,on=:fnum)
    sort!(fdf,:enum)
    select!(fdf,[:enum,:fnum,:fid,:nb])

    # return
    return fdf

end

### Plots results with GMM
function plot_neighbors(Rvec::Vector{Float64},Tvec::Vector{Float64},Nvec::Vector{Float64},
    thresh::Float64,gmix::GMM{Float64})

    # figure setup
    fig, axi = plt.subplots(1,2,figsize=(10,5),layout="tight")

    # data points to plot
    idx = isfinite.(Nvec)

    # plot histogram
    ax = axi[1]
    xx = collect(range(-15.0,0.0,step=0.2))
    ax.set_facecolor("gainsboro")
    ax.hist(Nvec[idx],xx,color="black",density=true)
    ax.set_xlabel("Nearest Neighbor Distance [logN]",fontsize=14)
    ax.set_ylabel("Density",fontsize=14)
    ax.set_xticks(range(-15,0,step=5))
    yL = ax.get_ylim()
    ax.set_yticks(range(0.0,yL[2],step=0.1))

    # add in mixture model
    g1 = Normal(gmix.μ[1],sqrt(gmix.Σ[1]))
    g2 = Normal(gmix.μ[2],sqrt(gmix.Σ[2]))
    ax.plot(xx,gmix.w[1]*pdf.(g1,xx),"--",color="deepskyblue",label="Gaussian Mixture")
    ax.plot(xx,gmix.w[2]*pdf.(g2,xx),"--",color="deepskyblue")
    yL = ax.get_ylim()
    ax.plot([thresh,thresh],[yL[1],yL[2]],"-r",lw=2,label=@sprintf("threshold = %.2f",thresh))
    ax.legend(loc="upper left",fontsize=12)

    # plot joint distribution
    ax = axi[2]
    sc = ax.hexbin(Tvec[idx],Rvec[idx],cmap="magma",gridsize=60)
    xL = [-8.5,-0.0]
    yL = thresh.-xL
    ax.plot(xL,thresh.-xL,"-w")
    ax.set_xlabel("Rescaled Time [logT]",fontsize=14)
    ax.set_ylabel("Rescaled Distance [logR]",fontsize=14)
    ax.set_xlim(xL[1],xL[2])
    ax.set_ylim(yL[2],yL[1])


    # return
    return fig

end

### Plots results without GMM
function plot_neighbors(Rvec::Vector{Float64},Tvec::Vector{Float64},Nvec::Vector{Float64},thresh::Float64)

    # figure setup
    fig, axi = plt.subplots(1,2,figsize=(10,5),layout="tight")

    # data points to plot
    idx = isfinite.(Nvec)

    # plot histogram
    ax = axi[1]
    xx = collect(range(-15.0,0.0,step=0.2))
    ax.set_facecolor("gainsboro")
    ax.hist(Nvec[idx],xx,color="black",density=true)
    ax.set_xlabel("Nearest Neighbor Distance [logN]",fontsize=14)
    ax.set_ylabel("Density",fontsize=14)
    ax.set_xticks(range(-15,0,step=5))
    yL = ax.get_ylim()
    ax.set_yticks(range(0.0,yL[2],step=0.1))

    # add in threshold
    yL = ax.get_ylim()
    ax.plot([thresh,thresh],[yL[1],yL[2]],"-r",lw=2,label=@sprintf("threshold = %.2f",thresh))
    ax.legend(loc="upper left",fontsize=12)

    # plot joint distribution
    ax = axi[2]
    sc = ax.hexbin(Tvec[idx],Rvec[idx],cmap="magma",gridsize=60)
    xL = [-8.5,-0.0]
    yL = thresh.-xL
    ax.plot(xL,thresh.-xL,"-w")
    ax.set_xlabel("Rescaled Time [logT]",fontsize=14)
    ax.set_ylabel("Rescaled Distance [logR]",fontsize=14)
    ax.set_xlim(xL[1],xL[2])
    ax.set_ylim(yL[2],yL[1])


    # return
    return fig

end