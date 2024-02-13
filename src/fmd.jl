##### Function for FMD and B-Value Analysis
# Daniel Trugman, 2024

### Frequency-Magnitude Distribution
#  inputs: vector of magnitudes, magnitude increment (e.g., 0.1, 0.01)
#  returns: event counts and associated magnitude bin centers
function fmd(mags::Vector{Float64}, dmag::Float64)

    # compute upper and lower bounds
    mmin, mmax = minimum(mags), maximum(mags)

    # find bin edges, given these values and target increment
    mround = -Int64(log10(dmag))
    mbinL = round(mmin,digits=mround) - dmag/2.0
    mbinH = round(mmax,digits=mround) + dmag/2.0

    # adjust these bins to fit data bounds, if necessary
    if mbinL < mmin - dmag; mbinL += dmag; end
    if mbinH > mmax + dmag; mbinH -= dmag; end

    # now fit a histogram
    mbinE = collect(mbinL:dmag:mbinH)
    h = fit(Histogram,mags,mbinE)

    # return results
    counts = h.weights # event counts
    mbinC = 0.5*(mbinE[1:end-1] .+ mbinE[2:end]) # bin centers
    return counts, mbinC
end


### Maximum Curvature Estimate of Mc 
#  (recommended to apply shift later to correct bias)
# mags: vector of magnitudes
# mprec: magnitude precision (0.1, 0.01, etc)
function mc_maxc(mags::Vector{Float64},mprec::Float64=0.1)

    # get fmd
    counts, mbins = fmd(mags, mprec)

    # find peak
    imax = argmax(counts)

    # return this value
    return mbins[imax]
end


#### Maximum Likelihood B-value Estimator (Bender, 1983) ####
# mags: vector of magnitudes
# mprec: magnitude precision (0.1, 0.01, etc)
# mc: magnitude of completeness
function bval_maxl(mags::Vector{Float64},mprec::Float64,mc::Float64)

    # calculate mean magnitude
    mean_mag = mean(mags[mags.>=mc])

    # estimate b-value
    return 1.0 / (log(10.0)*(mean_mag-(mc-mprec/2.0)))

end

#### B-positive implementation from EQ 6 of Main Text
# dmags: vector of differential magnitudes
# dmc: differential magnitude of completeness 
function bval_bpos(dmags::Vector{Float64},dmc::Float64)
    
    # calculate mean magnitude above dmc
    mean_dmag = mean(dmags[dmags.>=dmc])
    
    # return b-value
    beta = 1.0/(mean_dmag - dmc)
    return beta / log(10.0)
end

#### B-positive implementation from Eq A3 in van der elst 2021
# dmags: vector of differential magnitudes
# mprec: magnitude precision (0.1, 0.01, etc)
# dmc: differential magnitude of completeness 
function bval_bpos(dmags::Vector{Float64},mprec::Float64,dmc::Float64)
    
    # calculate mean magnitude
    mean_dmag = mean(dmags[dmags.>=dmc])
    
    # argument for inverse coth
    delta = mprec/2.0
    x = (1.0/delta)*(mean_dmag-dmc+delta)
    
    # maximum likelihood estimator
    coth_inv = 0.5*(log((x+1.0)/(x-1.0)))
    beta = (1.0/delta)*coth_inv
    
    # return b-value
    b = beta/log(10.0)
    return b
end