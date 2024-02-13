### Customized function to read the scsn catalog
function load_catalog_scsn(fcat,minmag)

    # load file
    df = DataFrame(CSV.File(fcat,delim=' ',ignorerepeated=true,
        select=["DATE","TIME","MAG","LAT","LON","DEPTH","EVID"]))
    rename!(df,Dict(zip(["EVID","LAT","LON","DEPTH","MAG"],["evid","lat","lon","dep","mag"])))

    # select rows
    df = df[df[!,:mag].>=minmag,:]

    # fix date times
    otimes = []
    for ii in 1:nrow(df)
        y,m,d = split(df[ii,:DATE],"/")
        year, month, day = parse(Int64,y),parse(Int64,m),parse(Int64,d)
        h, m, sm = split(df[ii,:TIME],":")
        s,ms = split(sm,".")
        hour, minute, second = parse(Int64,h),parse(Int64,m),parse(Int64,s)
        msec = 10*parse(Int64,ms)
        if second >= 60
            second-=60
            minute+=1
        end
        push!(otimes,DateTime(year,month,day,hour,minute,second,msec))
    end

    # origin time
    df[!,:otime] .= otimes

    # unix datetime
    df[!,:tepoch] = datetime2unix.(df[!,:otime])

    # sort by time
    sort!(df,:tepoch)

    # duplicates
    ikeep = df.tepoch .> 0.0
    for ii in 2:nrow(df)
        if (df.tepoch[ii]-df.tepoch[ii-1])<0.1
            ikeep[ii-1] = false
        end
    end
    println("Duplicates removed, keeping: ",sum(ikeep), "/",nrow(df))
    df = df[ikeep,:]

    # compile final catalog
    select!(df,[:evid,:otime,:tepoch,:lat,:lon,:dep,:mag])
    unique!(df)
    df[!,:enum] = 1:nrow(df)

    # return
    return df

end


### Customized function to read the relocated socal catalog from 2011
function load_catalog_hs2011(fcat,minmag)

    # define columns to use
    acols = ["year","month","day","hour","minute","second","evid",
    "lat","lon","dep","mag","nph","csta","rms","dflag","lflag","cid","nb","ndt",
        "aherr","azerr","rherr","rzerr","et","lt","poly"]
    scols = ["year","month","day","hour","minute","second","evid","mag","lat","lon","dep"]

    # load file
    df = DataFrame(CSV.File(fcat,delim=' ',silencewarnings=true,
        ignorerepeated=true,header=acols,select=scols))
        
    # select rows
    df = df[df[!,:mag].>=minmag,:]

    # fix timing for select events
    for ii in 1:nrow(df)
        if df[ii,:minute]>=60
            df[ii,:minute]-=60
            df[ii,:hour]+=1
        elseif df[ii,:minute]<0
            df[ii,:minute]+=60
            df[ii,:hour]-=1
        end
        if df[ii,:second]>=60.0
            df[ii,:second]-=60.0
            df[ii,:minute]+=1
        elseif df[ii,:second]<0
            df[ii,:second]+=60
            df[ii,:minute]-=1
        end
    end

    # compute timing
    df[!,:otime] = DateTime.(df[!,:year],df[!,:month],df[!,:day],
            df[!,:hour],df[!,:minute],convert.(Int64,floor.(df[!,:second])), 
            convert.(Int64,round.(1000.0*(df[!,:second].-floor.(df[!,:second])))))

    # unix datetime
    df[!,:tepoch] = datetime2unix.(df[!,:otime])

    # sort by time
    sort!(df,:tepoch)

    # duplicates
    ikeep = df.tepoch .> 0.0
    for ii in 2:nrow(df)
        if (df.tepoch[ii]-df.tepoch[ii-1])<0.1
            ikeep[ii-1] = false
        end
    end
    println("Duplicates removed, keeping: ",sum(ikeep), "/",nrow(df))
    df = df[ikeep,:]

    # compile final catalog
    select!(df,[:evid,:otime,:tepoch,:lat,:lon,:dep,:mag])
    unique!(df)
    df[!,:enum] = 1:nrow(df)

    # return
    return df
end

### Customized function to read the relocated socal catalog from 2022
function load_catalog_hs2022(fcat,minmag)

    # define columns to use
    scols = ["year","month","day","hour","minute","second","evid","lat","lon","dep","mag"]
    ncol = length(scols)

    # load file
    df = DataFrame(CSV.File(fcat,delim=' ',silencewarnings=true,
        ignorerepeated=true,header=false,select=range(1,ncol)))
    rename!(df,Dict(zip(1:ncol,scols)))
        
    # select rows
    df = df[df[!,:mag].>=minmag,:]

    # fix timing for select events
    for ii in 1:nrow(df)
        if df[ii,:minute]>=60
            df[ii,:minute]-=60
            df[ii,:hour]+=1
        elseif df[ii,:minute]<0
            df[ii,:minute]+=60
            df[ii,:hour]-=1
        end
        if df[ii,:second]>=60.0
            df[ii,:second]-=60.0
            df[ii,:minute]+=1
        elseif df[ii,:second]<0
            df[ii,:second]+=60
            df[ii,:minute]-=1
        end
    end

    # compute timing
    df[!,:otime] = DateTime.(df[!,:year],df[!,:month],df[!,:day],
            df[!,:hour],df[!,:minute],convert.(Int64,floor.(df[!,:second])), 
            convert.(Int64,round.(1000.0*(df[!,:second].-floor.(df[!,:second])))))

    # unix datetime
    df[!,:tepoch] = datetime2unix.(df[!,:otime])

    # sort by time
    sort!(df,:tepoch)

    # duplicates
    ikeep = df.tepoch .> 0.0
    for ii in 2:nrow(df)
        if (df.tepoch[ii]-df.tepoch[ii-1])<0.1
            ikeep[ii-1] = false
        end
    end
    println("Duplicates removed, keeping: ",sum(ikeep), "/",nrow(df))
    df = df[ikeep,:]

    # compile final catalog
    select!(df,[:evid,:otime,:tepoch,:lat,:lon,:dep,:mag])
    unique!(df)
    df[!,:enum] = 1:nrow(df)

    # return
    return df
end

### Customized function to read the ncsn catalog
function load_catalog_ncsn(fcat,minmag)

    # load file -only select columns to avoid mangling...
    df = DataFrame(CSV.File(fcat,delim=' ',silencewarnings=true,ignorerepeated=true,
     select=["Date","Time","Lat","Lon","Depth","Mag"], # avoid mangled columns
     types=Dict("Time"=>String,"Lat"=>Float64,"Lon"=>Float64,"Depth"=>Float64)))

    # get event ids
    f = open(fcat,"r")
    evids = zeros(Int64,nrow(df))
    for (ii, line) in enumerate(eachline(f))
        if ii == 1; continue; end
        evids[ii-1] = parse(Int64,line[88:end])
    end
    close(fcat)
    df[!,:EventID] .= evids
        
    # select rows
    df = df[df[!,:Mag].>=minmag,:]

    # datetime from string
    df[!,:Time] .*= "0"
    df[!,:otime] = DateTime.(df.Date .* " " .* df.Time,dateformat"Y/m/d H:M:S.s")

    # unix datetime
    df[!,:tepoch] = datetime2unix.(df[!,:otime])

    # sort by time
    sort!(df,:tepoch)

    # duplicates
    ikeep = df.tepoch .> 0.0
    for ii in 2:nrow(df)
        if (df.tepoch[ii]-df.tepoch[ii-1])<0.1
            ikeep[ii-1] = false
        end
    end
    println("Duplicates removed, keeping: ",sum(ikeep), "/",nrow(df))
    df = df[ikeep,:]

    # compile final catalog
    select!(df,[:EventID,:otime,:tepoch,:Lat,:Lon,:Depth,:Mag])
    unique!(df)
    df[!,:enum] = 1:nrow(df)
    rename!(df,Dict(zip(["EventID","Lat","Lon","Depth","Mag"],["evid","lat","lon","dep","mag"])))

    # return
    return df
end

### Customized function to read nsl compilation
function load_catalog_nsl(fcat,minmag)

    # load
    df = DataFrame(CSV.File(fcat,delim=','))

    # select rows
    df = df[df[!,:mag].>=minmag,:]

    # datetime
    df[!,:otime] = DateTime.(SubString.(df.otime,1,21),dateformat"Y-m-d H:M:S.s")

     # unix datetime
     df[!,:tepoch] = datetime2unix.(df[!,:otime])

     # sort by time
     sort!(df,:tepoch)
 
     # duplicates
     ikeep = df.tepoch .> 0.0
     for ii in 2:nrow(df)
         if (df.tepoch[ii]-df.tepoch[ii-1])<0.1
             ikeep[ii-1] = false
         end
     end
     println("Duplicates removed, keeping: ",sum(ikeep), "/",nrow(df))
     df = df[ikeep,:]

    # compile final catalog
    select!(df,[:evid,:otime,:tepoch,:lat,:lon,:dep,:mag])
    unique!(df)
    df[!,:enum] = 1:nrow(df)

    # return
    return df

end

### Customized function to read merged CA/NV (preprocessed)
function load_catalog_merged(fcat,minmag)
    
    # load
    df = DataFrame(CSV.File(fcat,delim=','))

    # select rows
    df = df[df[!,:mag].>=minmag,:]

    # datetime
    df[!,:otime] = DateTime.(df.otime,dateformat"Y-m-d H:M:S.s")

    # unix datetime
    df[!,:tepoch] = datetime2unix.(df[!,:otime])

    # sort by time
    sort!(df,:tepoch)

    # compile final catalog
    select!(df,[:evid,:otime,:tepoch,:lat,:lon,:dep,:mag])
    unique!(df)
    df[!,:enum] = 1:nrow(df)

    # return
    return df

end

### Customized function to read relocated NV catalog
function load_catalog_nvreloc(fcat,minmag)
    
    # load
    df = DataFrame(CSV.File(fcat,delim=' '))

    # select rows
    df = df[df[!,:mag].>=minmag,:]

    # datetime
    df[!,:otime] = DateTime.(df.otime,dateformat"Y-m-d H:M:S.s")

    # unix datetime
    df[!,:tepoch] = datetime2unix.(df[!,:otime])

    # sort by time
    sort!(df,:tepoch)

    # compile final catalog
    select!(df,[:evid,:otime,:tepoch,:lat,:lon,:dep,:mag])
    unique!(df)
    df[!,:enum] = 1:nrow(df)

    # return
    return df

end