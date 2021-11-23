import Seis, Indicators, FFTW, Plots, Statistics

function hvplot((hvresult, hvfreq), name)
    Plots.plot(hvfreq, hvresult,
    label = "H/V Ratio",
    title = string("H/V of ", name),
    xlabel="Frequency (Hz)",
    ylabel="H/V Ratio",
    yscale=:log10,
    xscale=:log10,
    xminorgrid=true,
    xlims = (-Inf, 40))
    #Plots.plot!(hvfreq[5:end], Indicators.sma(hvresult, n=5)[5:end])
    Plots.savefig(string(name, ".png"))
    return nothing
end

function movingaverage(arr, mvlength)
    averagearray = []
    for v in arr
        push!(averagearray, Indicators.sma(Seis.trace(v), n=mvlength))
        deleteat!(averagearray[end], 1:mvlength-1)
    end
    return averagearray
end

function fft(arr)
    hcat([broadcast(abs, FFTW.fft(v)) for v in arr]...)
end


# Input 3 .sac files and 
function hv(arr, timedelta)
    hvresult = .√(view(arr, :, 1).^2 + view(arr, :, 2).^2) ./ view(arr, :, 3)
    hvfreq = [1:length(hvresult);] ./ (timedelta * length(hvresult))
    hvresult, hvfreq
end

function findlocalmaxhv((hvresult, hvfreq))
    index01hz = 0 # index of the maximum frequency smaller than or equal to 0.1Hz
    index10hz = 0 # index of the maximum frequency smaller than or equal to 10Hz
    for i in 1:length(hvfreq)
        if hvfreq[i] <= 5
            index10hz = i
            if hvfreq[i] <= 0.1
                index01hz = i
            end
        else
            break
        end
    end
    findmaxhv((hvresult[index01hz:index10hz], hvfreq[index01hz:index10hz]))
end

function findmaxhv((hvresult, hvfreq))
    index = findmax(hvresult)
    return hvresult[index[2]], hvfreq[index[2]]
end

function floortoheight(floor)
    if floor == "B3"
        return 0
    elseif floor == "B2"
        return 3
    elseif floor == "B1"
        return 6
    else
        return (parse(Int, floor) + 2) * 3
    end
end

function plotfloorhv(heighthvfreqarray)
    Plots.plot(view(heighthvfreqarray, :, 1), view(heighthvfreqarray, :, 2),
    title = "Height to H/V",
    xlabel="Height (m)",
    ylabel="H/V Ratio",
    seriestype = :scatter,
    ylim = (0, 1000))
    averageerror = Array{Float64}(undef, 0, 3)
    for i in 0:3:39
        indexes = findall(x->x==i, view(heighthvfreqarray, :, 1))
        #=while true
            mean = Statistics.mean(view(heighthvfreqarray, :, 2)[indexes])
            std = Statistics.std(view(heighthvfreqarray, :, 2)[indexes])
            newindex = findall(x->abs(x-mean)<=(std*3),view(heighthvfreqarray, :, 2)[indexes])
            if newindex == indexes
                averageerror = vcat(averageerror, [i mean std])
                break
            end
            indexes = newindex
        end=#
        mean = Statistics.mean(view(heighthvfreqarray, :, 2)[indexes])
        std = Statistics.std(view(heighthvfreqarray, :, 2)[indexes])
        averageerror = vcat(averageerror, [i mean std])
    end
    Plots.plot!(view(averageerror, :, 1), view(averageerror, :, 2), yerror=view(averageerror, :, 3))
    Plots.savefig(string("FloortoHV", ".png"))
    return nothing
end

function plotfloorfreq(heighthvfreqarray)
    Plots.plot(view(heighthvfreqarray, :, 1), view(heighthvfreqarray, :, 3),
    title = "Height to Frequency",
    xlabel="Height (m)",
    ylabel="HFrequency (Hz)",
    seriestype = :scatter)
    Plots.savefig(string("FloortoFreq", ".png"))
    return nothing
end

# MAIN PROGRAM

movingaveragelength = 5
startpath = pwd()

heighthvfreqarray = Array{Float64}(undef, 0, 3)

for (root, dirs, files) in walkdir(".")
    println("CD $root")
    quakefilename = [filter(x->occursin(r".*BN1.*\.sac", x), files),
        filter(x->occursin(r".*BN2.*\.sac", x), files),
        filter(x->occursin(r".*BN3.*\.sac", x), files)]
    
    if !(isempty(quakefilename[1]) || isempty(quakefilename[2]) || isempty(quakefilename[3]))
        cd(root)
        seisarray = hcat([Seis.read_sac(v[1]) for v in quakefilename]...);
        hvdata = hv(fft(movingaverage(seisarray, movingaveragelength)), seisarray[1].delta)
        cd(startpath)

        hvplot(hvdata, splitpath(root)[end])
        localmaxhv = findlocalmaxhv(hvdata)
        global heighthvfreqarray = vcat(heighthvfreqarray, [floortoheight(match(r"(?<floor>^[B]?[1-9][0-9]*)",splitpath(root)[end])[:floor]) localmaxhv[1] localmaxhv[2]])
    end
end
plotfloorhv(heighthvfreqarray)
plotfloorfreq(heighthvfreqarray)