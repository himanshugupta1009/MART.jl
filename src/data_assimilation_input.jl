using CSV, DataFrames, Dates
using DataStructures: OrderedDict

"""
Write a batch of samples to CSV.
Inputs are equal-length vectors.
t can be DateTime or Float64 (POSIX seconds).
"""
function write_samples_csv(path::AbstractString;
    t::AbstractVector,
    x::AbstractVector, y::AbstractVector, z::AbstractVector,
    P::AbstractVector, T::AbstractVector,
    U::AbstractVector, V::AbstractVector, W::AbstractVector,
    overwrite::Bool = true,
)
    n = length(t)
    @assert n == length(x) == length(y) == length(z) == length(P) == length(T) ==
            length(U) == length(V) == length(W)

    # normalize time to POSIX float seconds
    tposix = map(t) do ti
        ti isa DateTime ? Dates.datetime2unix(ti) : float(ti)
    end

    # fixed column order
    cols = (:t, :x, :y, :z, :P, :T, :U, :V, :W)

    # ordered units
    units = OrderedDict(
        :t => "seconds",
        :x => "km",
        :y => "km",
        :z => "km",
        :P => "Pa",
        :T => "K",
        :U => "m/s",
        :V => "m/s",
        :W => "m/s",
    )

    df = DataFrame(
        t = tposix,
        x = Float64.(x), y = Float64.(y), z = Float64.(z),
        P = Float64.(P), T = Float64.(T),
        U = Float64.(U), V = Float64.(V), W = Float64.(W),
    )

    # --- Write CSV with units header ---
    newfile = overwrite || !isfile(path)
    open(path, newfile ? "w" : "a") do io
        if newfile
            # 1) units comment
            println(io, "# Units: ", join([string(sym, "=", units[sym]) for sym in cols], ", "))
            # 2) header row (explicit)
            println(io, join(cols, ","))
        end
        # 3) data rows (no header)
        CSV.write(io, df; append=true, writeheader=false)
    end

end


#=
using Random

path = "data_assimilation_input_samples.csv"
t = [(8*3600) + 1800+i*15 for i in 0:180] #Similar to the format in nature run
start_x = 8.0 #8 kms
x = [start_x + 0.2*i for i in 0:180]
start_y = 8.0 #8 kms
y = [start_y + 0.2*i for i in 0:180]
start_z = 1.0 #1.0 kms
z = [start_z + 0.0*i for i in 0:180]
rng = MersenneTwister(42)
P = [900.0 + randn(rng)*20 for i in 0:180]
T = [300 + randn(rng)*5 for i in 0:180]
U = [7.0 + randn(rng)*2 for i in 0:180]
V = [9.0 + randn(rng)*2 for i in 0:180]
W = [0.0 + randn(rng)*0.01 for i in 0:180]

write_samples_csv(path;
    t = t,
    x = x, y = y, z = z,
    P = P, T = T,
    U = U, V = V, W = W,
    overwrite=true,
)

=#