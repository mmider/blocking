using CSV, DataFrames, Statistics, JLD2, MCMCDiagnostics

OUTPUT = joinpath("/", "home", "ebullient-llama", "Desktop", "output")

T2xT = Dict(
    0.2 => 0.1,
    0.4 => 0.85,
    0.5 => 0.85,
    1.0 => 0.95,
    2.0 => 2.5,
    4.0 => 4.85,
)

function load_file(filename, T)
    filename = joinpath("T_$(T)", filename)
    file = joinpath(OUTPUT, filename)
    isfile(file) && return CSV.File(file) |> DataFrame
    nothing
end

function load_elapsed(num_knots, T, xT = T2xT[T])
    filename = (
        num_knots == 0 ?
        "elapsed_no_blocking_T_$(T)_xT_$(xT).csv" :
        "elapsed_blocking_num_knots_$(num_knots)_T_$(T)_xT_$(xT).csv"
    )
    load_file(filename, T)
end

function load_midpt_data(num_knots, T, xT = T2xT[T])
    filename = "end_pts_blocking_num_knots_$(num_knots)_T_$(T)_xT_$(xT).csv"
    load_file(filename, T)
end

all_Ts = sort(collect(keys(T2xT)))

elapsed = map(all_Ts) do T
    elapsed_for_T = Pair{Int64, Tuple{Float64, Float64}}[]
    for num_knots in 0:100
        data = load_elapsed(num_knots, T)
        data === nothing && continue
        d = data.elapsed[1000:end]
        push!(elapsed_for_T, (num_knots => (mean(d), std(d))))
    end
    T => elapsed_for_T
end

@save "elapsed_time.jld2" elapsed

ess = map(all_Ts) do T
    println("------------")
    println("T: $T")
    println("------------")
    println()
    ess_for_T = Pair{Int64, Float64}[]
    for num_knots in 0:100
        data = load_midpt_data(num_knots, T)
        data === nothing && continue
        println("number of knots: $num_knots...")

        if num_knots % 4 == 0
            left_num_knots = div(num_knots, 2)
            right_num_half_knots = div(num_knots, 4)
            col_idx = left_num_knots + 3 + right_num_half_knots
        elseif num_knots % 4 == 1
            left_num_half_knots = div(num_knots-1, 4)+1
            col_idx = left_num_half_knots + 1
        elseif num_knots % 4 == 2
            left_num_knots = div(num_knots, 2)
            right_num_half_knots = div(num_knots-2, 4)+1
            col_idx = left_num_knots + 3 + right_num_half_knots
        else
            left_num_knots = div(num_knots+1, 2)
            right_num_half_knots = div(num_knots-3, 4)+1
            col_idx = left_num_knots + 3 + right_num_half_knots
        end

        _ess = effective_sample_size(data[:,col_idx])
        push!(ess_for_T, (num_knots=>_ess))
    end
    T => ess_for_T
end

@save "effective_sample_size.jld2" ess




T_num_knots = [
    [(0.4, 3), (0.5, 4), (1.0, 9), (2.0, 19), (4.0, 39)],
    [(0.5, 2), (1.0, 5), (2.0, 11), (4.0, 23)]
]

cost_rel = map(T_num_knots) do line
    map(line) do Tnk
        T, num_knots = Tnk
        _el = filter(x->x[1] == num_knots, filter(e->e[1] == T, elapsed)[1][2])[1][2][1]
        _ess = filter(x->x[1] == num_knots, filter(e->e[1] == T, ess)[1][2])[1][2]
        Tnk => (_el, _ess)
    end
end

@save "cost_relation.jld2" cost_rel


using JLD2, Plots
gr()
@load "elapsed_time.jld2" elapsed
markers = Dict(
	0.2 => :circle,
	0.4 => :rect,
	0.5 => :star5,
	1.0 => :x,
	2.0 => :diamond,
	4.0 => :dtriangle,
)

linestyles = Dict(
	0.2 => :solid,
	0.4 => :dash,
	0.5 => :dot,
	1.0 => :dashdot,
	2.0 => :solid,
	4.0 => :solid,
)

begin
	p = plot(yaxis=:log)
	for elapsed_for_T in elapsed
        T, el = elapsed_for_T
		num_knots = getindex.(el, 1)
		el_mean_and_std = getindex.(el, 2)
        el_mean = getindex.(el_mean_and_std, 1)
        
        ms = ( T == 2.0 ? :cross : (T == 4.0 ? :circle : :pixel) )

        additional = (
            T == 2.0 ?
            (markershape = :x,) :
            (
                T == 4.0 ?
                (markershape = :+,) :
                tuple()
            )
        )

		plot!(
            p, num_knots, el_mean,
            label="T=$(T), xT=$(T2xT[T])",
            linestyle=linestyles[T],
            color="black";
            additional...
            #markershape = :markers[T],
        )
    end
    xlabel!(p, "Number of knots")
    ylabel!(p, "Elapsed time (sec)")
end
p


@load "effective_sample_size.jld2" ess


begin
	p = plot(yaxis=:log)
	for (elapsed_for_T, ess_for_T) in zip(elapsed, ess)
        T, el = elapsed_for_T
        T2, es = ess_for_T
        @assert T == T2
        
		num_knots = getindex.(el, 1)
		el_mean_and_std = getindex.(el, 2)
        el_mean = getindex.(el_mean_and_std, 1)
        if num_knots[1] == 0
            num_knots = num_knots[2:end]
            el_mean_and_std = el_mean_and_std[2:end]
            el_mean = el_mean[2:end]
        end

        num_knots2 = getindex.(es, 1)
        _ess = getindex.(es, 2)

        @assert all(num_knots .== num_knots2)
        
        ms = ( T == 2.0 ? :cross : (T == 4.0 ? :circle : :pixel) )

        additional = (
            T == 2.0 ?
            (markershape = :x,) :
            (
                T == 4.0 ?
                (markershape = :+,) :
                tuple()
            )
        )

		plot!(
            p, T./(num_knots.+1), _ess./el_mean,
            label="T=$(T), xT=$(T2xT[T])",
            linestyle=linestyles[T],
            color="black";
            additional...
            #markershape = :markers[T],
        )
    end
    xlabel!(p, "T/(number of knots+1)")
    ylabel!(p, "time-adjusted ESS [ESS/(Elapsed time in sec)]")
end
p

@load "cost_relation.jld2" cost_rel

markers=[:circle, :star]
linestyles=[:solid, :dash]

begin
    p = plot(yaxis=:log, xaxis=:log)
    for (i, rel) in enumerate(cost_rel)
        xdata = getindex.(rel, 1)
        ydata = getindex.(rel, 2)
        xaxis = getindex.(xdata, 1)
        label = "T/(num knots+1)=$(round(xdata[1][1]/(xdata[1][2]+1), digits=3))"
        yaxis = getindex.(ydata, 1) ./ getindex.(ydata, 2) .* 100000
        plot!(p, xaxis, yaxis, label=label, markershape=markers[i], linestyle=linestyles[i], color=:black)
    end
    k = 0.05
    plot!(p, [0.4, 4.0], [k*0.4^3, k*4.0^3], label="line with slope 3", linestyle=:dashdot, color=:black)
    xlabel!("T")
    ylabel!("cost per sample [(Elapsed time in sec)/ESS]")
end

xticks!([0.4, 0.5, 1.0, 2.0, 4.0], ["0.4", "0.5", "1", "2", "4"])

p

ess[1]

