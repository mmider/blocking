### A Pluto.jl notebook ###
# v0.11.9

using Markdown
using InteractiveUtils

# ╔═╡ 3357e5dc-eaae-11ea-0ad3-31ea5622f7cd
using JLD2, Plots, LaTeXStrings

# ╔═╡ 49db2846-eaae-11ea-060b-0b7ea3c77fb7
begin
	@load "elapsed_time.jld2" elapsed
	@load "effective_sample_size.jld2" ess
	@load "cost_relation.jld2" cost_rel
end

# ╔═╡ ea3730b6-eab1-11ea-1a4e-6f5b44515319
T2xT = Dict(
    0.2 => 0.1,
    0.4 => 0.85,
    0.5 => 0.85,
    1.0 => 0.95,
    2.0 => 2.5,
    4.0 => 4.85,
)

# ╔═╡ 5a9c766c-eaae-11ea-131a-2fc8688c4ff1
begin
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
end

# ╔═╡ fc9233c0-eaaf-11ea-1080-ddae46c87a13
begin
	p1 = plot(yaxis=:log, size=(700, 350))
	for elapsed_for_T in elapsed
        T, el = elapsed_for_T
		num_knots = getindex.(el, 1)
		el_mean_and_std = getindex.(el, 2)
        el_mean = getindex.(el_mean_and_std, 1)
        
        ms = ( T == 2.0 ? :cross : (T == 4.0 ? :circle : :pixel) )

        additional = (
            T == 2.0 ?
            (markershape = :x, markersize=3) :
            (
                T == 4.0 ?
                (markershape = :+, markersize=3) :
                tuple()
            )
        )

		plot!(
            p1, num_knots, el_mean,
            label=latexstring("\$T=$(T), x_T=$(T2xT[T])\$"),
            linestyle=linestyles[T],
            color="black",
			linewidth=2.0;
            additional...
            #markershape = :markers[T],
        )
    end
    xlabel!(p1, "Number of knots")
    ylabel!(p1, "Elapsed time (sec)")
end

# ╔═╡ e1d1ee2a-eaaf-11ea-3685-2bad1de12009
begin
	p2 = plot(yaxis=:log, size=(700, 350))
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
            (markershape = :x, markersize=3) :
            (
                T == 4.0 ?
                (markershape = :+, markersize=3) :
                tuple()
            )
        )

		plot!(
            p2, T./(num_knots.+1), _ess./el_mean,
            label=latexstring("\$T=$(T), x_T=$(T2xT[T])\$"),
            linestyle=linestyles[T],
            color="black",
			linewidth=2.0;
            additional...
            #markershape = :markers[T],
        )
    end
    xlabel!(p2, "T/(number of knots+1)")
    ylabel!(p2, "time-adjusted ESS [ESS/(Elapsed time in sec)]")
end

# ╔═╡ 4bd01f4a-eab2-11ea-257f-736e941822f8
begin
	markers3=[:circle, :star]
	linestyles3=[:solid, :dash]
end

# ╔═╡ df4a82ae-eaaf-11ea-0833-755423cd9474
begin
    p3 = plot(yaxis=:log, xaxis=:log, legend=:topleft, size=(700, 350))
    for (i, rel) in enumerate(cost_rel)
        xdata = getindex.(rel, 1)
        ydata = getindex.(rel, 2)
        xaxis = getindex.(xdata, 1)
        label = "T/(num knots+1)=$(round(xdata[1][1]/(xdata[1][2]+1), digits=3))"
        yaxis = getindex.(ydata, 1) ./ getindex.(ydata, 2) .* 100000
        plot!(
			p3, xaxis, yaxis,
			label=label,
			markershape=markers3[i],
			linestyle=linestyles3[i],
			color=:black
		)
    end
    k = 0.05
    plot!(
		p3, [0.4, 4.0], [k*0.4^3, k*4.0^3],
		label="line with slope 3",
		linestyle=:dashdot,
		color=:black
	)
    xlabel!("T")
    ylabel!("cost per sample [(Elapsed time in sec)/ESS]")
	xticks!([0.4, 0.5, 1.0, 2.0, 4.0], ["0.4", "0.5", "1", "2", "4"])
end

# ╔═╡ d2d52f60-eaaf-11ea-2eda-c5a974a630bf
savefig(p1, "blocking_sin_elapsed_bw.pdf")

# ╔═╡ ba4b5208-eaaf-11ea-22c7-1dc8ea7bb333
savefig(p2, "blocking_sin_ta_ess_bw.pdf")

# ╔═╡ 9f56de68-eaaf-11ea-3812-b1610bc50568
savefig(p3, "blocking_sin_cost_pred_bw.pdf")

# ╔═╡ Cell order:
# ╠═3357e5dc-eaae-11ea-0ad3-31ea5622f7cd
# ╠═49db2846-eaae-11ea-060b-0b7ea3c77fb7
# ╠═ea3730b6-eab1-11ea-1a4e-6f5b44515319
# ╠═5a9c766c-eaae-11ea-131a-2fc8688c4ff1
# ╠═fc9233c0-eaaf-11ea-1080-ddae46c87a13
# ╠═e1d1ee2a-eaaf-11ea-3685-2bad1de12009
# ╠═4bd01f4a-eab2-11ea-257f-736e941822f8
# ╠═df4a82ae-eaaf-11ea-0833-755423cd9474
# ╠═d2d52f60-eaaf-11ea-2eda-c5a974a630bf
# ╠═ba4b5208-eaaf-11ea-22c7-1dc8ea7bb333
# ╠═9f56de68-eaaf-11ea-3812-b1610bc50568
