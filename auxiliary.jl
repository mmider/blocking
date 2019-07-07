function linInterpol(x0, xT, Ts)
    T = sum(Ts)
    Tsc = [t for t in Ts]
    Tsc = cumsum(Tsc)
    x0 .* (1 .-Tsc./T) .+ xT .* Tsc ./ T
end
