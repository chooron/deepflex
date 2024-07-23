function Base.show(io::IO, flux::AbstractSimpleFlux)
    compact = get(io, :compact, false)

    if compact
        for (flux_expr, output) in zip(flux.flux_exprs, collect(flux.output_info))
            println(io, output ~ flux_expr)
        end
    else
        for (flux_expr, output) in zip(flux.flux_exprs, collect(flux.output_info))
            println(io, output ~ flux_expr)
        end
    end
end

function Base.show(io::IO, flux::AbstractStateFlux)
    compact = get(io, :compact, false)

    if compact
        println(io, first(collect(flux.output_info)) ~ flux.state_expr)
    else
        println(io, first(collect(flux.output_info)) ~ flux.state_expr)
    end
end

function Base.show(io::IO, flux::AbstractLagFlux)
    compact = get(io, :compact, false)
    # todo 创建用于展示lag flux的函数
    if compact
        println(io, first(collect(flux.output_info)) ~ flux.state_expr)
    else
        println(io, first(collect(flux.output_info)) ~ flux.state_expr)
    end
end

function Base.show(io::IO, flux::AbstractNeuralFlux)
    compact = get(io, :compact, false)
    # todo 创建用于展示nn flux的函数
    if compact
        println(io, first(collect(flux.output_info)) ~ flux.state_expr)
    else
        println(io, first(collect(flux.output_info)) ~ flux.state_expr)
    end
end