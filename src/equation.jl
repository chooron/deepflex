struct HydroEquation{input_names,output_names,param_names} <: AbstractEquation
    inputs::Vector{Num}
    outputs::Vector{Num}
    params::Vector{Num}

    function HydroEquation(
        input_names::Vector{Symbol},
        output_names::Vector{Symbol},
        param_names::Vector{Symbol},
    )
        inputs = vcat([@variables $var(t) = 0.0 for var in input_names]...)
        outputs = vcat([@variables $var(t) = 0.0 for var in output_names]...)
        params = vcat([@parameters $p = 0.0 [tunable = true] for p in param_names]...)
        return new{Tuple(input_names),Tuple(output_names),Tuple(param_names)}(inputs, outputs, params)
    end

    function HydroEquation(
        inputs::Vector{Num},
        outputs::Vector{Num},
        params::Vector{Num},
    )
        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        param_names = Symbolics.tosymbol.(params, escape=false)
        return new{Tuple(input_names),Tuple(output_names),Tuple(param_names)}(inputs, outputs, params)
    end
end
