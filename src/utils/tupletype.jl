
function gen_nt_type(vars::Vector{Symbol}, dtype)
    NamedTuple{tuple(vars...),Tuple{[dtype for _ in vars]...}}
end

function f1(v::gen_nt_type([:a,:b],T)) where {T<:Number}
    @info "ok"
end

function f2(v::(@NamedTuple{a::T,b::T}))  where {T<:Number}
    @info "ok"
end