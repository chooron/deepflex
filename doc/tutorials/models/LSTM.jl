using Lux
using Lux: AbstractExplicitContainerLayer
using LuxCore
using StableRNGs

struct SpiralClassifier{L,C} <: AbstractExplicitContainerLayer{(:lstm_cell, :classifier)}
    lstm_cell::L
    classifier::C
end

function SpiralClassifier(in_dims, hidden_dims, out_dims)
    return SpiralClassifier(
        LSTMCell(in_dims => hidden_dims), Dense(hidden_dims => out_dims, sigmoid))
end

function (s::SpiralClassifier)(
    x::AbstractArray{T,3}, ps::NamedTuple, st::NamedTuple) where {T}
    # First we will have to run the sequence through the LSTM Cell
    # The first call to LSTM Cell will create the initial hidden state
    # See that the parameters and states are automatically populated into a field called
    # `lstm_cell` We use `eachslice` to get the elements in the sequence without copying,
    # and `Iterators.peel` to split out the first element for LSTM initialization.
    x_init, x_rest = Iterators.peel(LuxOps.eachslice(x, Val(2)))
    (y, carry), st_lstm = s.lstm_cell(x_init, ps.lstm_cell, st.lstm_cell)
    # Now that we have the hidden state and memory in `carry` we will pass the input and
    # `carry` jointly
    for x in x_rest
        (y, carry), st_lstm = s.lstm_cell((x, carry), ps.lstm_cell, st_lstm)
    end
    # After running through the sequence we will pass the output through the classifier
    y, st_classifier = s.classifier(y, ps.classifier, st.classifier)
    # Finally remember to create the updated state
    st = merge(st, (classifier=st_classifier, lstm_cell=st_lstm))
    return vec(y), st
end

classify = SpiralClassifier(3, 10, 1)
classify_func = (x, ps, st) -> LuxCore.apply(classify, x, ps, st)

classify_func(ones(3, 10, 1),)

lce = Lux.LSTMCell(3 => 10)
lps, lst = Lux.setup(StableRNGs.LehmerRNG(1234), lce)

lce_func = (x, ps, st) -> LuxCore.apply(lce, x, ps, st)[1]

re = lce_func(ones(3, 10), lps, lst)
length(re[2])