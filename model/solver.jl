using Flux, JLD2

cd(@__DIR__)
@load "model.jld2" nn

neural_collision(x) = nn(x)
