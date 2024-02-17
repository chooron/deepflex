using QuadGK

integral, error = quadgk(x -> cos(200x), 0, 1)