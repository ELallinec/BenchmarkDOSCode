using LinearAlgebra
println("Arguments received:")
for (i, arg) in enumerate(ARGS)
	println("Arg $i: $arg")
end

# You can also parse the arguments, e.g.:
x = parse(Int, ARGS[1])
y = parse(Float64, ARGS[2])
println("x + y = ", x + y)