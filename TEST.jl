using Pardiso
MKLPardisoSolver()
get_matrixtype(ps)
set_nprocs!(ps, 5) # Sets the number of threads to use
get_nprocs(ps) # Gets the number of threads being used
A = sparse(rand(10, 10))
B = rand(10, 2)
X = zeros(10, 2)
solve!(ps, X, A, B)
get_matrixtype(ps)