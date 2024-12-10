from fenics import *

# Create mesh and define function space
nx = 32
mesh = UnitIntervalMesh(nx)
V = FunctionSpace(mesh, "P", 1)

# Define boundary condition
u_D = Expression("1 + x[0]*x[0]", degree=2)
bc = DirichletBC(V, u_D, "on_boundary")

# Define problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v)) * dx
L = f * v * dx

# Solve
u = Function(V)
solve(a == L, u, bc)

# Plot solution
import matplotlib.pyplot as plt
plot(u)
plt.title("Solution of Poisson Equation (FEM)")
plt.show()