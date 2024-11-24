from fenics import *
import matplotlib.pyplot as plt

# Step 1: Define the geometry and mesh
length, width = 4.0, 1.0  # Length and width of the artery
mesh = RectangleMesh(Point(0, 0), Point(length, width), 40, 10)

# Step 2: Define the function spaces
V = VectorFunctionSpace(mesh, 'P', 2)  # Velocity space
Q = FunctionSpace(mesh, 'P', 1)        # Pressure space
#W = MixedFunctionSpace([V, Q])         # Combined space for velocity and pressure
W = FunctionSpace([V, Q])   

# Step 3: Define the problem variables
u, p = TrialFunctions(W)  # Trial functions (for solving)
v, q = TestFunctions(W)   # Test functions (for testing)
w = Function(W)           # Solution function

# Fluid properties (blood parameters)
rho = 1060.0      # Density (kg/m^3)
mu = 0.004        # Dynamic viscosity (Pa·s)

# Step 4: Define boundary conditions
inflow = Expression(('4.0 * (1.0 - pow(x[1]/width, 2))', '0.0'), width=width, degree=2)  # Parabolic inflow
noslip = Constant((0, 0))  # No-slip condition
bc1 = DirichletBC(W.sub(0), inflow, 'near(x[0], 0)')  # Inflow
bc2 = DirichletBC(W.sub(0), noslip, 'on_boundary && !(near(x[0], 0) || near(x[0], length))')  # Walls
bc3 = DirichletBC(W.sub(1), Constant(0), 'near(x[0], length)')  # Zero pressure at the outlet
bcs = [bc1, bc2, bc3]

# Step 5: Define the Navier-Stokes equations (weak form)
F = (
    rho * dot((u - v) / Constant(1.0), v) * dx
    + mu * inner(grad(u), grad(v)) * dx
    - div(v) * p * dx
    + q * div(u) * dx
)

# Step 6: Solve the problem
solve(lhs(F) == rhs(F), w, bcs)
u, p = w.split()  # Extract velocity and pressure solutions

# Step 7: Visualize the results
plt.figure()
plot(u, title="Velocity Field")
plt.figure()
plot(p, title="Pressure Field")
plt.show()
