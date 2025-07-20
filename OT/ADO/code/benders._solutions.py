from gurobipy import *

# ---------------------
# Problem Data
# ---------------------
facilities = [0, 1, 2, 3, 4]
customers = [0, 1, 2]

f = [50, 60, 55, 70, 65]
u = [60, 50, 80, 70, 55]
d = [40, 30, 30]
c = [
    [4, 6, 9],
    [5, 4, 8],
    [3, 5, 7],
    [6, 3, 6],
    [7, 5, 4]
]

# ---------------------
# Set the y value here
# ---------------------
y_fixed = [0, 0, 0, 1, 0]
#y_fixed = [1, 0, 1, 0, 0]
#y_fixed = [0, 0, 1, 0, 1]

# ---------------------
# Subproblem
# ---------------------
sub = Model("Subproblem")
sub.setParam("OutputFlag", 0)

x = sub.addVars(facilities, customers, lb=0, name="x")

# Demand constraints
demand_constrs = sub.addConstrs(
    (quicksum(x[i, j] for i in facilities) == 1 for j in customers),
    name="Demand"
)

# Capacity constraints
capacity_constrs = sub.addConstrs(
    (quicksum(d[j] * x[i, j] for j in customers) <= u[i] * y_fixed[i] for i in facilities), 
    name="Cap"
)

# Objective
sub.setObjective(quicksum(c[i][j] * d[j] * x[i, j] for i in facilities for j in customers), GRB.MINIMIZE)

sub.optimize()

# ---------------------
# Feasibility Cut via Farkas Dual
# ---------------------
if sub.status == GRB.INFEASIBLE:
    print("Subproblem infeasible. Using FarkasDual to extract feasibility cut...")

    # Compute Farkas certificate
    sub.computeIIS()

    # Get Farkas duals
    varlambda = [demand_constrs[j].FarkasDual for j in customers]
    varmu = [capacity_constrs[i].FarkasDual for i in facilities]

    y_coeffs = [-u[i] * varmu[i] for i in facilities]
    const = -sum(varlambda)

    print("Feasibility cut:")
    print(f"  0 ≥ {const:.2f}")
    for i in facilities:
        print(f"    + {y_coeffs[i]:.2f} * y_{i+1}")

# ---------------------
# Optimality Cut via Dual Solution
# ---------------------
else:
    print("Subproblem feasible. You can extract an optimality cut from the dual.")

    # Get duals
    varlambda = [demand_constrs[j].Pi for j in customers]
    varmu = [0] * len(facilities)
    for i in capacity_constrs:
        varmu[i] = capacity_constrs[i].Pi  # only for open facilities

    print(f"Dual variables (lambda): {varlambda}")
    print(f"Dual variables (mu): {varmu}")

    y_coeffs = [u[i] * varmu[i] for i in facilities]
    const = sum(varlambda)

    print("Optimality cut:")
    print(f"phi ≥ {const:.2f}")
    for i in facilities:
        print(f"    + {y_coeffs[i]:.2f} * y_{i+1}")
