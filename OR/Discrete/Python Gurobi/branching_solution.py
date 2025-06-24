import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import itertools
from gurobipy import Model, GRB, quicksum
import time
import os

np.random.seed(42)

# ---------------------
# Parameters
# ---------------------
m = 100  # Universe size (number of target audiences)
n = 150  # Number of subsets (number of celebrities)

subsets = []
costs = []

# ---------------------
# Random instance generation
# ---------------------
for _ in range(n):
    center = np.random.randint(0, m)
    num_covered = np.random.choice([3, 4, 5, 10, 20], p=[0.2, 0.3, 0.3, 0.15, 0.05])
    coverage = set([center])
    while len(coverage) < num_covered:
        coverage.add(np.random.randint(0, m))
    subsets.append(coverage)
    cost = 1.5 * len(coverage) + np.random.normal(0, 2)
    costs.append(max(1, round(cost, 2)))

# Ensure full coverage
all_covered = set().union(*subsets)
for elem in set(range(m)) - all_covered:
    subsets.append({elem})
    costs.append(1.0)

n = len(subsets)

# ---------------------
# Auxiliary Functions
# ---------------------
class BBTreeNode:
    def __init__(self, fixed_vars, depth, parent_id=None, branch=None):
        self.fixed_vars = fixed_vars
        self.depth = depth
        self.node_id = next(node_counter)
        self.parent_id = parent_id
        self.branch = branch  # (var, value)

node_counter = itertools.count()

# Visualization helper
def hierarchy_pos(G, root=None, width=1.0, vert_gap=0.3, vert_loc=0, xcenter=0.5):
    if root is None:
        root = list(nx.topological_sort(G))[0]
    def _hierarchy_pos(G, root, leftmost, width, vert_gap, vert_loc, pos, parent=None):
        children = list(G.successors(root))
        if not children:
            pos[root] = (leftmost[0], vert_loc)
            leftmost[0] += width
        else:
            start = leftmost[0]
            for child in children:
                _hierarchy_pos(G, child, leftmost, width / len(children), vert_gap, vert_loc - vert_gap, pos, root)
            mid = (start + leftmost[0] - width / len(children)) / 2
            pos[root] = (mid, vert_loc)
        return pos
    return _hierarchy_pos(G, root, [0], width, vert_gap, vert_loc, {})

# Draw the branch-and-bound tree
def draw_branch_and_bound_tree(tree, title="Branch-and-Bound Tree"):
    pos = hierarchy_pos(tree)
    plt.figure(figsize=(18, 10))
    nx.draw(tree, pos, with_labels=False, node_size=800, node_color='skyblue', edge_color='gray')
    plt.title(title)
    plt.axis('off')
    plt.margins(0.1)
    plt.subplots_adjust(left=0, right=1, top=0.95, bottom=0)
    # Create a safe filename by replacing spaces and parentheses
    filename = title.replace(" ", "_").replace("(", "").replace(")", "") + ".png"
    # Get directory where this script resides
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, filename)
    plt.savefig(filepath, bbox_inches='tight')
    plt.close()
    print(f"Saved branch-and-bound tree image: {filepath}")

# This function might be helpful...
def get_uncovered_elements(fixed_vars):
    covered = set()
    for j, val in fixed_vars.items():
        if val == 1:
            covered.update(subsets[j])
    return set(range(m)) - covered

# ---------------------
# Implement LP relaxation (with support for fixed variables)
# ---------------------
def build_lp_model_scp(fixed_vars):
    model = Model()
    model.setParam('OutputFlag', 0)
    x = {j: model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, obj=costs[j], name=f"x_{j}")
         for j in range(n)}
    for j, val in fixed_vars.items():
        x[j].LB = val
        x[j].UB = val
    for i in range(m):
        model.addConstr(quicksum(x[j] for j in range(n) if i in subsets[j]) >= 1)
    model.modelSense = GRB.MINIMIZE
    model.update()
    return model, x

# ---------------------
# Implement custom branch-and-bound algorithm
# ---------------------
def manual_branch_and_bound_scp_with_tree(branching_rule="max_coverage"):
    tree = nx.DiGraph()
    best_obj = float('inf')
    best_sol = None
    explored_nodes = 0

    root = BBTreeNode(fixed_vars={}, depth=0)
    stack = [root]
    tree.add_node(root.node_id, label="root")

    while stack:
        node = stack.pop()
        explored_nodes += 1

        model, x = build_lp_model_scp(node.fixed_vars)
        model.optimize()

        if model.Status != GRB.OPTIMAL:
            tree.nodes[node.node_id]['label'] = "INFEASIBLE"
            continue

        obj_val = model.ObjVal
        if obj_val >= best_obj:
            tree.nodes[node.node_id]['label'] = f"Pruned (obj={obj_val:.1f})"
            continue

        sol = model.getAttr('X', x)
        fractional = [(j, val) for j, val in sol.items() if 1e-6 < val < 1 - 1e-6]

        if not fractional:
            if obj_val < best_obj:
                best_obj = obj_val
                best_sol = [j for j, val in sol.items() if val > 0.5]
            tree.nodes[node.node_id]['label'] = f"Leaf (obj={obj_val:.1f})"
            continue

        if branching_rule == "naive":
            branch_var = fractional[0][0]
        elif branching_rule == "most_infeasible":
            branch_var = min(fractional, key=lambda t: abs(t[1] - 0.5))[0]
        elif branching_rule == "max_coverage":
            uncovered = get_uncovered_elements(node.fixed_vars)
            branch_var = max(fractional, key=lambda t: len(subsets[t[0]] & uncovered))[0]
        else:
            raise ValueError("Unknown rule")

        for val in [1, 0]:
            child_vars = node.fixed_vars.copy()
            child_vars[branch_var] = val
            child = BBTreeNode(child_vars, node.depth + 1,
                               parent_id=node.node_id,
                               branch=(branch_var, val))
            stack.append(child)
            label = f"x[{branch_var}]={val}"
            tree.add_node(child.node_id, label=label)
            tree.add_edge(node.node_id, child.node_id)

    return best_obj, best_sol, explored_nodes, tree

# Run and compare strategies
if __name__ == "__main__":
    for rule in ["naive", "most_infeasible", "max_coverage"]:
        print(f"\nRunning with branching rule: {rule}")
        start = time.time()
        best_obj, best_sol, nodes, tree = manual_branch_and_bound_scp_with_tree(rule)
        duration = time.time() - start
        print(f"Best objective: {best_obj:.2f}")
        print(f"Solution sets used: {len(best_sol)}")
        print(f"Nodes explored: {nodes}")
        print(f"Time taken: {duration:.2f}s")
        draw_branch_and_bound_tree(tree, title=f"Branch-and-Bound Tree ({rule})")
