# Author : Jacob M. Aguirre
# 23rd April, 2024

# This code solves the following cutting stock model:

# Master problem:
#     min     \sum_{p in P} x_p
#     s.t.    \sum_{p in P} patterns_{ip} * x_p ≥ d_i, for i in I
#             x_p ≥ 0 and integer, for p in P

# Subproblem:
#     min     1 - \sum_{i in I} price_i * use_i
#     s.t.    \sum_{i in I} w_i * use_i ≤ W_roll
#             use_i ≥ 0 and integer, for i in I

# =========================== Julia Code ======================================

using JuMP, Gurobi, Random, Logging, LinearAlgebra, MosekTools, Mosek

function keyboard_terminate(model, where) # Enable pause m.optimize() by 'ctrl + c'
    try
        pass
    catch
        model.terminate()
    end
end

# Set up logger
logger = logging.getLogger(__name__)
if isempty(logger.handlers)
    logger.addHandler(logging.StreamHandler())
    file_handler = logging.FileHandler("RunLog.log")
    formatter = logging.Formatter(
        fmt="%(asctime)s %(filename)s:%(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG)
end

logger.info("Begin")

# =========================== Parameters ======================================
TOL = 1e-6
W_roll = 100
I = 1:5
w = rand(1:50, length(I))
d = rand(1:50, length(I))
patterns = [Diagonal([W_roll // w[i] for i in I])...]

# ========================= Master Problem ====================================
m = Model(Gurobi.Optimizer)
set_optimizer_attribute(m, "OutputFlag", 0)
@variable(m, x[1:length(patterns)] >= 0)
@objective(m, Min, sum(x))
@constraint(m, c1[i in I], patterns[i][i] * x[i] >= d[i])

# ======================= Subproblem and Iteration ============================
for iter_count in count()
    optimize!(m)
    price = [dual(c1[i]) for i in I]
    println("Price = ", price)

    sp = Model(Gurobi.Optimizer)
    set_optimizer_attribute(sp, "OutputFlag", 0)
    @variable(sp, use[i in I] >= 0, Int)
    @objective(sp, Max, sum(price[i] * use[i] for i in I))
    @constraint(sp, sum(w[i] * use[i] for i in I) <= W_roll)
    optimize!(sp)
    min_rc = 1 - objective_value(sp)
    if min_rc < -TOL
        push!(patterns, [value(use[i]) for i in I])
        logger.info("min reduced cost = $(round(min_rc, digits=4)); new pattern: $(patterns[end])")
        @variable(m, x[iter_count + length(I)] >= 0)
        set_start_value(x[iter_count + length(I)], 0)
        set_start_value(x[iter_count + length(I)], patterns[end])
    else
        break
    end
end

# Relaxed Model Result
min_rc = 1 - objective_value(sp)
if min_rc >= 0
    logger.info("min reduced cost = $(round(min_rc, digits=4)) ≥ 0")
    relaxed_result = ["$(round(value(v), digits=4)) * $(patterns[p])" for (p, v) in enumerate(all_variables(m)) if value(v) > TOL]
    unshift!(relaxed_result, "Relaxed result = $(round(objective_value(m), digits=4)) rolls")
    logger.info(join(relaxed_result, "\n\t"))
end

# ====================== Integer Model Result =================================
# Integer Model Result
for v in all_variables(m)
    set_integer(v)
end
optimize!(m)
integer_result = ["$(round(value(v), digits=4)) * $(patterns[p])" for (p, v) in enumerate(all_variables(m)) if value(v) > TOL]
unshift!(integer_result, "Integer result = $(round(objective_value(m), digits=4)) rolls")
logger.info(join(integer_result, "\n\t"))

logger.info("End")
"""


