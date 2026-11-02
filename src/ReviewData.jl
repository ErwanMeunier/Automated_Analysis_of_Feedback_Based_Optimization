"""Solver-free loading, correction overlays, and auditing of reviewer tables."""
module ReviewData
using Arrow, DataFrames
export ROOT, load_reference, paper_domain, point_key, audit_reference, overlay

const ROOT = normpath(joinpath(@__DIR__, ".."))
point_key(row) = (Int(row.T), round(Float64(row.norm_A); digits=10),
                  round(Float64(row.kappa); digits=10))
paper_domain(df) = df[(0 .<= df.norm_A .<= 0.9 + 1e-10) .&
                      in.(df.kappa, Ref((0.0, 0.1))), :]

function overlay(original, corrections)
    result = copy(original)
    index = Dict(point_key(row) => i for (i, row) in enumerate(eachrow(result)))
    for row in eachrow(corrections)
        # A failed retry never replaces an existing reference result.
        String(row.status) == "OK" || continue
        key = point_key(row)
        if haskey(index, key)
            for field in (:n_steps,:mu_f,:L_f,:norm_B,:R_x,:R_u)
                hasproperty(result,field) && hasproperty(row,field) || continue
                isapprox(result[index[key],field],row[field];atol=1e-12,rtol=1e-12) ||
                    error("correction changes physical parameters for $key: $field")
            end
            for field in propertynames(result)
                result[index[key], field] = row[field]
            end
        else
            push!(result, NamedTuple(row); cols=:setequal)
            index[key] = nrow(result)
        end
    end
    return sort!(result, [:T, :norm_A, :kappa])
end

function load_reference(; directory=joinpath(ROOT, "data", "reference"),
                          corrections=joinpath(ROOT, "data", "corrections"))
    tables = Dict{Symbol,DataFrame}()
    for (kind, file) in ((:sweep, "merged_sweep.arrow"),
                         (:design, "merged_design.arrow"),
                         (:analytic, "merged_analytic_theorem1.arrow"))
        path = joinpath(directory, file)
        isfile(path) || error("missing reference table: $path")
        table = DataFrame(Arrow.Table(path); copycols=true)
        required = [:T, :norm_A, :kappa, :mu_f, :L_f, :norm_B, :R_x, :R_u, :status]
        all(in(propertynames(table)), required) || error("invalid schema: $file")
        all(isapprox.(table.kappa, table.mu_f ./ table.L_f; atol=1e-12)) ||
            error("inconsistent kappa=mu_f/L_f in $file")
        correction = corrections === nothing ? "" : joinpath(corrections, file)
        kind != :sweep && isfile(correction) &&
            (table = overlay(table, DataFrame(Arrow.Table(correction))))
        tables[kind] = paper_domain(table)
    end
    return (sweep=tables[:sweep], design=tables[:design], analytic=tables[:analytic])
end

"""One row per expected stable-domain (T, norm_A, kappa) design group."""
function audit_reference(tables)
    design = Dict(point_key(row) => row for row in eachrow(tables.design))
    analytic = Dict(point_key(row) => row for row in eachrow(tables.analytic))
    keys_all = sort!(collect(union(keys(design), keys(analytic))))
    rows = map(keys_all) do key
        d, a = get(design, key, nothing), get(analytic, key, nothing)
        d_ok = d !== nothing && String(d.status) == "OK" &&
               isfinite(d.alpha_best) && isfinite(d.pep_best)
        a_ok = a !== nothing && String(a.status) == "OK" && isfinite(a.pep_analytic)
        (T=key[1], norm_A=key[2], kappa=key[3],
         design_status=d === nothing ? "MISSING" : String(d.status),
         analytic_status=a === nothing ? "MISSING" : String(a.status),
         design_ok=d_ok, analytic_ok=a_ok, needs_repair=!(d_ok && a_ok))
    end
    return DataFrame(rows)
end
end
