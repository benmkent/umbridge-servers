DarcyModel = UMBridge.Model(
    name = "darcy",
    inputSizes = [50], #[config["n_basis_input"]],
    outputSizes = [50, 50, 50], #[config["n_basis_output"]],
    supportsGradient = false,
    evaluate = (y, config) -> evaluate_darcy(y, config),
)

AdvDiffModel = UMBridge.Model(
    name = "advdiff",
    inputSizes = [1,50, 50, 50], #[config["n_basis_input"]],
    outputSizes = [50], #[config["n_basis_output"]],
    supportsGradient = false,
    evaluate = (y, config) -> evaluate_advdiff(y, config),
)

"""
    evaluate_darcy(y, config)

Evaluates the Darcy model.

# Arguments
- `y`: Input data.
- `config`: Configuration dictionary.

# Returns
- Encoded pressure and velocity components.
"""
function evaluate_darcy(y, config)
    config = complete_config(config)
    y_trimmed = y[1][1:prod(config["n_basis_input"])]
    println("Evaluating with y: ", y_trimmed)
    println("Config: ", config)

    x_u, f_u, x_p, f_p = darcy_solve(y_trimmed, config)
    
    println("x_u: ", x_u)
    println("f_u: ", f_u)
    println("x_p: ", x_p)
    println("f_p: ", f_p)

    println("Encoding u,p")
    enc_u = [encoder(x_u, f_u_i, config["ndim"], config["n_basis_output"]) for f_u_i in eachcol(f_u) ]
    enc_p = encoder(x_p, f_p, config["ndim"], config["n_basis_output"])
    println("Encoded u: ", enc_u)
    println("Encoded p: ", enc_p)

    # Put extra zeros on end to make length 50
    enc_p_full = zeros(50)
    enc_ux_full = zeros(50)
    enc_uy_full = zeros(50)
    enc_p_full[1:length(enc_p)] = enc_p
    enc_ux_full[1:length(enc_u[1])] = enc_u[1]
    if length(enc_u) > 1
        enc_uy_full[1:length(enc_u[2])] = enc_u[2]
    end

    return enc_p_full, enc_ux_full, enc_uy_full
end

"""
    evaluate_advdiff(y, config)

Evaluates the advection-diffusion model.

# Arguments
- `y`: Input data.
- `config`: Configuration dictionary.

# Returns
- Encoded solution.
"""
function evaluate_advdiff(y, config)
    config = complete_config(config)
    diffusion_coefficient = y[1][1]
    p_trimmed = y[2][1:prod(config["n_basis_output"])]
    ux_trimmed = y[3][1:prod(config["n_basis_output"])]
    uy_trimmed = y[4][1:prod(config["n_basis_output"])]
    
    ux = decoder(ux_trimmed, config["ndim"])
    uy = decoder(uy_trimmed, config["ndim"])

    # Solve the PDE
    x, u = advdiff_solve(ux, uy, diffusion_coefficient, config)
    println("x: ", x)
    println("u: ", u)

    # Encode the solution
    enc_u = encoder(x, u, config["ndim"], config["n_basis_output"])
    println("Encoded u: ", enc_u)

    # Put extra zeros on end to make length 50
    enc_u_full = zeros(50)
    enc_u_full[1:length(enc_u)] = enc_u
    return enc_u_full
end

"""
    advdiff_solve(adv_u1, adv_u2, ε, config)

Solves the advection-diffusion PDE.

# Arguments
- `adv_u1`: Advection in x-direction.
- `adv_u2`: Advection in y-direction.
- `ε`: Diffusion coefficient.
- `config`: Configuration dictionary.

# Returns
- Solution grid and values.
"""
function advdiff_solve(adv_u1,adv_u2, ε, config)
    pyimport("sys").path.append(".")
    AdvDiffPDE = pyimport("AdvectionDiffusionPDE")

    py_ux = Py((x, y) -> begin
            x = PyArray(x)
            y = PyArray(y)
            ux_out = Py(adv_u1.(x, y))
            Py(ux_out)
        end)
    py_uy = Py((x, y) -> begin
        x = PyArray(x)
        y = PyArray(y)
        uy_out = Py(adv_u2.(x, y))
        Py(uy_out)
    end)
    advdiff = AdvDiffPDE.AdvDiffPDE()
    advdiff.setup_function_spaces(ndim=config["ndim"])
    advdiff.setup_problem(py_ux, py_uy, ε)
    advdiff.solve()
    x, u = advdiff.get_u()
    x_out = pyconvert(Array,x)
    u_out = pyconvert(Array,u)
    return x_out,u_out
end

"""
    darcy_solve(y, config)

Solves the Darcy PDE.

# Arguments
- `y`: Input data.
- `config`: Configuration dictionary.

# Returns
- Solution grid and values for velocity and pressure.
"""
function darcy_solve(y, config)
    # Import Python module
    pyimport("sys").path.append(".")
    DarcyPDE = pyimport("DarcyPDE")
    
    # Instantiate and call methods
    d = DarcyPDE.DarcyPDE()
    d.setup_function_spaces(ndim=config["ndim"])
    d.setup_permeability(y, config["basis_input"], config["n_basis_input"])
    d.setup_problem()
    d.solve()
    u = d.get_u()
    p = d.get_p()

    x_u = pyconvert(Array,u[0])
    f_u = pyconvert(Array,u[1])
    x_p = pyconvert(Array,p[0])
    f_p = pyconvert(Array,p[1])
    return x_u, f_u, x_p, f_p
end

"""
    encoder(x, u, ndim, m=20)

Encodes the solution using Chebyshev basis.

# Arguments
- `x`: Grid points.
- `u`: Solution values.
- `ndim`: Number of dimensions.
- `m`: Number of basis functions (default: 20).

# Returns
- Encoded coefficients.
"""
function encoder(x,u, ndim, m=20)
    S = Chebyshev(-1..1)^ndim;
    n = size(x,1);
    V = Array{Float64}(undef,n,m);
    if ndim == 1
        x_data = [xii[1] for xii = eachrow(x[:,1:ndim])]
    else
        x_data = [xii for xii = eachrow(x[:,1:ndim])]
    end
    for k = 1:m

        V[:,k] = Fun(S,[zeros(k-1);1]).(x_data)
    end

    f = Fun(S,V\u)
    return f.coefficients
end

"""
    decoder(coeffs, ndim)

Decodes coefficients into a function.

# Arguments
- `coeffs`: Coefficients to decode.
- `ndim`: Number of dimensions.

# Returns
- Decoded function.
"""
function decoder(coeffs, ndim)
    S = Chebyshev(-1..1)^ndim;
    f = Fun(S,coeffs)
end

"""
    complete_config(config)

Completes the configuration dictionary with defaults.

# Arguments
- `config`: Configuration dictionary.

# Returns
- Completed configuration dictionary.
"""
function complete_config(config)
    config = isempty(config) ? Dict() : config

    defaults = Dict(
        "ndim" => 2,
        "basis_input" => "patch",
        "n_basis_input" => fill(4, get(config, "ndim", 2)),
        "n_basis_output" => 20
    )

    for (key, value) in defaults
        config[key] = get(config, key, value)
    end
    
    return config
end

UMBridge.serve_models([DarcyModel, AdvDiffModel], 4242)
