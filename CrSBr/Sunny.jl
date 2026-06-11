using Sunny, GLMakie, LinearAlgebra

units = Units(:meV, :angstrom) # Ensures correct units
cif = "CrSBr_mp-22998_symmetrized.cif" # Pulls the cif file from the directory
cryst = Crystal(cif; symprec=0.001) # Imports the crystal with symprec being a discrepancy for symmetry
cryst_cr = subcrystal(cryst, "Cr0") # Isolates the chromium atoms
save(joinpath(@__DIR__, "Crystal.png"), view_crystal(cryst; ndims=3, ghost_radius=5)) # Saves the crystal as a png

momentinfo = [1 => Moment(s=3/2, g=2)] # Initiates the moment of the chromium atoms
sys = System(cryst_cr, momentinfo, :dipole); # Initializes the system

# All the exchange constants and their bonds they are exchanging by. Sunny ensures all symmetry equivalent bonds are exchanged as well
J1 = -1.97 * I + dmvec([0, 0.4, 0])
J2 = -3.38
J3 = -1.67 * I + dmvec([0.4, 0, 0])
J4 = -0.09
J5 = J4
J7 = 0.37
J8 = -0.29
set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, J2, Bond(1, 2, [0, 0, -1]))
set_exchange!(sys, J3, Bond(1, 1, [0, 1, 0]))
set_exchange!(sys, J4, Bond(1, 1, [1, 1, 0]))
set_exchange!(sys, J5, Bond(1, 2, [1, 0, -1]))
set_exchange!(sys, J7, Bond(1, 1, [2, 0, 0]))
set_exchange!(sys, J8, Bond(2, 1, [0, 1, 1]))

D = 0.144*10^(-3)
set_onsite_coupling!(sys, S -> -D*S[2]^2, 1) # Sets the single ion anisotropy with the easy axis along y

randomize_spins!(sys) # Randomizes the spins
minimize_energy!(sys) # Finds the ground state
energy_per_site(sys)

sys_super = repeat_periodically(sys, (6, 6, 1)) # Repeats the system to make a supercell for better visualization

save(joinpath(@__DIR__, "Spins.png"), plot_spins(sys_super; ndims=3, ghost_radius=3)) # Saves a png of the spin configuration


function save_int(sys, h1,h2,k1,k2,E1,E2,Eres,name; fpath1=name, fpath2="Energies.txt",size=[500,1000], Elim=[E1,E2])
    # Save .txt file with intensities of spinwave colorplot
    formfactors = [1 => FormFactor("Cr0")]
    swt = SpinWaveTheory(sys; measure=ssf_perp(sys; formfactors))

    l_values = range(-1, 1, length=1000)
    acc = nothing

    for l in l_values
        path = q_space_path(cryst_cr, [[h1, k1, l], [h2, k2, l]], size[1])
        res = intensities(
            swt,
            path;
            energies=range(Elim[1], Elim[2], size[2]),
            kernel=gaussian(fwhm=Eres)
        )

        if acc === nothing
            acc = zeros(Base.size(res.data))
        end

        acc .+= res.data
    end

    res_data = acc ./ length(l_values)
    res_energies = collect(range(Elim[1], Elim[2], size[2]))

    open(fpath1, "w") do f
        for row in eachrow(res_data)
            println(f, join(row, ", "))
        end
    end
    open(fpath2, "w") do f
        println(f, join(res_energies, ", "))
    end
end

function save_pol_int_yy(sys,h1,h2,k1,k2,E1,E2,Eres,name; fpath1=name, fpath2="Energies.txt",size=[500,1000], Elim=[E1,E2])
    formfactors = [1 => FormFactor("Cr0")]

    measure = ssf_custom_bm(sys; formfactors=formfactors, u=[1, 0, 0], v=[0, 1, 0]) do q, ssf
        real(ssf[2,2])
    end
    swt = SpinWaveTheory(sys; measure)

    path = q_space_path(cryst_cr, [[h1, k1, 0], [h2, k2, 0]], size[1])
    res = intensities(
        swt,
        path;
        energies=range(Elim[1], Elim[2], size[2]),
        kernel=gaussian(fwhm=Eres)
    )
    res_data = res.data
    res_energies = collect(range(Elim[1], Elim[2], size[2]))

    open(fpath1, "w") do f
        for row in eachrow(res_data)
            println(f, join(row, ", "))
        end
    end
    open(fpath2, "w") do f
        println(f, join(res_energies, ", "))
    end
end

function save_pol_int_zz(sys,h1,h2,k1,k2,E1,E2,Eres,name; fpath1=name, fpath2="Energies.txt",size=[500,1000], Elim=[E1,E2])
    formfactors = [1 => FormFactor("Cr0")]

    measure = ssf_custom_bm(sys; formfactors=formfactors, u=[1, 0, 0], v=[0, 1, 0]) do q, ssf
        real(ssf[3,3])
    end
    swt = SpinWaveTheory(sys; measure)

    path = q_space_path(cryst_cr, [[h1, k1, 0], [h2, k2, 0]], size[1])
    res = intensities(
        swt,
        path;
        energies=range(Elim[1], Elim[2], size[2]),
        kernel=gaussian(fwhm=Eres)
    )
    res_data = res.data
    res_energies = collect(range(Elim[1], Elim[2], size[2]))

    open(fpath1, "w") do f
        for row in eachrow(res_data)
            println(f, join(row, ", "))
        end
    end
    open(fpath2, "w") do f
        println(f, join(res_energies, ", "))
    end
end

# function save_int_bands(sys, h1,h2,k1,k2,E1,E2,name; fpath=name)
#     # Save .txt file with intensities of spinwave colorplot
#     formfactors = [1 => FormFactor("Cr0")]
#     swt = SpinWaveTheory(sys; measure=ssf_perp(sys; formfactors))

#     path = q_space_path(cryst_cr, [[h1, k1, 1], [h2, k2, 1]], 500)
#     res = intensities_bands(swt, path)

#     save(joinpath(@__DIR__, "Bands.png"), plot_intensities(res; units=units))

#     open(fpath, "w") do f
#         for row in eachrow(res.data)
#             println(f, join(row, ", "))
#         end
#     end
#end