using Sunny, GLMakie # Imports the Sunny Ppackage and GLMakie for plotting

units = Units(:meV, :angstrom) # Ensures units of meV and angstrom

# Arbitrary lattice vectors and positions
# Creates the crystal and saves it as a png
latvecs = lattice_vectors(3, 3, 3, 90, 90, 90) 
cryst = Crystal(latvecs, [[0, 0, 0]])
save(joinpath(@__DIR__, "AFM Heisenberg chain crystal.png"), view_crystal(cryst; ndims=3, ghost_radius=3))

# Intiates the system using the crystal and defining the moment and g-factor
# :dipole means we are treating the moment as a dipole
# dims=(2, 1, 1) means we are treating the system as a 2D system with 1 site in the y and z direction, so it is a chain in the x direction
sys = System(cryst, [1 => Moment(s=3/2, g=2)], :dipole; dims=(2, 1, 1))

J = 1 # Sets the exchange constant to be positive which means it is antiferromagnetic
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0])) # Sets the exchange of all relevant symmetries
randomize_spins!(sys) # Randomizes the spins to start the energy minimization from a non-trivial state
minimize_energy!(sys) # Minimizes the energy of the system to find the ground state configuration
energy_per_site(sys) # Prints the energy per site of the system
save(joinpath(@__DIR__, "AFM Heisenberg chain spins.png"),plot_spins(sys; ndims=3, ghost_radius=3)) # Saves a png of the spin configuration

# Calculates the spin wave spectrum and saves it as a png
# Saves the data as a .txt file
function save_int(sys; fpath=joinpath(@__DIR__, "AFM Heisenberg chain intensities.txt"),size=[500,1000], Elim=[0,5], Eres=0.3, qs=[[0,0,0],[1,0,0]])
    # Save .txt file with intensities of spinwave colorplot
    swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
    path = q_space_path(cryst, qs, size[1])

    res = intensities(swt, path; energies=range(Elim[1],Elim[2],size[2]), kernel=gaussian(fwhm=Eres))
    save(joinpath(@__DIR__, "AFM Heisenberg chain intensities.png"), plot_intensities(res; units=units))

    open(fpath, "w") do f
        for row in eachrow(res.data)
            println(f, join(row, ", "))
        end
    end
end
save_int(sys)

# Calculates the spin wave spectrum and saves the energies as a .txt file
function save_int_energies(sys; fpath=joinpath(@__DIR__, "AFM Heisenberg chain energies.txt"),size=[500,1000], Elim=[0,5], Eres=0.3, qs=[[0,0,0],[1,0,0]])
    # Save .txt file with intensities of spinwave colorplot
    swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
    path = q_space_path(cryst, qs, size[1])

    res = intensities(swt, path; energies=range(Elim[1],Elim[2],size[2]), kernel=gaussian(fwhm=Eres))

    open(fpath, "w") do f
        for row in eachrow(res.energies)
            println(f, join(row, ", "))
        end
    end
end
save_int_energies(sys)