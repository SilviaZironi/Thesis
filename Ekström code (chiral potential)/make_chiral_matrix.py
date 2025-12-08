import numpy as np
import constants as const           # modulo del repo Ekström
from NN_potential import chiral_LO as potential
from potential import V

def main():
    #leggi la griglia di TRIMOD (in fm^-1)
    with open("p_grid.dat") as f:
        first = f.readline()
        nptot = int(first.split()[0])
    p_fm = np.loadtxt("p_grid.dat", skiprows=1)
    p_fm = np.asarray(p_fm, dtype=float).flatten()
    assert len(p_fm) == nptot

    print(f"Loaded grid with {nptot} points")

    # --- converto in MeV/c per il potenziale di Ekström ---
    try:
        hbarc = const.hbarc
    except AttributeError:
        hbarc = 197.3269804  # MeV*fm

    p_MeV = p_fm * hbarc

    #canale 3S1, T=0 (deuterone np)
    L  = 0
    LL = 0
    S  = 1
    J  = 1
    T  = 0
    Tz = 0

    Vmat = np.zeros((nptot, nptot))

    for i, pi in enumerate(p_MeV):
        for j, pj in enumerate(p_MeV):
            Vmat[i, j] = V(pi, pj, L, LL, S, J, T, Tz)

    # debug: guarda ordini di grandezza
    print("Vmat min, max =", Vmat.min(), Vmat.max())
    
    with open("v_chiral_matrix.dat", "w") as f:
        f.write(f"{nptot}\n")
        for i in range(nptot):
            for j in range(nptot):
                f.write(f"{Vmat[i,j]: .16E} ")
            f.write("\n")

    print("Saved v_chiral_matrix.dat")

if __name__ == "__main__":
    main()


