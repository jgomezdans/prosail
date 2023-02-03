import numpy as np
import matplotlib.pyplot as plt
import prosail

if __name__ == "__main__":
    # Th = 37; Tc=29;Ts=50;Td=26;
    lam = 9.5  # um
    tveg = 33.0  # degC
    tsoil = 42.0  # degC
    t_atm = -14.0  # degC
    emv = 0.98  #
    ems = 0.94  #

    lidfa = -1.0
    tveg_sunlit = 42.0
    tsoil_sunlit = 58.0
    hspot = 0.05
    tts = 24.2
    for lai in [0.5, 1, 2, 4]:
        BT = []
        angs = list(range(-90, 90, 1))
        for tto in angs:
            if tto <= 0:
                psi = 0.0
            else:
                psi = 180.0
            retval = prosail.run_thermal_sail(
                lam,
                tveg + 273.15,
                tsoil + 273.15,
                tveg_sunlit + 273.15,
                tsoil_sunlit + 273.15,
                t_atm + 273.15,
                lai,
                lidfa,
                hspot,
                tts,
                np.abs(tto),
                psi,
                emv=emv,
                ems=ems,
                lidfb=0.0,
                typelidf=1,
            )
            BT.append(retval[1] - 273.15)
        plt.plot(angs, BT, "-", label=f"LAI={lai}")
