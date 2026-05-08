import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from PA3Py.PebbleAccretion3 import PebbleAccretionModule3

DATADIR = r"C:\astro\Codigos practica + docs + papers\codigos\data_gaps_pipeline\t_5e6\sinusoidal_A10_fuerte"
pam3 = PebbleAccretionModule3.from_datadir(DATADIR, M_star=1.0)
M = pam3.M_EARTH * 1e-4

for r_au in [1.0, 3.0, 5.0, 10.0, 100.0]:
    r = r_au * pam3.AU
    p = pam3._local(20, r)
    St = p['St']
    eta = p['eta']
    M_onset = St * eta**3 * pam3.M_star
    
    print(f"\nr = {r_au} AU")
    print(f"M_onset = {M_onset/pam3.M_EARTH:.2e} Me")
    
    if M < M_onset:
        print("M < M_onset: REGIMEN SAFRONOV")
    else:
        print("M >= M_onset: PEBBLE ACCRETION")
        
    rate_saf = pam3._accretion_rate(M, r, 20)
    print(f"Rate en M0 = {rate_saf:.3e} g/s")
    
    rate_peb = pam3._accretion_rate(M_onset * 1.5, r, 20)
    print(f"Rate justo supra M_onset (Pebble) = {rate_peb:.3e} g/s")
