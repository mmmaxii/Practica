@echo off
setlocal
set PY=C:\astro\tripodpy\env_tripod\Scripts\python.exe
set SCRIPT=heatmap_2d_embryo.py

echo.
echo ============================================================
echo  Generando todos los heatmaps 2D (6 figuras)
echo ============================================================

echo.
echo [1/6]  1 Myr  --  fraccion H2O  (inferno)
%PY% %SCRIPT% --dataroot data/1myr --savedir figures/1myr --quantity fh2o --cmap inferno

echo.
echo [2/6]  1 Myr  --  masa final    (viridis)
%PY% %SCRIPT% --dataroot data/1myr --savedir figures/1myr --quantity mass --cmap viridis

echo.
echo [3/6]  5 Myr  --  fraccion H2O  (inferno)
%PY% %SCRIPT% --dataroot data/5myr --savedir figures/5myr --quantity fh2o --cmap inferno

echo.
echo [4/6]  5 Myr  --  masa final    (viridis)
%PY% %SCRIPT% --dataroot data/5myr --savedir figures/5myr --quantity mass --cmap viridis

echo.
echo [5/6]  Comparativa 1Myr vs 5Myr  --  fraccion H2O  (inferno)
%PY% %SCRIPT% --datasets data/1myr:1Myr data/5myr:5Myr --savedir figures/comparativa --quantity fh2o --cmap inferno

echo.
echo [6/6]  Comparativa 1Myr vs 5Myr  --  masa final    (viridis)
%PY% %SCRIPT% --datasets data/1myr:1Myr data/5myr:5Myr --savedir figures/comparativa --quantity mass --cmap viridis

echo.
echo ============================================================
echo  LISTO. Figuras en:
echo    figures\1myr\
echo    figures\5myr\
echo    figures\comparativa\
echo ============================================================
pause
