# Estructura del Proyecto

Este repositorio contiene los códigos, datos y documentos relacionados con el estudio de la evolución de discos protoplanetarios, dinámica de polvo, acreción de pebbles y la formación de *waterworlds*. 

A continuación se presenta un mapa de los directorios y archivos principales del proyecto (excluyendo archivos temporales, cachés o datos pesados ignorados en `.gitignore`), junto con una breve explicación de su propósito.

```text
.
├── Primeros acercamientos y archivos antiguos/  ---- Códigos iniciales y versiones antiguas (como PebbleAccretion 1 y 2).
│
├── informe_practica/                            ---- Código fuente en LaTeX para el informe final de práctica.
│   ├── chapters/                                ---- Capítulos individuales del informe en formato .tex.
│   ├── figures/                                 ---- Gráficos, diagramas y mosaicos usados en el documento.
│   ├── main.tex                                 ---- Archivo principal de LaTeX para compilar todo el informe.
│   └── biblio.bib                               ---- Referencias bibliográficas.
│
├── extracted_papers.json                        ---- Base de datos extraída sobre literatura relevante (metadata).
│
└── codigos/                                     ---- Directorio principal del desarrollo del modelo físico.
    │
    ├── PA3Py/                                   ---- Módulo avanzado de Acreción de Pebbles.
    │   └── PebbleAccretion3.py                  ---- Implementación de la física de acreción (Ormel 2017 / Drążkowska 2023).
    │
    ├── pipeline_methods/                        ---- Submódulos que definen la física y estructura del disco.
    │   ├── disk_setup.py                        ---- Configura perfiles iniciales de gas, polvo y parámetros de la grilla.
    │   ├── snowline_physics.py                  ---- Modela la termodinámica y evolución de la línea de nieve de agua.
    │   ├── disk_chemistry.py                    ---- Administra las fracciones de masa y abundancias (H2O, silicatos).
    │   ├── pressure_bumps.py                    ---- Funciones para inyectar perfiles de presión y perturbaciones (gaps).
    │   └── oka_interpolation.py                 ---- Métodos de interpolación basados en el modelo de Oka et al.
    │
    ├── pipeline_analysis/                       ---- Herramientas para leer y analizar los resultados de las simulaciones.
    │   ├── snapshot_analyzer.py                 ---- Extrae y procesa datos de los archivos HDF5 de las simulaciones.
    │   ├── plotters.py                          ---- Funciones base para el ploteo estandarizado de resultados.
    │   └── mcross_analysis.py                   ---- Analizador del cruce de masas de aislamiento y pebble flux.
    │
    ├── scratch/                                 ---- Entorno de pruebas, experimentación y scripts misceláneos.
    │   ├── geryon/                              ---- Scripts "lanzadores" y "simuladores" para correr código en el clúster Geryon.
    │   └── Post oka and hartmann/               ---- Pruebas específicas para reproducir resultados de modelos teóricos previos.
    │
    ├── pipeline_snowlines.py                    ---- Script orquestador principal que corre las simulaciones de evolución del disco.
    ├── generar_benchmarks.py                    ---- Genera los casos base de simulación (smooth, gaps, etc.).
    ├── generar_mosaicos_paper.py                ---- Ensambla gráficos compuestos (mosaicos) para usarlos directamente en el informe.
    │
    ├── plot_*.py                                ---- (Ej: plot_bubble_v3.py, plot_dynamic_snowline.py) Scripts creados para generar figuras específicas y visualizaciones del análisis.
    ├── resultados*.py                           ---- (Ej: resultados_bubble.py, resultados_delayed.py) Procesan las salidas de simulación para calcular y extraer la masa y composición final de los embriones.
    │
    ├── extraer_datos.py                         ---- Extrae información tabular masiva desde los archivos generados.
    └── summary_master_todos_casos.csv           ---- Tabla maestra recopilatoria con el resumen estadístico de las corridas.
```
