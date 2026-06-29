# Checklist del Informe

## Introducción
- [X] TODO BIEN POR AHORA
## Marco Teórico
- [X] TODO BIEN POR AHORA

## Métodos
- [ ] Separar los parametros globales del discos de los parametros explorados. 

## Resultados
- [ ] 4.1: No explicar la figura fuera del caption,
- [ ] Fig 4.1: Mejorar la figura, quitar $\alpha = 0.01$ y dejar solo algunos, hacer un grafico mas grande. Arreglar su caption igual porque no refleja bien
- [X] Fig 4.2: Es mejor quitarla, ya sabesmos efecto tiene el gradiente de presion, por la parte de metodo computacional, sera mejor quitarla, no aporta. Quizas sea mejor dejarlo, pero el efecto del stoke number no se ve bien, estaria bien dejar un r_gap para distintas amplitudes y dejar solo uno, cierto? Para bencharmk esta bien todo, pero para el paper no creo.
- [ ] Reorganizar la estructura.
  - Como lo hicimos:
    - Hicimso el pipeline con tripodpy
    - Tenemos nuestro modelo de acreción
    - Tenemos bencharmks para las simulaciones -hovoller, flux, eta, stoke number y eso.
  - Que datos explarmos:
    - Para todos los casos de alpha desde 1e-4 hasta 1e-2 tenemos los casos smooth, consideramos tambien gaps que se generan mas tarde que otro.
    - Gap unicos:
        - Tenemos casos de gaps unicos varianod r_gap y m_gap. para todos los alphas en v_frag 10
        - Se aprecio que para alpha 0.003 cesa el crecimiento, nada alcanza la masa de exito, a no ser que se aumente la masa INICIAL del embrion, este es un parametro imporante.
        - Notamos que casos de v_frag menores, para ciertos alphas y la misma masa del embrion, se produce el mismo efecto que para v_frag mayor. Al igual que si variamos la masa del embrion es tambien un parametro degenerado. Pero no se explora que parametros son los equivalentes, por ejemplo que alpha 0.0005 v_frag 3ms es lo mismo que alhpa 0.001 a v_frag con la misma masa del embrion. Pero es un efecto que se observa, yo creo que deberiamos agregar algunso graficos mas adescripitvos de estos casos, los bubble puede ser? no lo se.
    - Para gaps sinusoidales
        - Se exploran tambien los hovover, demeran agregarse mas pronto en la parte de resultados, muy impornate, al igual que tambien quizas el eta y a_max a stokes number
        - hicismo tambien los bubble de sinusodal, se parecen un poco a los de gaps simples, deberiamos explorar igual que tanto impacata la cantidad de gaps que se tienen.
  - Luego los graficos que hicimos en primera instacia:
    - Teniamos los de el historila de crecimiento, pero teniamos mmuchos casos, quizas deberiamos explorar 1 caso base y variar el v_frag, y a su vez variar el alpha. quisieramos saber por ejemplo en el historial de evolicion de un planeta que tanto varia su vairmo el v_frag y variamos el alpha, podemos hcerlo para 3 alphas y 3 v_frag, por ejemplo 0.0001 (aqui casi todos estan secos, se evidencia que en ninugn caso se acreta despues del snowline, lo cual es relevante para decir que una turbulencia muy baja genera planetas muy secos.) uno intermedio, 0.0005 y 0.001 (aun que aqui podriamos hacer tambien 0.0003 y 0.0005 porque tambien es relevante en un sentido como ese alpha 0.0003 es el que genera los earth-like) 
    - Los heatmapas categoricos, la verdad ocupan muchos espacio, sun utiles visualmente, pero son bastante pesados de leer, el unico que considero relevante podria ser como se ve que hay una zona privilegiada para el vf 10ms porque se ve claramente que una profundidad optima para el crecimiento "excescivo" pero todos generan waterworld, lo cual puede ser relevante, para mostrar que hay como un sweet spot para le crecimeinto, lo cual no estamos inteoduciendo al paper
  - De los graficos poblacionales:
    - Aqui esta toda la historia condesnada de nuestros datos, un grafico bastnate hermoso, pero muy denso.
    - Seria mejro reestructurar nuestros graficos. Para dar el mismo resultado.
- 
- [ ] 
- [ ]
- [ ] 
- [ ]
  
## Discusiones
- [ ]


## Conclusiones
- [ ]

