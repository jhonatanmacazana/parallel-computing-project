# Oscilador Armónico Acoplado

> Autores: Jhonatan Macazana, Samir Muñoz

## Compilación

Asegurarse de contar con el paquete de librerias de C++. 

``` bash
# en Ubuntu
sudo apt install build-essentials make
```

Compilar el código

``` bash
make all
```

## Ejecución en local

Ejecutar el código secuencial

``` bash
make run
```

## Ejecución en cluster

## Enunciado

1. Analice el caso de un defecto en la cadena. Escriba un programa en paralelo que resuelva la matriz y calcule los valores y vectores propios para N=99, n0=60 y n0=30 con ayuda de la librería `tqli` de Numerical Recipees.

2. Escoja condiciones iniciales y grafique x(t) para un intervalo adecuado de tiempo.

Se propone lo siguiente:

a) Desarrolle puntos 1 y 2 de la descripción. Elevar la dimensión de la cadena
lo suficiente para obtener mediciones significativas de tiempos.

b) Registrar el desarrollo de las simulaciones de manera ordenada en por lo
menos 3 pasos (versiones beta), hasta llegar al código final.

c) Estimar la complejidad del código secuencial y en paralelo, y compararla con
métricas de tiempo presentadas en gráficas adecuadas

d) Medir la velocidad (en FLOPs) del algoritmo y representarla en gráficas. Medir la escalabilidad del software

e) Optimizar el software con un desarrollo orientado al paralelismo, y presentarlo en las conclusiones del paper.
