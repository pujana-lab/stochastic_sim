# moran
Basic Moran Process Simulator

Authors:
 - Luis Palomero <lpalomerol@gmail.com>
 - Victor Manso <victor.mansov@gmail.com>

Este repositorio busca crear un sistema de Moran. Lo estamos montando siguiendo TDD.

Las estructuras de datos y las funciones se desarrollan de manera incremental, asegurando que cada componente pase sus pruebas antes de avanzar.

Fases:

- Fase 0: Construir un sistema que reciba una lista con poblaciones
- Fase 1: Implementar en la simulación un 'paso' en el que se actualicen las poblaciones según las reglas de Moran. Vamos a dos componentes, uno que determina el individuo que se reproduce, y otro en el que se determina el que fallece. Vamos a hacer el sistema determinista. Para ello vamos a hacer que las poblaciones tengan un fitness determinista, que se le pasará en el constructor.