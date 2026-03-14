# moran
Basic Moran Process Simulator

Authors:
 - Luis Palomero <lpalomerol@gmail.com>
 - Victor Manso <victor.mansov@gmail.com>

Este repositorio busca crear un sistema de Moran. Lo estamos montando siguiendo TDD.

Las estructuras de datos y las funciones se desarrollan de manera incremental, asegurando que cada componente pase sus pruebas antes de avanzar.

Fases (Katas):

- KATA 1: Construir un sistema que reciba una lista con poblaciones
- KATA 2: Implementar en la simulación un 'paso' en el que se actualicen las poblaciones según las reglas de Moran. Vamos a dos componentes, uno que determina el individuo que se reproduce, y otro en el que se determina el que fallece. Vamos a hacer el sistema determinista. Para ello vamos a hacer que las poblaciones tengan un fitness determinista, que se le pasará en el constructor.
- KATA 3: Implementar el seguimiento de la simulación, registrando el estado de la población en cada paso y permitiendo convertirlo a un DataFrame para su análisis.
- KATA 4: Implementar un resumen del seguimiento de la simulación, mostrando la evolución de cada grupo a lo largo del tiempo.
- KATA 5: REFACTOR: Vamos a hacer que la población siga la estructura de "tell, don't ask". De esta manera, separamos la población de las reglas de evolución, encapsulando el posible comportamiento no determinista en un componente separado y testable.