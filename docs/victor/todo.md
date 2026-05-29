# TODO

- quitar la logica de los rates de dentro de los Clone Type 
- Crear objeto reacciones 
- Crear objeto pseudomatriz de interaccion
Esto nos permite generalizar los procesos de propensities y hacerlo escalable, solo tenemos que definir los tipos poblacionales con los que queremos trabajar y una matriz o pseudo matriz (podria ser un excel que podamos editar facil, lo que sea). por defecto todo es cero y simplemente le decimos al programa cual interactua con tal para que solo coja y sume/multiplique lo que toca donde toca. Parecido al planteamiento matricial.

Asi no hay qye programar a mano las interacciones d etodas con todas a mano dentro de las propensities.

Por ejemplo creamos un nuevo tipo de cell y para cada propiedad (que esta ya programada en la expresion comun de birth_rate ) vamos marcando casillas. por ejemplo nuestro nuevo tipo le afecta el crowding de las de cancer y de las de immune pero no las de wild type y las exhaysted tienen accion citolitica sobre ella (ejemplo no biologico). Marcamos las casillas correspondientes.

No se si tiene sentido pero lo hablamos el lunes 

