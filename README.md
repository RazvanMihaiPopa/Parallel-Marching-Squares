# tema1

    Tema 1 APD

  Scopul temei a fost paralelizarea si optimizarea unui program care genereaza
contururi pentru harti topologice folosind algoritmul Marching Squares.
  Pentru asta, am declarat o structura numita 'function_arguments' care contine
toate variabilele de care am nevoie in functia pasata threadurilor. Am declarat
structura, am creat threadurile si am alocat memorie unde a fost nevoie.
  In functia f, intai am paralelizat rescale_image. M-am folosit de impartirea
thread-urilor cu start si end, asa cum am facut si la laborator. Dupa ce se
executa for-urile, urmeaza o bariera pentru a proteja 'args->image', deoarece
urma apoi sa ii dezaloc memoria. Am pus o bariera si dupa dezalocare, pentru
ca urma ca image sa fie initializat cu new_image.
  In continuare, pe acelasi principiu cu start si end si bariere functioneaza
si urmatoarele 2 functii. Am pus si o bariera intre sample_grid si march
deoarece march foloseste grid-ul abia calculat in functia anterioara.
  La sfarsit, se da join la thread-uri, se distruge bariera si extrag imaginea
pe care vrem sa o scriem din structura.
