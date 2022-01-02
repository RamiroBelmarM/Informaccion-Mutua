# Informaccion Mutua
Información mutua para estimar el time delay de una serie de tiempo

```cpp
#include "muttinf.h"
int main() {
    //datos bajo estudio
    //debe de tener al menos 1000 observaciones
    std::vector<double> dat; 

    //dt, tiempo de muestreo, se ignoran los primeras observaciones 
    //que cumplan t < 50   

    //se crea el objeto Mutual Information
    //a sus argumentos se pasa los datos y dt
    MuttInf test(dat, dt);     

    //se inicia el test, limite superior es el maximo valor de tus datos
    // y limite inferior, lo contrario.
    test.init(limite_sup, limite_inf);

    //se hace el test
    test.make();

    //por si quieres los datos de la informacion mutua
    std::vector<double> resultado=test.get_dat();

    //por si solo quieres obtener el time delay
    //ten cuidado es mas recomendable inspeccionar primero lo anterior
    double tau=test.get_pos();

    return 0;
}
```

# Referencia
Kliková, B., & Raidl, A. (2011, May). Reconstruction of phase space of dynamical systems using method of time delay. In Proceedings of the 20th Annual Conference of Doctoral Students-WDS (pp. 83-87).

Nichols, J. M., & Nichols, J. D. (2001). Attractor reconstruction for non-linear systems: a methodological note. Mathematical biosciences, 171(1), 21-32.