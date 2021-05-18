file LABO.m

Corrisponde al main

PART 1
Lanciandolo riproduce i risultati del test numerico selezionato nella variabile Test (Test1 reference case Test2 problema con soluzione esatta) per valori di refinment 3, 4, 5, e 6.
Decidendo di attivare la variabile debug_plot è possibile vedere anche graficamente l'azione del penalty operator.
Decidendo di attivare la variabile incr_debug è possibile vedere l'errore incrementale
Decidendo di attivare la variabile cavitat_area_debug è possibile vedere l'evoluzione della cavitation area

NOTA I file delle simulazioni sono già salvati, se si volessero ripetere simulazioni con altri parametri o con altri test e si volessero salvara, decommentare linea 139

PART 2

Usando i file già salvati inclusi in directory attuale, esegue test convergenza con epsilon = h^2

PART 3

Usando i file già salvati inclusi in directory attuale, esegue test convergenza per vari epsilon

PART 4

Materiale EXTRA sul Lagrangian Augmented:utile se si volessero usare BC quadratiche o per stabilizzare problemi i cui autovettori del obstacle problem fossero molto vicini tra loro 
dando luogo ad instabilità

Cartella Numerical Error

Test dei risultati di convergenza usando approccio numerico (sufficiente runnare file estensione_err_H1.m)

Cartella Exponential interpolation

Nel file max_cavitationarea.m viene caricato il file max_contact.mat, risultato di funzione commentata sopra, e restituisce errore puntuale massimo su area cavitazione




