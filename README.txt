Progetto di Linguaggi di Programmazione - PROLOG

*** Lorenzo Bini  ***  Matr.806980 ***

Tutte le funzioni richieste, eccetto as_monomial e as_polynomial,
sono state implementate in modo tale da accettare in input sia monomi e
polinomi in forma tradizionale che monomi e polinomi già sottoposti a parsing.


Monomio tradizionale:    5*y
Monomio parsato:         m(5, 1, [v(1, y)])
Polinomio tradizionale:  3*x^2+3
Polinomio parsato:       poly([m(3, 0, []), m(3, 2, [v(2, x)])])


Le funzioni as_monomial e as_polynomial richiedono in input rispettivamente
un monomio ed un polinomio in forma tradizionale, che verrà restituito 
parsato, normalizzato ed ordinato.

Le funzioni polyplus, polyminus e polytimes gestiscono qualsiasi tipo di 
input se e solo se entrambi i valori passati sono della medesima tipologia.
(i.e. Passare un polinomio tradizionale ed un monomio parsato o qualsiasi altra
combinazione produce risultato FALSE).

Tutte le funzioni richieste gestiscono il polinomio normalizzando  
e ordinando secondo i criteri previsti sia le varibili di ogni 
monomio sia monomi di un polinomio, in entrambi i tipi di input.

Si assume che i polinomi forniti in input in forma già parsata siano corretti, 
al più non normalizzati e non ordinati.

Nel progetto lo 0 ("zero") viene considerato a sé stante come:
m(0, 0, [])
Se lo zero è presente all'interno del polinomio esso viene ignorato.
Uguale condizione si verifica se durante la normalizzazione si ottengono
monomi con valore matematico uguale a zero.

Esempi di input:
0   --->  poly([m(0, 0, [])])
x+0  --->  poly([m(1, 1, [v(1, x)])])
3*x-3*x+2*y  --->  poly([m(2, 1, [v(1, y)])])
3*x*0+2  --->  poly([m(2, 0, [])])


Nel progetto i simboli di variabile elevati a zero vengono interpretati come 1.
Uguale condizione si verifica se durante la normalizzazione si ottengono come
risultato simboli di variabile elevati a zero.

Esempi di input:
x^0  --->  poly([m(1, 0, [])])
3*x^0+2  --->  poly([m(5, 0, [])])
x^2*x^(-2)*a  --->  poly([m(1, 1, [v(1, a)])]


****Informazioni sull'ordinamento****

Per l'ordinamento è stato importato un algoritmo Merge Sort profondamente
modificato al fine di essere adattato alle esigenze del programma.
L'algoritmo merge sort prevede codice dedicato all'ordinamento delle variabili
di un monomio e codice dedicato all'ordinamento dei monomi di un polinomio.

I simboli di variabile di un monomio vengono ordinati in ordine lessicografico
crescente.
I monomi di un polinomio vengono inizialmente ordinati in ordine crescente 
in base al grado complessivo degli stessi con spareggi determinati dai  
simboli di variabile.
L'ordinamento di due monomi con stessi simboli di variabile viene effettuato
in maniera crescente rispetto alle combinazioni variabile/esponente.


****Informazioni su funzioni secondarie****

Sono state implementate delle funzioni secondarie per facilitare la stesura
di quelle primarie e garantire una maggiore facilità di lettura.
Le funzioni secondarie tendenzialmente richiedono in input una lista di monomi
già sottoposti a ordinamento e talvolta a normalizzazione (requisito fondamentale).
Sottoporre alle funzioni secondarie una lista di monomi non ordinata potrebbe
portare ad un output inatteso od errato.

Tutte le funzioni che possono ricevere in input un polinomio in forma 
tradizionale si servono della funzione is_trad_poly\1 che permette, appunto,
il loro riconoscimento.

Le funzioni che trattano sottrazioni tra polinomi si servono della funzione
invert_coefficient\2 per invertire il segno del coefficiente e permettere una
corretta normalizzazione.
 