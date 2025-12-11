PseudoBlast 

PseudoBlast è un tool avanzato per l'allineamento e la ricerca di sequenze biologiche, ispirato all'algoritmo BLAST.
Questa nuova versione offre un'architettura modulare e due interfacce di utilizzo: via terminale (TUI) e via browser (Web).

Struttura del Progetto

`run.py`: Script principale per avviare l'applicazione.
`src/`: Contiene il core logico dell'algoritmo di allineamento.
`tui/`: Interfaccia Utente Testuale (per uso da terminale).
`web/`: Interfaccia Web (per uso tramite browser).
`data/`: Contiene i dataset genomici (es. *T. alexandrinum*).

Installazione e Utilizzo

1.  Clona il repository:

    ```bash
    git clone https://github.com/~Directoryutente~/PseudoBlast.git
    cd PseudoBlast
    ```

2.  Installa le dipendenze:
    ```bash
    pip install -r requirements.txt
    ```

3.  Avvia l'applicazione:
    Puoi avviare il programma principale con:
    ```bash
    python run.py
    ```
    *Segui le istruzioni a schermo per scegliere tra modalità Web o TUI.*

---
*Creato da [debusercccp](https://github.com/debusercccp)*
