#! /usr/bin/env python3

import os, csv, sys, gc

# ============================================================
# PARSER FASTA
# ============================================================
def parse_fasta(path):
    """Legge un file FASTA e ritorna un dizionario {header: sequenza}."""
    records = {}
    header = None
    seq_lines = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seq = "".join(seq_lines).upper()
                    seq = "".join([c for c in seq if c in "ATGC"])
                    records[header] = seq
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

        if header:
            seq = "".join(seq_lines).upper()
            seq = "".join([c for c in seq if c in "ATGC"])
            records[header] = seq

    return records

# ============================================================
# REVERSE COMPLEMENT (Sequenza complementare inversa)
# ============================================================
def revcomp(seq):
    tab = str.maketrans("ATGC", "TACG")
    return seq.translate(tab)[::-1]

# ============================================================
# COSTRUZIONE INDICE K-MER
# ============================================================
def build_kmer_index(genomes, K):
    index = {}
    for filename, genome in genomes.items():
        for header, seq in genome.items():
            L = len(seq)

            # forward
            for i in range(L - K + 1):
                kmer = seq[i:i+K]
                index.setdefault(kmer, []).append((filename, header, i, False))

            # reverse complement
            rc = revcomp(seq)
            for i in range(L - K + 1):
                kmer = rc[i:i+K]
                # pos convertita nella coordinata forward
                pos = len(seq) - i - K
                index.setdefault(kmer, []).append((filename, header, pos, True))
    return index

##############################################################
# VEDE SE UNA SEQUENZA È GIUSTA
##############################################################

def sequenza_valida(seq):
    """Controlla che la sequenza contenga solo basi valide."""
    for c in seq:
        if c not in "ACGT":
            return False
    return True


# ============================================================
# FUNZIONI DI SUPPORTO
# ============================================================
def mismatch_count(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def mismatch_positions(a, b):
    """Ritorna le POSIZIONI dei mismatch."""
    return [i for i, (x, y) in enumerate(zip(a, b)) if x != y]       


def highlight_mismatches(a, b):
    out = []
    for x, y in zip(a, b):
        if x == y:
            out.append(x)
        else:
            out.append(f"[{x}/{y}]")
    return "".join(out)


def choose_k(query_len):
    return max(8, min(15, query_len // 3))

# ============================================================
# SCELTA DI 2 GENOMI DA CONFRONTARE
# ============================================================

def scegli_due_genomi(folder):
    """
    Permette di scegliere due genomi tra quelli disponibili nella cartella.
    Ritorna i due path selezionati.
    """
    specie = [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]

    if len(specie) < 2:
        print("Errore: servono almeno 2 genomi nella cartella.")
        exit()

    print("\nGenomi disponibili:\n")
    for i, s in enumerate(specie, 1):
        print(f"{i}) {s}")

    print("\nSeleziona due genomi da confrontare:\n")

    try:
        s1 = int(input("Primo genoma (numero): "))
        s2 = int(input("Secondo genoma (numero): "))
    except:
        print("Errore di input: devi inserire numeri.")
        exit()

    if s1 == s2:
        print("Errore: devi scegliere due genomi differenti.")
        exit()

    if not (1 <= s1 <= len(specie)) or not (1 <= s2 <= len(specie)):
        print("Errore: scelta fuori range.")
        exit()

    g1 = os.path.join(folder, specie[s1 - 1])
    g2 = os.path.join(folder, specie[s2 - 1])

    print(f"\nHai scelto:\n - {specie[s1 - 1]}\n - {specie[s2 - 1]}\n")

    return g1, g2

# ============================================================
# CONFRONTO TRA DUE GENOMI (fuzzy search)
# ============================================================

def get_windows(seq, window_size=500, overlap=50):
    """
    Generatore che taglia una sequenza lunga in finestre più piccole.
    Restituisce: (start_index, chunk_sequence)
    """
    seq_len = len(seq)
    if seq_len <= window_size:
        yield 0, seq
        return

    # Passo di avanzamento (step)
    step = window_size - overlap
    
    for i in range(0, seq_len, step):
        # Estrai il blocco
        chunk = seq[i : i + window_size]
        yield i, chunk
        
        # Se siamo alla fine, usciamo
        if i + window_size >= seq_len:
            break

def confronta_genomi(genA, genB, score_threshold=100):
    """
    Versione a BLOCCHI (Sliding Window).
    Non calcola mai matrici giganti, ma fa migliaia di piccoli confronti.
    """
    risultati = []
    
    # PARAMETRI DI TUNING PER RASPBERRY
    WINDOW_SIZE = 400   # Lunghezza del pezzetto da analizzare (matrice 400x400 è piccolissima per la CPU)
    OVERLAP = 50        # Sovrapposizione per non perdere match sui bordi
    SEED_K = 6          # Filtro rapido: se non hanno 6 lettere uguali, salta
    
    print(f"\n--- Inizio Ricerca a Blocchi (Window: {WINDOW_SIZE}bp) ---")
    
    # Contatori per statistiche
    confronti_totali = 0
    confronti_saltati = 0
    match_trovati = 0

    # Iteriamo sui file del Genoma A
    for fileA, recA in genA.items():
        for headerA, full_seqA in recA.items():
            
            # 1. Spezzettiamo la SeqA in finestre
            for startA, chunkA in get_windows(full_seqA, WINDOW_SIZE, OVERLAP):
                
                # Pre-calcolo seeds per il blocco A (velocissimo su stringhe piccole)
                try:
                    seeds_A = {chunkA[i:i+SEED_K] for i in range(len(chunkA) - SEED_K + 1)}
                except:
                    continue # Se il chunk è troppo piccolo o vuoto

                # Iteriamo sui file del Genoma B
                for fileB, recB in genB.items():
                    for headerB, full_seqB in recB.items():
                        
                        # 2. Spezzettiamo la SeqB in finestre
                        for startB, chunkB in get_windows(full_seqB, WINDOW_SIZE, OVERLAP):
                            
                            confronti_totali += 1
                            
                            # --- FILTRO RAPIDO (SEED) ---
                            # Controlliamo se c'è almeno un seed in comune
                            has_seed = False
                            # Scansioniamo chunkB a salti per velocità
                            for k in range(0, len(chunkB) - SEED_K + 1, 5):
                                kmer = chunkB[k : k+SEED_K]
                                if kmer in seeds_A:
                                    has_seed = True
                                    break
                            
                            if not has_seed:
                                confronti_saltati += 1
                                continue
                            
                            # --- ALLINEAMENTO REALE (Smith-Waterman) ---
                            # Ora la matrice è piccola (es. 400x400), il Raspberry la fa in millisecondi
                            try:
                                score, qa, sa = smith_waterman(chunkA, chunkB)
                                
                                if score >= score_threshold:
                                    # Salviamo il risultato aggiustando la posizione globale
                                    match_trovati += 1
                                    risultati.append({
                                        "fileA": fileA,
                                        "headerA": headerA,
                                        "posA_start": startA, # Dove inizia nel genoma intero
                                        "fileB": fileB,
                                        "headerB": headerB,
                                        "posB_start": startB,
                                        "score": score,
                                        "q_aln": qa, # Frammento allineato
                                        "s_aln": sa
                                    })
                                    
                                    # Feedback immediato se trovo qualcosa di buono
                                    if score > score_threshold * 1.5:
                                        sys.stdout.write(f"\n[!] MATCH SCORE {score} | {fileA}:{startA} <-> {fileB}:{startB}")

                            except Exception as e:
                                print(f"Errore catturato: {e}")
                                # Se qualcosa va storto in un blocco, non fermare tutto
                                pass

                            # Feedback progressivo (sovrascrive la riga)
                            if confronti_totali % 2000 == 0:
                                perc_skip = (confronti_saltati / confronti_totali) * 100
                                sys.stdout.write(f"\rAnalisi blocchi: {confronti_totali} | Skip: {perc_skip:.1f}% | RAM Clean...")
                                sys.stdout.flush()
                                
                                # --- PULIZIA RAM ---
                                # Fondamentale su Raspberry: forziamo la pulizia ogni 2000 cicli
                                gc.collect()

    print("\n\nAnalisi completata.")
    print(f"Totale blocchi analizzati: {confronti_totali}")
    print(f"Match trovati: {len(risultati)}")
    
    # Ordiniamo per score decrescente
    return sorted(risultati, key=lambda x: -x["score"])

# ============================================================
# RICERCA ESATTA + MISMATCH
# ============================================================
def search_sequence(query, index, genomes, K, max_mismatch):
    query = query.upper()
    LQ = len(query)

    if LQ < K:
        raise ValueError(f"La query deve essere almeno lunga K={K} basi.")

    seed1 = query[:K]
    mid = (LQ // 2) - (K // 2)
    seed2 = query[mid:mid+K]

    candidates = []
    if seed1 in index:
        candidates.extend(index[seed1])
    if seed2 in index:
        candidates.extend(index[seed2])

    candidates = list(set(candidates))
    results = []

    for filename, header, pos, is_rev in candidates:
        seq = genomes[filename][header]
        working_seq = revcomp(seq) if is_rev else seq

        if pos + LQ <= len(working_seq):
            segment = working_seq[pos:pos + LQ]
            mismatches = mismatch_count(segment, query)
            if mismatches <= max_mismatch:
                results.append({
                    "file": filename,
                    "header": header,
                    "pos": pos,
                    "rev": is_rev,
                    "mismatches": mismatches,
                    "segment": segment,
                    "mm_positions": mismatch_positions(segment, query)
                })

    return sorted(results, key=lambda x: x["mismatches"])

# ============================================================
# SMITH–WATERMAN (Fuzzy Search)
# ============================================================
def smith_waterman(query, target, match=2, mismatch=-1, gap=-1):
    """
    Algoritmo Smith-Waterman per allineamento locale.
    """
    m, n = len(query), len(target)

    # matrice punteggi
    H = [[0]*(n+1) for _ in range(m+1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m+1):
        for j in range(1, n+1):
            score_diag = H[i-1][j-1] + (match if query[i-1] == target[j-1] else mismatch)
            score_up = H[i-1][j] + gap
            score_left = H[i][j-1] + gap

            H[i][j] = max(0, score_diag, score_up, score_left)

            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = (i, j)

    # traceback
    q_align = []
    t_align = []
    i, j = max_pos

    while i > 0 and j > 0 and H[i][j] > 0:
        if H[i][j] == H[i-1][j-1] + (match if query[i-1] == target[j-1] else mismatch):
            q_align.append(query[i-1])
            t_align.append(target[j-1])
            i -= 1
            j -= 1
        elif H[i][j] == H[i-1][j] + gap:
            q_align.append(query[i-1])
            t_align.append("-")
            i -= 1
        else:
            q_align.append("-")
            t_align.append(target[j-1])
            j -= 1

    return max_score, "".join(reversed(q_align)), "".join(reversed(t_align))

#
#  RICERCA FUZZY
#

def fuzzy_search(query, genomes, min_score=10):
    """
    Cerca la 'query' dentro il genoma usando finestre scorrevoli.
    Ottimizzato per Raspberry Pi (evita crash di memoria).
    """
    results = []
    
    # Parametri di sicurezza
    WINDOW_SIZE = 600   # Analizziamo il genoma a pezzi di 600bp
    OVERLAP = 100       # Sovrapposizione
    SEED_K = 5          # Basta un match di 5 lettere per attivare il controllo approfondito
    
    # Pre-calcoliamo i k-mer della query (lo facciamo UNA volta sola)
    query_seeds = set()
    if len(query) >= SEED_K:
        query_seeds = {query[i:i+SEED_K] for i in range(len(query) - SEED_K + 1)}
    
    print(f"\n--- Ricerca Fuzzy Ottimizzata (Query: {len(query)}bp) ---")
    
    count = 0
    skipped = 0
    
    # Iteriamo su tutto il genoma
    for file_name, records in genomes.items():
        for header, sequence in records.items():
            
            # Usiamo get_windows che hai già aggiunto prima!
            for start_pos, chunk in get_windows(sequence, WINDOW_SIZE, OVERLAP):
                count += 1
                
                # --- 1. FILTRO RAPIDO (Seed) ---
                # Se la query ha dei seed, controlliamo se il chunk ne contiene almeno uno.
                if query_seeds:
                    has_match = False
                    # Scansioniamo il chunk a passi di 5 per velocità
                    for k in range(0, len(chunk) - SEED_K + 1, 5):
                        kmer = chunk[k : k+SEED_K]
                        if kmer in query_seeds:
                            has_match = True
                            break
                    
                    if not has_match:
                        skipped += 1
                        continue # Salta il calcolo pesante
                
                # --- 2. CALCOLO SW (Solo se passa il filtro) ---
                try:
                    score, qa, sa = smith_waterman(query, chunk)
                    
                    if score >= min_score:
                        results.append({
                            "file": file_name,
                            "header": header,
                            "pos_start": start_pos, # Salviamo dove inizia il blocco
                            "score": score,
                            "q_aln": qa,
                            "s_aln": sa
                        })
                        # Feedback immediato se trovi un buon match
                        if score > min_score * 2:
                             sys.stdout.write(f"\n[!] MATCH: {header[:15]}... (Pos {start_pos}) Score: {score}")

                except Exception as e:
                    print(f"Errore nel blocco: {e}")
                    pass

                # Feedback visivo e pulizia RAM ogni 1000 blocchi
                if count % 1000 == 0:
                    sys.stdout.write(f"\rAnalizzati {count} blocchi... (Skip: {skipped})")
                    sys.stdout.flush()
                    gc.collect()

    print(f"\nAnalisi finita. Match trovati: {len(results)}")
    return sorted(results, key=lambda x: -x["score"])

# ==========================================================
# CARICAMENTO GENOMI
# ============================================================
def load_genomes(folder):
    print(f"\n--- DEBUG: Analisi cartella '{os.path.basename(folder)}' ---")
    
    if not os.path.exists(folder):
        print(f"ERRORE: La cartella {folder} non esiste sul disco.")
        return {}

    genomes = {}
    files_in_folder = os.listdir(folder)
    print(f"📂 File trovati nella cartella: {files_in_folder}")

    # Estensioni accettate
    valid_extensions = (".fa", ".fasta", ".txt", ".fna", ".faa")
    loaded_count = 0

    for file in files_in_folder:
        # Ignora file nascosti (che iniziano con punto)
        if file.startswith("."):
            continue

        if file.endswith(valid_extensions):
            path = os.path.join(folder, file)
            print(f"  Provo a leggere: {file} ... ", end="")
            
            try:
                records = parse_fasta(path)
                if records:
                    genomes[file] = records
                    print(f"OK! ({len(records)} sequenze)")
                    loaded_count += 1
                else:
                    print(" VUOTO o formato errato (nessuna sequenza >Header trovata).")
            except Exception as e:
                print(f" ERRORE LETTURA: {e}")
        else:
            print(f" Salto '{file}' (estensione non valida o file compresso)")

    if loaded_count == 0:
        print(" ATTENZIONE: Nessun genoma caricato. Controlla le estensioni dei file!")
        print(f"   (Il programma accetta solo: {valid_extensions})")
        print("   (Se hai file .gz, devi prima estrarli/decomprimerli)")

    return genomes
# ============================================================
# SALVATAGGIO CSV
# ============================================================
def save_csv(results, filename="risultati.csv"):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "header", "pos", "reverse", "mismatches", "segment", "mismatch_positions"])
        for r in results:
            writer.writerow([r["file"], r["header"], r["pos"], r["rev"], r["mismatches"], r["segment"], r["mm_positions"]])

    print(f"\nRisultati salvati in {filename}")


# ============================================================
# SCELTA SPECIE
# ============================================================
def lista_specie(folder):
    specie = [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]

    if not specie:
        print("Nessuna specie trovata.")
        exit()

    print("\nSpecie disponibili:\n")
    for i, s in enumerate(specie, 1):
        print(f"{i}) {s}")

    scelta = int(input("\nScegli la specie (numero): "))
    return os.path.join(folder, specie[scelta - 1])

# ============================================================
# COMPARAZIONE GENOMA vs GENOMA 
# ============================================================

def compare_genomes(genA, genB, score_threshold=20):
    """
    Confronta due genomi completi usando Smith-Waterman.
    genA e genB sono dizionari generati da load_genomes().
    Ritorna una lista di allineamenti significativi.
    """

    results = []

    for fileA, genomeA in genA.items():
        for headerA, seqA in genomeA.items():

            for fileB, genomeB in genB.items():
                for headerB, seqB in genomeB.items():

                    score, qa, sa = smith_waterman(seqA, seqB)

                    if score >= score_threshold:
                        results.append({
                            "fileA": fileA,
                            "headerA": headerA,
                            "fileB": fileB,
                            "headerB": headerB,
                            "score": score,
                            "q_aln": qa,
                            "s_aln": sa
                        })

    return sorted(results, key=lambda x: -x["score"])

