#!/usr/bin/env python3
import sys
import os
import subprocess
import argparse

def launch_tui():
    """Lancia l'interfaccia terminale."""
    print("\n[INFO] Avvio interfaccia Terminale (TUI)...")
    # Usa sys.executable per assicurarsi di usare lo stesso interprete Python
    script_path = os.path.join("tui", "main.py")
    try:
        subprocess.run([sys.executable, script_path])
    except KeyboardInterrupt:
        print("\n[INFO] Chiusura TUI.")

def launch_web():
    """Lancia l'interfaccia Web Flask."""
    print("\n[INFO] Avvio interfaccia Web...")
    print("[INFO] Apri il browser su: http://127.0.0.1:5000")
    print("[INFO] Premi CTRL+C per fermare il server.\n")
    
    script_path = os.path.join("web", "app.py")
    try:
        subprocess.run([sys.executable, script_path])
    except KeyboardInterrupt:
        print("\n[INFO] Chiusura Server Web.")

def show_menu():
    """Mostra un menu interattivo se non vengono passati argomenti."""
    while True:
        print("\n" + "="*40)
        print("      PSEUDOBLAST LAUNCHER")
        print("="*40)
        print("1) Avvia Interfaccia Terminale (TUI)")
        print("2) Avvia Interfaccia Web")
        print("Q) Esci")
        print("-" * 40)
        
        scelta = input("Scegli un'opzione: ").strip().lower()
        
        if scelta == '1':
            launch_tui()
        elif scelta == '2':
            launch_web()
        elif scelta == 'q':
            print("Uscita.")
            sys.exit(0)
        else:
            print("Scelta non valida.")

if __name__ == "__main__":
    # Configurazione argomenti da linea di comando (opzionale)
    parser = argparse.ArgumentParser(description="Launcher per PseudoBlast")
    parser.add_argument("mode", nargs="?", choices=["tui", "web"], help="Modalità di avvio (tui o web)")
    
    args = parser.parse_args()

    # Se l'utente specifica un modo (es: python run.py web), lancia quello
    if args.mode == "tui":
        launch_tui()
    elif args.mode == "web":
        launch_web()
    else:
        # Altrimenti mostra il menu interattivo
        show_menu()
