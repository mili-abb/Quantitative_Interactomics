import os
from file_op import read_csv_file, extract_species_and_bait
from calcul import calculate_Kapp, get_bait_concentration, add_complexome_column

def main():
    # === Étape 1 : Charger le fichier de données ===
    data_path = "essai/input/human_hck_data.csv"  # Nom du fichier à adapter
    df = read_csv_file(data_path)
    if df is None:
        return

    # === Étape 2 : Extraire espèce et bait ===
    species, bait = extract_species_and_bait(data_path)
    print(f"Espèce détectée : {species} | Bait : {bait}")

    # === Étape 3 : Calculer Kapp ===
    df = calculate_Kapp(df)

    # === Étape 4 : Récupérer la concentration du bait ===
    bait_id = df['Uniprot_ID'].iloc[0]  # Hypothèse : bait en première ligne
    bait_conc = get_bait_concentration(df, bait_id)
    print(f"Concentration du bait ({bait_id}) : {bait_conc} µM")

    if bait_conc is None:
        print("❌ Concentration du bait introuvable.")
        return

    # === Étape 5 : Calcul du complexome ===
    df = add_complexome_column(df, bait_conc)

    # === Étape 7 : Sauvegarde des résultats ===
    output_file = f"{species}_{bait}_analyzed.csv"
    df.to_csv(output_file, index=False)
    print(f"✅ Analyse terminée. Fichier sauvegardé : {output_file}")

if __name__ == "__main__":
    main()

