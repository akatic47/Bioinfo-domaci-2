import gzip
import matplotlib.pyplot as plt

def analiziraj_varijante(vcf_fajl):
    raw_snps = 0
    raw_indels = 0
    pass_snps = 0
    pass_indels = 0
    failed_variants = 0

    print(f"[INFO] Čitam fajl: {vcf_fajl}...")
    
    # Otvaramo komprimovani VCF fajl kao običan tekst
    with gzip.open(vcf_fajl, "rt") as f:
        for linija in f:
            # Preskačemo zaglavlje fajla
            if linija.startswith("#"):
                continue
                
            delovi = linija.strip().split("\t")
            ref = delovi[3]
            alt = delovi[4]
            filter_status = delovi[6]
            
            # Klasifikacija: Ako su i REF i ALT dužine 1, to je SNP. Inače je INDEL.
            is_snp = (len(ref) == 1 and len(alt) == 1)
            
            # Brojimo sirove varijante (sve koje postoje u fajlu)
            if is_snp:
                raw_snps += 1
            else:
                raw_indels += 1
                
            # Brojimo filtrirane (samo one sa statusom PASS)
            if filter_status == "PASS" or filter_status == ".":
                if is_snp:
                    pass_snps += 1
                else:
                    pass_indels += 1
            else:
                failed_variants += 1

    # Matematika za izveštaj
    total_raw = raw_snps + raw_indels
    total_pass = pass_snps + pass_indels
    failed_percent = (failed_variants / total_raw) * 100 if total_raw > 0 else 0

    # Ispis tražene statistike
    print("\n=== ZADATAK 2: STATISTIKA VARIJANTI ===")
    print(f"Ukupan broj sirovih varijanti: {total_raw}")
    print(f"Ukupan broj PASS varijanti: {total_pass}")
    print(f"Broj SNP (pre / posle filtriranja): {raw_snps} / {pass_snps}")
    print(f"Broj INDEL (pre / posle filtriranja): {raw_indels} / {pass_indels}")
    print(f"Broj varijanti koje nisu prošle filter: {failed_variants}")
    print(f"Procenat odbačenih varijanti: {failed_percent:.2f}%")

    # Pravljenje grafikona koji nosi bodove
    labele = ['Sirove SNP', 'PASS SNP', 'Sirove INDEL', 'PASS INDEL']
    vrednosti = [raw_snps, pass_snps, raw_indels, pass_indels]
    boje = ['#3498db', '#2ecc71', '#e74c3c', '#f1c40f']

    plt.figure(figsize=(9, 6))
    plt.bar(labele, vrednosti, color=boje)
    plt.title('Poređenje SNP i INDEL varijanti pre i posle hard filtriranja')
    plt.ylabel('Broj detektovanih varijanti')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Čuvanje slike
    ime_slike = "grafikon_varijanti.png"
    plt.savefig(ime_slike)
    print(f"\n[USPEH] Grafikon je sačuvan u radnom folderu kao '{ime_slike}'")

if __name__ == "__main__":
    vcf_putanja = "04_variant_calling/sample_0.filtered.vcf.gz"
    analiziraj_varijante(vcf_putanja)