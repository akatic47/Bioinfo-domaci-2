import os
import argparse
import subprocess
import json

def setup_directories():
    """Kreira sve neophodne foldere za izlazne podatke, osim ako već ne postoje."""
    dirs = [
        "01_fastqc", 
        "02_alignment", 
        "03_processing", 
        "04_variant_calling", 
        "logs"
    ]
    for d in dirs:
        # exist_ok=True je magija koja sprečava pucanje ako folder već postoji!
        os.makedirs(d, exist_ok=True)
        print(f"[INFO] Proveren/kreiran folder: {d}")

def run_command(command, log_file):
    """Pokreće terminalsku komandu i sav ispis (stdout i stderr) čuva u log fajl."""
    with open(log_file, "w") as log:
        # check=True znači da ako alat pukne, Python odmah prekida skriptu
        subprocess.run(command, shell=True, stdout=log, stderr=subprocess.STDOUT, check=True)

def run_fastqc(fastq1, fastq2):
    """Pokreće FastQC analizu samo ako izveštaji već ne postoje."""
    print("\n--- Korak 1: Kontrola kvaliteta (FastQC) ---")
    
    # FastQC po defaultu uzima ime fajla i dodaje '_fastqc.html' na kraj
    # npr. sample_0.chrom11.exome.pe1.fq.gz postaje sample_0.chrom11.exome.pe1_fastqc.html
    base1 = os.path.basename(fastq1).replace(".fq.gz", "_fastqc.html")
    base2 = os.path.basename(fastq2).replace(".fq.gz", "_fastqc.html")
    
    out1 = os.path.join("01_fastqc", base1)
    out2 = os.path.join("01_fastqc", base2)
    
    # IDEMPOTENTNOST: Ako fajlovi već postoje, ne radi ništa!
    if os.path.exists(out1) and os.path.exists(out2):
        print("[PRESKOČENO] FastQC izveštaji već postoje. Idemo dalje!")
        return

    print("[POKRETANJE] Generišem FastQC izveštaje...")
    cmd = f"fastqc {fastq1} {fastq2} -o 01_fastqc"
    log_path = "logs/01_fastqc.log"
    
    try:
        run_command(cmd, log_path)
        print(f"[USPEH] FastQC završen. Pogledaj izveštaje u folderu 01_fastqc.")
    except subprocess.CalledProcessError:
        print(f"[GREŠKA] FastQC je prijavio problem! Detalji su u fajlu: {log_path}")
        exit(1)

def run_bwa(fastq1, fastq2, ref_genome):
    """Indeksira referentni genom (ako nije) i mapira read-ove."""
    print("\n--- Korak 2: Mapiranje (BWA & Samtools) ---")
    
    # 1. INDEKSIRANJE REFERENTNOG GENOMA
    # BWA alat traži fajl koji se završava na .bwt da bi znao da li je indeksiranje gotovo
    ref_bwt = ref_genome + ".bwt"
    if os.path.exists(ref_bwt):
        print(f"[PRESKOČENO] Genom {ref_genome} je već indeksiran.")
    else:
        print("[POKRETANJE] Indeksiranje genoma (ovo može potrajati sat-dva za ceo hg38!)...")
        cmd_index = f"bwa index {ref_genome}"
        run_command(cmd_index, "logs/bwa_index.log")
        print("[USPEH] Indeksiranje referentnog genoma završeno.")

    # 2. MAPIRANJE I SORTIRANJE
    out_bam = os.path.join("02_alignment", "sample_0.sorted.bam")
    
    # IDEMPOTENTNOST: Ako krajnji BAM fajl već postoji, preskačemo mapiranje
    if os.path.exists(out_bam):
        print("[PRESKOČENO] Mapiranje je već urađeno. Idemo dalje!")
        return

    print("[POKRETANJE] BWA mapiranje i Samtools sortiranje...")
    # Trik sa 'pipe': bwa mem izbacuje SAM format -> samtools ga pretvara u BAM -> pa ga sortira
    cmd_mem = f"bwa mem -R \"@RG\\tID:1\\tSM:sample_0\\tPL:ILLUMINA\" {ref_genome} {fastq1} {fastq2} | samtools view -Sb - | samtools sort -o {out_bam}"
    
    try:
        run_command(cmd_mem, "logs/bwa_mem.log")
        
        # GATK alatima će kasnije trebati i indeks ovog našeg novog BAM fajla
        print("[POKRETANJE] Indeksiranje BAM fajla (.bai)...")
        run_command(f"samtools index {out_bam}", "logs/samtools_index.log")
        
        print("[USPEH] Mapiranje uspešno završeno. Fajl je u folderu 02_alignment.")
    except subprocess.CalledProcessError:
        print("[GREŠKA] BWA mapiranje je puklo! Proveri fajl: logs/bwa_mem.log")
        exit(1)

def run_processing(ref_genome, known_sites):
    """Priprema indekse za GATK, uklanja duplikate i kalibriše kvalitet (BQSR)."""
    print("\n--- Korak 3: Čišćenje i kalibracija (GATK) ---")
    
    bam_in = os.path.join("02_alignment", "sample_0.sorted.bam")
    dedup_bam = os.path.join("03_processing", "sample_0.dedup.bam")
    bqsr_bam = os.path.join("03_processing", "sample_0.bqsr.bam")
    
    # IDEMPOTENTNOST: Gledamo da li je finalni BQSR BAM već napravljen
    if os.path.exists(bqsr_bam):
        print("[PRESKOČENO] GATK čišćenje i kalibracija su već urađeni. Idemo dalje!")
        return

    # --- 1. PRIPREMA REFERENCI ZA GATK ---
    ref_dict = ref_genome.replace(".fa", ".dict")
    if not os.path.exists(ref_dict):
        print("[POKRETANJE] Pravljenje FASTA indeksa i rečnika za GATK...")
        run_command(f"samtools faidx {ref_genome}", "logs/samtools_faidx.log")
        run_command(f"./gatk-4.5.0.0/gatk CreateSequenceDictionary -R {ref_genome}", "logs/gatk_dict.log")
        
    known_idx = known_sites + ".csi" # Ako je BCF, indeks je obično .csi
    if not os.path.exists(known_idx) and not os.path.exists(known_sites + ".tbi"):
        print("[POKRETANJE] Indeksiranje fajla sa poznatim varijantama...")
        run_command(f"./gatk-4.5.0.0/gatk IndexFeatureFile -I {known_sites}", "logs/gatk_index_bcf.log")

    # --- 2. MARK DUPLICATES ---
    if not os.path.exists(dedup_bam):
        print("[POKRETANJE] Uklanjanje PCR duplikata (MarkDuplicates)...")
        metrics = os.path.join("03_processing", "sample_0.metrics.txt")
        cmd_markdup = f"./gatk-4.5.0.0/gatk MarkDuplicates -I {bam_in} -O {dedup_bam} -M {metrics}"
        run_command(cmd_markdup, "logs/gatk_markdup.log")
        run_command(f"samtools index {dedup_bam}", "logs/samtools_index_dedup.log")

    # --- 3. BQSR (Base Quality Score Recalibration) ---
    print("[POKRETANJE] Kalibracija kvaliteta (BQSR)...")
    recal_table = os.path.join("03_processing", "sample_0.recal_data.table")
    
    # Prvi deo: Računanje grešaka
    cmd_recal = f"./gatk-4.5.0.0/gatk BaseRecalibrator -I {dedup_bam} -R {ref_genome} --known-sites {known_sites} -O {recal_table}"
    run_command(cmd_recal, "logs/gatk_baserecal.log")
    
    # Drugi deo: Primena korekcije
    cmd_apply = f"./gatk-4.5.0.0/gatk ApplyBQSR -I {dedup_bam} -R {ref_genome} --bqsr-recal-file {recal_table} -O {bqsr_bam}"
    run_command(cmd_apply, "logs/gatk_applybqsr.log")
    
    # Konačno indeksiranje čistog fajla
    run_command(f"samtools index {bqsr_bam}", "logs/samtools_index_bqsr.log")
    
    print("[USPEH] Čišćenje i kalibracija završeni. Čist fajl je u folderu 03_processing.")

def run_variant_calling(ref_genome, bam_file):
    print("\n--- Korak 4: Pronalaženje varijanti (Variant Calling) ---")
    raw_vcf = "04_variant_calling/sample_0.raw.vcf.gz"
    filtered_vcf = "04_variant_calling/sample_0.filtered.vcf.gz"
    pass_vcf = "04_variant_calling/sample_0.pass.vcf.gz"

    # 1. HaplotypeCaller (Pronalaženje apsolutno svih mutacija)
    if not os.path.exists(raw_vcf):
        print("[POKRETANJE] GATK HaplotypeCaller (traženje mutacija)...")
        cmd_hc = f"./gatk-4.5.0.0/gatk HaplotypeCaller -R {ref_genome} -I {bam_file} -O {raw_vcf}"
        run_command(cmd_hc, "logs/gatk_haplotypecaller.log")
    else:
        print("[PRESKOČENO] HaplotypeCaller je već odrađen.")

    # 2. VariantFiltration (Hard filtering - obeležavanje loših varijanti)
    if not os.path.exists(filtered_vcf):
        print("[POKRETANJE] GATK VariantFiltration (Hard Filtering)...")
        # Standardni GATK filteri (Napomena: ako vam je na vežbama profa dao 
        # neke specifične brojeve za QD, FS, MQ, slobodno ih promeni ovde)
        filter_expr = '"QD < 2.0 || FS > 60.0 || MQ < 40.0"'
        cmd_filter = f"./gatk-4.5.0.0/gatk VariantFiltration -R {ref_genome} -V {raw_vcf} -O {filtered_vcf} --filter-expression {filter_expr} --filter-name \"HardFilter\""
        run_command(cmd_filter, "logs/gatk_variantfiltration.log")
    else:
        print("[PRESKOČENO] Filtriranje je već odrađeno.")

    # 3. SelectVariants (Izdvajanje konačnih, sigurnih mutacija)
    if not os.path.exists(pass_vcf):
        print("[POKRETANJE] Izdvajanje PASS varijanti...")
        cmd_pass = f"./gatk-4.5.0.0/gatk SelectVariants -R {ref_genome} -V {filtered_vcf} -O {pass_vcf} --exclude-filtered"
        run_command(cmd_pass, "logs/gatk_selectvariants.log")
        print("[USPEH] Sve mutacije su uspešno pronađene i filtrirane!")
    else:
        print("[PRESKOČENO] Izdvajanje PASS varijanti je već odrađeno.")

    return raw_vcf, filtered_vcf, pass_vcf

def main():
    # 1. Čitanje argumenata iz komandne linije
    parser = argparse.ArgumentParser(description="Automatizovani Variant Calling Pipeline")
    parser.add_argument("--fastq1", required=True, help="Putanja do prvog (R1) FASTQ fajla")
    parser.add_argument("--fastq2", required=True, help="Putanja do drugog (R2) FASTQ fajla")
    parser.add_argument("--ref", required=True, help="Putanja do referentnog genoma (.fa)")
    parser.add_argument("--known", required=True, help="Putanja do baze poznatih varijanti (.bcf)")
    
    args = parser.parse_args()

    print("=== Pokretanje pipeline-a ===")
    
    # 2. Priprema strukture
    setup_directories()
    
    run_fastqc(args.fastq1, args.fastq2) #1
    run_bwa(args.fastq1, args.fastq2, args.ref) #2
    run_processing(args.ref, args.known) #3
    
    #4
    bam_bqsr = "03_processing/sample_0.bqsr.bam"
    raw_vcf, filtered_vcf, pass_vcf = run_variant_calling(args.ref, bam_bqsr)

    manifest = {
        "fastqc_1": "01_fastqc/sample_0.chrom11.exome.pe1_fastqc.html",
        "fastqc_2": "01_fastqc/sample_0.chrom11.exome.pe2_fastqc.html",
        "bam_aligned": "02_alignment/sample_0.sorted.bam",
        "bam_dedup": "03_processing/sample_0.dedup.bam",
        "bam_bqsr": bam_bqsr,
        "vcf_raw": raw_vcf,
        "vcf_filtered": filtered_vcf,
        "vcf_pass": pass_vcf
    }
    
    with open("results_manifest.json", "w") as f:
        json.dump(manifest, f, indent=4)
    print("\n[USPEH] Generisan results_manifest.json fajl! Pajplajn je GOTOV!")

    print("\n[INFO] Svi pokrenuti koraci su uspešno završeni!")

if __name__ == "__main__":
    main()