import os
import argparse
import subprocess
import json
import shutil

def setup_directories(output_dir):
    """Kreira sve neophodne foldere unutar zadatog izlaznog direktorijuma."""
    dirs = [
        "01_fastqc", 
        "02_alignment", 
        "03_processing", 
        "04_variant_calling", 
        "logs"
    ]
    for d in dirs:
        os.makedirs(os.path.join(output_dir, d), exist_ok=True)
    print(f"[INFO] Provereni/kreirani folderi unutar izlaznog direktorijuma: {output_dir}")

def run_command(command, log_file):
    """Pokreće terminalsku komandu i sav ispis (stdout i stderr) čuva u log fajl."""
    with open(log_file, "w") as log:
        subprocess.run(command, shell=True, stdout=log, stderr=subprocess.STDOUT, check=True)

def check_step_done(output_dir, step_name):
    """Proverava da li postoji skriveni marker fajl koji potvrđuje da je korak uspešno završen."""
    return os.path.exists(os.path.join(output_dir, f".{step_name}.done"))

def mark_step_done(output_dir, step_name):
    """Kreira skriveni marker fajl nakon što se korak izvrši bez greške."""
    with open(os.path.join(output_dir, f".{step_name}.done"), "w") as f:
        f.write("DONE")

def clean_failed_outputs(files_list):
    """Uklanja nedovršene privremene fajlove ako alat pukne na pola posla."""
    for f in files_list:
        if os.path.exists(f):
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)

# --- KORAK 1: KONTROLA KVALITETA (FastQC) ---
def run_fastqc(fastq1, fastq2, output_dir):
    print("\n--- Korak 1: Kontrola kvaliteta (FastQC) ---")
    if check_step_done(output_dir, "fastqc"):
        print("[PRESKOČENO] FastQC izveštaji već postoje. Idemo dalje!")
        return

    out_dir_fq = os.path.join(output_dir, "01_fastqc")
    log_path = os.path.join(output_dir, "logs", "01_fastqc.log")
    
    print("[POKRETANJE] Generišem FastQC izveštaje...")
    cmd = f"fastqc {fastq1} {fastq2} -o {out_dir_fq}"
    
    try:
        run_command(cmd, log_path)
        mark_step_done(output_dir, "fastqc")
        print("[USPEH] FastQC završen uspešno.")
    except subprocess.CalledProcessError:
        print(f"[GREŠKA] FastQC je otkazao! Detalji u fajlu: {log_path}")
        exit(1)

# --- KORAK 2: MAPIRANJE (BWA-MEM & Samtools) ---
def run_bwa(fastq1, fastq2, ref_genome, output_dir):
    print("\n--- Korak 2: Mapiranje (BWA-MEM & Samtools) ---") 
    
    required_index_files = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    missing = [ext for ext in required_index_files if not os.path.exists(ref_genome + ext)]
    if missing:
        print(f"[GREŠKA] Nedostaju BWA indeks fajlovi: {missing}")
        print(f"[INFO] Pokreni ručno: bwa index {ref_genome}")
        exit(1)
    else:
        print("[INFO] BWA indeks pronađen. Preskačem indeksiranje.")
    
    if check_step_done(output_dir, "alignment"):
        print("[PRESKOČENO] Mapiranje je već ranije uspešno završeno. Idemo dalje!")
        return

    final_bam = os.path.join(output_dir, "02_alignment", "sample_0.sorted.bam")
    final_bai = final_bam + ".bai"
    tmp_bam = os.path.join(output_dir, "02_alignment", "tmp_sorted.bam")
    tmp_bai = tmp_bam + ".bai"
    tmp_sam = os.path.join(output_dir, "02_alignment", "tmp_aligned.sam")
    log_bwa = os.path.join(output_dir, "logs", "bwa_mem.log")
    log_sort = os.path.join(output_dir, "logs", "samtools_sort.log")
    log_index = os.path.join(output_dir, "logs", "samtools_index.log")

    try:
        print("[POKRETANJE] BWA mapiranje...")
        cmd1 = f"bwa mem -t 2 -R \"@RG\\tID:1\\tSM:sample_0\\tPL:ILLUMINA\" {ref_genome} {fastq1} {fastq2} -o {tmp_sam}"
        run_command(cmd1, log_bwa)

        print("[POKRETANJE] Samtools konverzija i sortiranje...")
        cmd2 = f"samtools view -Sb {tmp_sam} | samtools sort -m 1G -o {tmp_bam}"
        run_command(cmd2, log_sort)
        os.remove(tmp_sam)

        print("[POKRETANJE] Indeksiranje BAM fajla...")
        run_command(f"samtools index {tmp_bam}", log_index)
        
        os.rename(tmp_bam, final_bam)
        os.rename(tmp_bai, final_bai)
        
        mark_step_done(output_dir, "alignment")
        print("[USPEH] Mapiranje i indeksiranje BAM fajla završeno.")
    except subprocess.CalledProcessError:
        print("[GREŠKA] BWA mapiranje je puklo! Čistim nepotpune izlaze...")
        clean_failed_outputs([tmp_sam, tmp_bam, tmp_bai, final_bam, final_bai])
        exit(1)

# --- KORAK 3: GATK OBRADA (MarkDuplicates & BQSR) ---
def run_processing(ref_genome, known_sites, output_dir, interval=None):
    print("\n--- Korak 3: GATK Čišćenje i Kalibracija ---")
    
    # 1a. Pravljenje FASTA indeksa (.fai) ako nedostaje (rešava problem sa UCSC-om)
    ref_fai = ref_genome + ".fai"
    if not os.path.exists(ref_fai):
        print("[POKRETANJE] Pravljenje FASTA indeksa (.fai) preko Samtools-a...")
        run_command(f"samtools faidx {ref_genome}", os.path.join(output_dir, "logs", "samtools_faidx.log"))
    
    # 1b. Pravljenje sequence dictionary (.dict) ako nedostaje
    ref_dict = ref_genome.replace(".fa", ".dict").replace(".fasta", ".dict")
    if not os.path.exists(ref_dict):
        print("[POKRETANJE] Pravljenje sequence dictionary za GATK...")
        run_command(f"./gatk-4.5.0.0/gatk CreateSequenceDictionary -R {ref_genome}", os.path.join(output_dir, "logs", "gatk_dict.log"))

    # 1c. Indeksiranje baze poznatih varijanti ako nedostaje indeks fajl
    known_idx_csi = known_sites + ".csi"
    known_idx_tbi = known_sites + ".tbi"
    if not os.path.exists(known_idx_csi) and not os.path.exists(known_idx_tbi):
        print("[POKRETANJE] Indeksiranje fajla sa poznatim varijantama...")
        run_command(f"./gatk-4.5.0.0/gatk IndexFeatureFile -I {known_sites}", os.path.join(output_dir, "logs", "gatk_index_bcf.log"))

    if check_step_done(output_dir, "processing"):
        print("[PRESKOČENO] GATK obrada (BQSR) je već uspešno završena. Idemo dalje!")
        return

    bam_in = os.path.join(output_dir, "02_alignment", "sample_0.sorted.bam")
    dedup_bam = os.path.join(output_dir, "03_processing", "sample_0.dedup.bam")
    metrics = os.path.join(output_dir, "03_processing", "sample_0.metrics.txt")
    recal_table = os.path.join(output_dir, "03_processing", "sample_0.recal_data.table")
    bqsr_bam = os.path.join(output_dir, "03_processing", "sample_0.bqsr.bam")
    
    tmp_dedup = dedup_bam + ".tmp"
    tmp_bqsr = bqsr_bam + ".tmp"

    # Priprema interval argumenta
    interval_arg = f"-L {interval}" if interval else ""

    try:
        # 2. MarkDuplicates (Ostaje isto)
        if not os.path.exists(dedup_bam):
            print("[POKRETANJE] Uklanjanje PCR duplikata (MarkDuplicates)...")
            cmd_md = f"./gatk-4.5.0.0/gatk MarkDuplicates -I {bam_in} -O {tmp_dedup} -M {metrics}"
            run_command(cmd_md, os.path.join(output_dir, "logs", "gatk_markdup.log"))
            run_command(f"samtools index {tmp_dedup}", os.path.join(output_dir, "logs", "samtools_index_dedup.log"))
            os.rename(tmp_dedup, dedup_bam)
            os.rename(tmp_dedup + ".bai", dedup_bam + ".bai")
        else:
            print("[PRESKOČENO] MarkDuplicates je već urađen.")

        # 3. BaseRecalibrator (OVDE DODAJEMO INTERVAL)
        print("[POKRETANJE] GATK BaseRecalibrator (Računanje tabele grešaka)...")
        cmd_br = f"./gatk-4.5.0.0/gatk BaseRecalibrator -I {dedup_bam} -R {ref_genome} --known-sites {known_sites} {interval_arg} -O {recal_table}"
        run_command(cmd_br, os.path.join(output_dir, "logs", "gatk_baserecal.log"))

        # 4. ApplyBQSR (I ovde dodajemo interval kako bi BAM bio manji i brži za dalji rad)
        print("[POKRETANJE] GATK ApplyBQSR (Primena rekalibracije na BAM)...")
        cmd_apply = f"./gatk-4.5.0.0/gatk ApplyBQSR -I {dedup_bam} -R {ref_genome} --bqsr-recal-file {recal_table} {interval_arg} -O {tmp_bqsr}"
        run_command(cmd_apply, os.path.join(output_dir, "logs", "gatk_applybqsr.log"))
        run_command(f"samtools index {tmp_bqsr}", os.path.join(output_dir, "logs", "samtools_index_bqsr.log"))
        os.rename(tmp_bqsr, bqsr_bam)
        os.rename(tmp_bqsr + ".bai", bqsr_bam + ".bai")

        mark_step_done(output_dir, "processing")
        print("[USPEH] GATK obrada uspešno završena. Fajlovi su u folderu 03_processing.")
    except subprocess.CalledProcessError:
        print("[GREŠKA] GATK obrada je prekinuta! Čistim privremene fajlove...")
        clean_failed_outputs([tmp_dedup, tmp_dedup + ".bai", tmp_bqsr, tmp_bqsr + ".bai"])
        exit(1)

# --- KORACI 4, 5, 6, 7: VARIANT CALLING & HARD FILTERING ---
def run_variant_calling(ref_genome, output_dir, interval=None):
    print("\n--- Korak 4: Pozivanje i filtriranje varijanti (GATK) ---")
    if check_step_done(output_dir, "variant_calling"):
        print("[PRESKOČENO] Pozivanje varijanti je već završeno. Idemo dalje!")
        return

    bqsr_bam = os.path.join(output_dir, "03_processing", "sample_0.bqsr.bam")
    raw_vcf = os.path.join(output_dir, "04_variant_calling", "sample_0.raw.vcf.gz")
    filtered_vcf = os.path.join(output_dir, "04_variant_calling", "sample_0.filtered.vcf.gz")
    pass_vcf = os.path.join(output_dir, "04_variant_calling", "sample_0.pass.vcf.gz")

    tmp_raw = raw_vcf + ".tmp.vcf.gz"
    tmp_filtered = filtered_vcf + ".tmp.vcf.gz"
    tmp_pass = pass_vcf + ".tmp.vcf.gz"

    # Dodavanje opcionog intervala (-L) za ubrzanje pretrage celog hg38 genoma
    interval_arg = f"-L {interval}" if interval else ""

    try:
        # 1. HaplotypeCaller (Pravljenje sirovog VCF-a)
        print("[POKRETANJE] GATK HaplotypeCaller...")
        cmd_hc = f"./gatk-4.5.0.0/gatk HaplotypeCaller -R {ref_genome} -I {bqsr_bam} -O {tmp_raw} {interval_arg}"
        run_command(cmd_hc, os.path.join(output_dir, "logs", "gatk_hc.log"))
        os.rename(tmp_raw, raw_vcf)
        os.rename(tmp_raw + ".tbi", raw_vcf + ".tbi")

        # 2. Hard Filtering (VariantFiltration prema standardnim GATK preporukama)
        print("[POKRETANJE] GATK VariantFiltration (Primenjujem Hard Filtering)...")
        cmd_vf = (
            f"./gatk-4.5.0.0/gatk VariantFiltration -R {ref_genome} -V {raw_vcf} -O {tmp_filtered} "
            f'--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" '
            f'--filter-name "GATK_Hard_Filter"'
        )
        run_command(cmd_vf, os.path.join(output_dir, "logs", "gatk_filter.log"))
        os.rename(tmp_filtered, filtered_vcf)
        os.rename(tmp_filtered + ".tbi", filtered_vcf + ".tbi")

        # 3. SelectVariants (Izdvajanje isključivo PASS varijanti)
        print("[POKRETANJE] GATK SelectVariants (Izdvajam samo PASS varijante)...")
        cmd_pass = f"./gatk-4.5.0.0/gatk SelectVariants -R {ref_genome} -V {filtered_vcf} --exclude-filtered -O {tmp_pass}"
        run_command(cmd_pass, os.path.join(output_dir, "logs", "gatk_select_pass.log"))
        os.rename(tmp_pass, pass_vcf)
        os.rename(tmp_pass + ".tbi", pass_vcf + ".tbi")

        mark_step_done(output_dir, "variant_calling")
        print("[USPEH] Pozivanje i filtriranje varijanti uspešno završeno!")
    except subprocess.CalledProcessError:
        print("[GREŠKA] Variant calling je prekinut! Čistim privremene fajlove...")
        clean_failed_outputs([tmp_raw, tmp_raw + ".tbi", tmp_filtered, tmp_filtered + ".tbi", tmp_pass, tmp_pass + ".tbi"])
        exit(1)

# --- KORAK 8: GENERISANJE MANIFEST FAJLA ---
def generate_manifest(output_dir):
    manifest_path = os.path.join(output_dir, "results_manifest.json")
    manifest = {
        "pipeline_name": "BWA-GATK Variant Calling Pipeline",
        "output_directory": output_dir,
        "results": {
            "fastqc_html_1": os.path.join(output_dir, "01_fastqc", "sample_0.chrom11.exome.pe1_fastqc.html"),
            "fastqc_html_2": os.path.join(output_dir, "01_fastqc", "sample_0.chrom11.exome.pe2_fastqc.html"),
            "sorted_bam": os.path.join(output_dir, "02_alignment", "sample_0.sorted.bam"),
            "dedup_bam": os.path.join(output_dir, "03_processing", "sample_0.dedup.bam"),
            "bqsr_bam": os.path.join(output_dir, "03_processing", "sample_0.bqsr.bam"),
            "raw_vcf": os.path.join(output_dir, "04_variant_calling", "sample_0.raw.vcf.gz"),
            "filtered_vcf": os.path.join(output_dir, "04_variant_calling", "sample_0.filtered.vcf.gz"),
            "pass_vcf": os.path.join(output_dir, "04_variant_calling", "sample_0.pass.vcf.gz")
        }
    }
    with open(manifest_path, "w") as jf:
        json.dump(manifest, jf, indent=4)
    print(f"\n[MANIFEST] Uspešno zapisan u: {manifest_path}")

def main():
    parser = argparse.ArgumentParser(description="Automatizovani BWA-GATK Pipeline")
    parser.add_argument("--fastq1", required=True, help="Putanja do R1 FASTQ fajla")
    parser.add_argument("--fastq2", required=True, help="Putanja do R2 FASTQ fajla")
    parser.add_argument("--reference", required=True, help="Putanja do hg38.fa")
    parser.add_argument("--known-sites", required=True, help="Putanja do baze poznatih varijanti (.vcf.gz)")
    parser.add_argument("--output-dir", required=True, help="Izlazni direktorijum za rezultate")
    parser.add_argument("--interval", required=False, default=None, help="Opcioni genomski region za restrikciju analize (npr. chr11)")

    args = parser.parse_args()

    print("=== Pokretanje kompletnog bioinformatičkog pipeline-a ===")
    
    # Kreiranje foldera na osnovu unetog --output-dir argumenta
    setup_directories(args.output_dir)
    
    # Pokretanje svih koraka redom
    run_fastqc(args.fastq1, args.fastq2, args.output_dir)
    run_bwa(args.fastq1, args.fastq2, args.reference, args.output_dir)
    run_processing(args.reference, args.known_sites, args.output_dir, args.interval)
    run_variant_calling(args.reference, args.output_dir, args.interval)
    
    # Pravljenje JSON manifesta za profesora
    generate_manifest(args.output_dir)
    
    print("\n[INFO] Pipeline je uspešno završio sve zadate faze rada!")

if __name__ == "__main__":
    main()