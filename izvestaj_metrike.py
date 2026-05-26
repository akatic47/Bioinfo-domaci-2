import subprocess
import pysam
import matplotlib.pyplot as plt

def get_flagstat_metrics(bam_path):
    """Pokreće samtools flagstat i izvlači broj mapiranih i uparenih readova."""
    print(f"Pokrećem samtools flagstat za {bam_path}...")
    result = subprocess.run(["samtools", "flagstat", bam_path], capture_output=True, text=True)
    
    total_reads = 0
    mapped_reads = 0
    properly_paired = 0
    
    for line in result.stdout.split('\n'):
        if "in total" in line:
            total_reads = int(line.split()[0])
        elif "mapped (" in line:
            mapped_reads = int(line.split()[0])
        elif "properly paired (" in line:
            properly_paired = int(line.split()[0])
            
    return total_reads, mapped_reads, properly_paired

def parse_picard_metrics(metrics_file):
    """Čita Picard metrics fajl da bi izvukao procenat duplikata."""
    dup_percent = 0.0
    with open(metrics_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("LIBRARY") and "PERCENT_DUPLICATION" in line:
                # Sledeća linija sadrži vrednosti
                values = lines[i+1].split('\t')
                # Kolona 8 u standardnom MarkDuplicates izlazu je PERCENT_DUPLICATION
                try:
                    dup_percent = float(values[8]) * 100
                except (IndexError, ValueError):
                    pass
                break
    return dup_percent

def plot_insert_size(bam_path, output_png):
    """Izvlači template_length iz BAM fajla i crta histogram."""
    print(f"Čitam BAM fajl za histogram fragmenata: {bam_path}...")
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    insert_sizes = []
    # Čitamo prvih milion readova radi brzine (dovoljno za reprezentativan uzorak)
    for i, read in enumerate(bam.head(1000000)):
        if read.is_proper_pair and read.is_read1 and read.template_length > 0:
            # Uzimamo samo pozitivne vrednosti do 1000bp (da izbegnemo outliere)
            if read.template_length < 1000:
                insert_sizes.append(read.template_length)
    
    print("Crtam histogram...")
    plt.figure(figsize=(8, 5))
    plt.hist(insert_sizes, bins=100, color='coral', edgecolor='black', alpha=0.7)
    plt.title('Histogram dužina sekvenciranih fragmenata (Template Length)')
    plt.xlabel('Dužina fragmenta (bp)')
    plt.ylabel('Broj readova')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.savefig(output_png)
    print(f"Grafikon sačuvan kao: {output_png}")

def main():
    # Putanje do fajlova u tvom workspace-u
    sorted_bam = "results/02_alignment/sample_0.sorted.bam"
    dup_metrics = "results/03_processing/sample_0.metrics.txt"
    dedup_bam = "results/03_processing/sample_0.dedup.bam"
    
    # 1. Statistika mapiranja (BWA korak)
    total, mapped, paired = get_flagstat_metrics(sorted_bam)
    print("\n=== BWA-MEM Statistika ===")
    print(f"Ukupno readova: {total}")
    print(f"Mapirani readovi: {mapped} ({(mapped/total)*100:.2f}%)")
    print(f"Pravilno upareni (properly paired): {paired} ({(paired/total)*100:.2f}%)")
    
    # 2. Statistika duplikata
    dup_pct = parse_picard_metrics(dup_metrics)
    print("\n=== MarkDuplicates Statistika ===")
    print(f"Procenat PCR/optičkih duplikata: {dup_pct:.2f}%")
    
    # 3. Generisanje histograma
    print("\n=== Crtanje histograma ===")
    plot_insert_size(dedup_bam, "template_length_histogram.png")

if __name__ == "__main__":
    main()