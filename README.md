# BWA-GATK Variant Calling Pipeline

Ovaj repozitorijum sadrži automatizovani Python pipeline za obradu uparenih FASTQ fajlova i identifikaciju varijanti prema GATK najboljim praksama (Best Practices). Pipeline je prilagođen za analizu egzomskih podataka (hromozom 11) i podržava prekid/nastavak rada (idempotentnost).

## Korišćeni alati i verzije
Za uspešno pokretanje ovog pipeline-a, neophodno je da u okruženju budu instalirani sledeći alati:
* **FastQC** (za kontrolu kvaliteta sekvenci)
* **BWA / BWA-MEM** (za mapiranje na referentni genom)
* **Samtools** (za manipulaciju SAM/BAM fajlovima i indeksiranje)
* **GATK** (verzija 4.5.0.0 - za procesiranje BAM fajlova, rekalibraciju i pozivanje varijanti)
* **Python 3** (sa bibliotekama `pysam` i `matplotlib` za dodatnu analitiku i crtanje grafika)

## Uputstvo za pokretanje (Usage)

Glavna skripta `main.py` prihvata putanje do ulaznih fajlova, baze poznatih varijanti i izlaznog direktorijuma. Opciono, analiza se može ograničiti na specifičan genomski region (npr. `chr11`).

**Primer komande za pokretanje kompletnog pipeline-a:**

```bash
python main.py \
  --fastq1 sample_0.chrom11.exome.pe1.fq.gz \
  --fastq2 sample_0.chrom11.exome.pe2.fq.gz \
  --reference reference/hg38.fa \
  --known-sites known_sites_chr11_fixed.vcf.gz \
  --output-dir results \
  --interval chr11