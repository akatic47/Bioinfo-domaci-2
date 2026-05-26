import pysam
import matplotlib.pyplot as plt

def get_vcf_stats(vcf_path):
    """Parsira VCF fajl i vraća osnovnu statistiku."""
    vcf = pysam.VariantFile(vcf_path)
    
    stats = {
        'total': 0,
        'snp': 0,
        'indel': 0
    }
    
    for rec in vcf:
        stats['total'] += 1
        # Jednostavna provera da li je SNP (sve alele su dužine 1)
        is_snp = all(len(allele) == 1 for allele in rec.alleles)
        if is_snp:
            stats['snp'] += 1
        else:
            stats['indel'] += 1
            
    return stats

def main():
    # Putanje do fajlova (prilagodite ako je potrebno)
    raw_vcf_path = "results/04_variant_calling/sample_0.raw.vcf.gz"
    pass_vcf_path = "results/04_variant_calling/sample_0.pass.vcf.gz"
    
    # 1. Prikupljanje statistike
    raw_stats = get_vcf_stats(raw_vcf_path)
    pass_stats = get_vcf_stats(pass_vcf_path)
    
    # Računanje dodatnih metrika
    failed_count = raw_stats['total'] - pass_stats['total']
    failed_percent = (failed_count / raw_stats['total']) * 100 if raw_stats['total'] > 0 else 0

    # 2. Ispis rezultata u terminal (korisno za izveštaj)
    print("=== Statistika varijanti ===")
    print(f"Ukupan broj varijanti (pre filtriranja): {raw_stats['total']}")
    print(f"Ukupan broj varijanti (PASS): {pass_stats['total']}")
    print(f"Broj varijanti koje nisu prošle filter: {failed_count}")
    print(f"Procenat odbačenih varijanti: {failed_percent:.2f}%")
    print("-" * 30)
    print(f"SNP pre filtriranja: {raw_stats['snp']} | SNP posle filtriranja (PASS): {pass_stats['snp']}")
    print(f"INDEL pre filtriranja: {raw_stats['indel']} | INDEL posle filtriranja (PASS): {pass_stats['indel']}")

    # 3. Pravljenje grafikona (Bar plot)
    categories = ['Total', 'SNP', 'INDEL']
    raw_counts = [raw_stats['total'], raw_stats['snp'], raw_stats['indel']]
    pass_counts = [pass_stats['total'], pass_stats['snp'], pass_stats['indel']]

    x = range(len(categories))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar([i - width/2 for i in x], raw_counts, width, label='Sirove (Raw)', color='skyblue')
    ax.bar([i + width/2 for i in x], pass_counts, width, label='Filtrirane (PASS)', color='lightgreen')

    ax.set_ylabel('Broj varijanti')
    ax.set_title('Poređenje varijanti pre i posle hard filtering-a')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend()

    # Čuvanje grafika kao slike (možete je posle ubaciti u PDF izveštaj)
    plt.savefig('variant_comparison.png')
    print("\nGrafikon je sačuvan kao 'variant_comparison.png'.")

if __name__ == "__main__":
    main()