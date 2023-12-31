
# awk -f chrom-sizes.awk Results/07-reference/chromosome-sizes.tsv

BEGIN {
    chrom_total=0;
    contig_total=0;
}

$1 ~ "chr" {
    chrom_total += $2;
    print;
}

$1 ~ "^C" {
    contig_total += $2;
}

END {
    printf("Chromosome bases: %d  Contig bases: %d\n", chrom_total, contig_total);
    printf("Contigs = %f%%\n", contig_total / (chrom_total + contig_total) * 100.0);
}
