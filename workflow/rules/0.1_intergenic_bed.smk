localrules:
    gff2intergenic,
    user_sites,


rule gff2intergenic:
    input:
        gff=config["gff"],
        genbed="results/ref/{ref}/beds/genome.bed",
        fai="results/ref/{ref}/{ref}.fa.fai",
    output:
        genes="results/datasets/{dataset}/filters/intergenic/{ref}_genes.gff3.gz",
        bed="results/datasets/{dataset}/filters/intergenic/{ref}_intergenic.bed",
    benchmark:
        "benchmarks/{dataset}/filters/gff2intergenic/{dataset}.{ref}_gff2intergenic.log"
    log:
        "logs/{dataset}/filters/gff2intergenic/{dataset}.{ref}_gff2intergenic.log",
    params:
        buffer=config["gene-buffer"],
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (zcat {input.gff} | grep "\tgene\t" | gzip >> {output.genes}

        bedtools slop -i {output.genes} -g <(cut -f1-2 {input.fai}) -b {params.buffer} \
            | bedtools subtract -a {input.genbed} -b - > {output.bed}) 2> {log}
        """
