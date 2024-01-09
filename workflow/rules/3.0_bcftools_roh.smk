localrules:
    combine_roh,
    gatk_list,
    bcftools_list,
    bcf2plink,


rule angsd_doBcf_likes:
    input:
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=angsd.get_bamlist_bams,
        bais=angsd.get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
    output:
        bcf=temp(
            "results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.bcf"
        ),
        maf=temp(
            "results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.mafs.gz"
        ),
        arg="results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.arg",
    log:
        "logs/{dataset}/angsd/doBcf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doBcf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.GLonly.log"
    container:
        angsd.angsd_container
    wildcard_constraints:
        population="|".join([i for i in angsd.samples.population.values.tolist()]),
        dp=".{0}",
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        snp_pval=config["params"]["angsd"]["snp_pval"],
        minmaf=config["params"]["angsd"]["min_maf"],
        minind=lambda w: math.ceil(len(angsd.get_samples_from_pop(w.population)) / 2),
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        """
        angsd -doBcf 1 -bam {input.bam} -GL {params.gl_model} -ref {input.ref} \
            -doMajorMinor 1 -doMaf 1 -SNP_pval {params.snp_pval} -minInd {params.minind} \
            -minMaf {params.minmaf} -nThreads {threads} {params.extra} \
            -minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
            -dopost 1 --ignore-RG 0 -dogeno 1 -docounts 1 -rf {input.regions} \
            -setMinDepthInd 10 -out {params.out} &> {log}
        """


rule concat_bcf:
    input:
        calls=lambda w: expand(
            "results/datasets/{{dataset}}/bcfs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_{{sites}}-filts.GLonly.bcf",
            chunk=angsd.chunklist,
        ),
    output:
        bcf="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.GLonly.bcf",
        pos="results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.GLonly.pos",
    log:
        "logs/{dataset}/bcftools/concat/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-bcf.GLonly.log",
    benchmark:
        "benchmarks/{dataset}/bcftools/concat/{dataset}.{ref}_{population}{dp}_{sites}-filts_merge-bcf.GLonly.log"
    conda:
        "../envs/bcftools.yaml"
    wildcard_constraints:
        population="|".join([i for i in angsd.samples.population.values.tolist()]),
        dp=".{0}",
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        """
        bcftools concat -Ob -o {output.bcf} {input} 2> {log}
        bcftools view {output.bcf} | grep -v "#" | cut -f1,2 > {output.pos} 2>> {log}
        """


rule bcftools_roh_likes:
    input:
        "results/datasets/{dataset}/bcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.GLonly.bcf",
    output:
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.GLonly.roh",
        sites="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.sites.GLonly.roh",
        regs="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.regs.GLonly.roh",
    log:
        "logs/{dataset}/bcftools/roh/{dataset}.{ref}_{population}{dp}_{sites}-filts_roh.GLonly.log",
    wildcard_constraints:
        population="|".join([i for i in angsd.samples.population.values.tolist()]),
        dp=".{0}",
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt
    params:
        recrate=config["recrate"],
    shell:
        """
        (bcftools roh -M {params.recrate} --AF-tag AF -o {output.roh} {input}
        awk '$1=="ST"' {output.roh} > {output.sites}
        awk '$1=="RG"' {output.roh} > {output.regs}) 2> {log}
        """


rule bcftools_list:
    input:
        bed="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.bed",
        reg="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
    output:
        bed="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_chunk{chunk}_{sites}-filts.bed",
    conda:
        "../envs/shell.yaml"
    shell:
        """
        > {output.bed}

        while read chr; do
            awk -v chr=$chr '$1 == chr' {input.bed} >> {output.bed}
        done < {input.reg}
        """


rule bcftools_mpileup_call:
    input:
        alignments=lambda w: expand(
            "results/datasets/{{dataset}}/bams/{sample}.{{ref}}.bam",
            sample=angsd.get_samples_from_pop(w.population),
        ),
        bai=lambda w: expand(
            "results/datasets/{{dataset}}/bams/{sample}.{{ref}}.bam.bai",
            sample=angsd.get_samples_from_pop(w.population),
        ),
        ref="results/ref/{ref}/{ref}.fa",
        index="results/ref/{ref}/{ref}.fa.fai",
        regions="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_chunk{chunk}_{sites}-filts.bed",
    output:
        bcf="results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.bcf",
    conda:
        "../envs/bcftools.yaml"
    params:
        extra=lambda w, input: f"--no-BAQ --min-MQ 30 --min-BQ 30 -a FORMAT/AD,FORMAT/DP,INFO/AD",
    log:
        "logs/{dataset}/bcftools/mpileup/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.log",
    resources:
        runtime=360,
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -R {input.regions} -Ou \
            {params.extra} {input.alignments} |
        bcftools call -m --variants-only --skip-variants indels -f GQ,GP -o {output.bcf}
        """


rule bcf_filter:
    input:
        "results/datasets/{dataset}/bcfs/chunks/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.bcf",
    output:
        "results/datasets/{dataset}/vcfs/chunks/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.filter.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/{dataset}/bcftools/filter/{dataset}.{ref}_{population}_chunk{chunk}_{sites}-filts.log",
    shell:
        """
        bcftools filter -e'QUAL<30' -Ou {input} | bcftools +setGT -Ou -- -t q -n . \
            -i'FMT/DP<10 | FMT/GQ<30' | bcftools +fill-tags -Ou -- -t all | \
            bcftools filter -i 'F_MISSING<0.5' -Oz -o {output} 2> {log}
        """


rule bcf_concat:
    input:
        expand(
            "results/datasets/{{dataset}}/vcfs/chunks/{{dataset}}.{{ref}}_{{population}}_chunk{chunk}_{{sites}}-filts.filter.vcf.gz",
            chunk=angsd.chunklist,
        ),
    output:
        vcf="results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}_{sites}-filts.filter.vcf.gz",
        pos="results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}_{sites}-filts.filter.pos",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/{dataset}/bcftools/concat/{dataset}.{ref}_{population}_{sites}-filts.log",
    shell:
        """
        bcftools concat -Oz -o {output.vcf} {input} 2> {log}
        zcat {output.vcf} | grep -v "#" | cut -f1,2 > {output.pos} 2>> {log}
        """


rule bcftools_roh_calls:
    input:
        "results/datasets/{dataset}/vcfs/{dataset}.{ref}_{population}{dp}_{sites}-filts.filter.vcf.gz",
    output:
        roh="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.calls.roh",
        sites="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.sites.calls.roh",
        regs="results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_{population}{dp}_{sites}-filts.regs.calls.roh",
    log:
        "logs/{dataset}/bcftools/roh/{dataset}.{ref}_{population}{dp}_{sites}-filts_roh.calls.log",
    wildcard_constraints:
        population="|".join([i for i in angsd.samples.population.values.tolist()]),
        dp=".{0}",
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt
    params:
        recrate=config["recrate"],
    shell:
        """
        (bcftools roh -G30 -M {params.recrate} --AF-tag AF -o {output.roh} {input}
        awk '$1=="ST"' {output.roh} > {output.sites}
        awk '$1=="RG"' {output.roh} > {output.regs}) 2> {log}
        """


rule combine_roh:
    input:
        expand(
            "results/datasets/{{dataset}}/analyses/roh/bcftools/{{dataset}}.{{ref}}_{population}{{dp}}_{{sites}}-filts.regs.{{type}}.roh",
            population=angsd.pop_list,
        ),
    output:
        "results/datasets/{dataset}/analyses/roh/bcftools/{dataset}.{ref}_all{dp}_{sites}-filts.regs.{type}.roh",
    log:
        "logs/{dataset}/bcftools/roh/{dataset}.{ref}_all{dp}_{sites}-filts_roh.{type}.log",
    conda:
        "../envs/shell.yaml"
    wildcard_constraints:
        dp=".{0}",
    shell:
        """
        cat {input} > {output} 2> {log}
        """
