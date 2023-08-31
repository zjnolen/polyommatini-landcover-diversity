rule downsample_doSaf:
    """
    Generate a site allele frequency file for a given downsampled population and genome
    chunk.
    """
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}.bamlist",
        ref="results/ref/{ref}/{ref}.fa",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
    output:
        sublist="results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.bamlist",
        saf=temp(
            "results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.saf.gz"
        ),
        safidx=temp(
            "results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.saf.idx"
        ),
        safpos=temp(
            "results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.saf.pos.gz"
        ),
        arg="results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.log"
    container:
        angsd.angsd_container
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    resources:
        runtime="120m",
    shell:
        """
        (shuf -n{wildcards.samplesize} {input.bamlist} > {output.sublist}

        angsd -doSaf 1 -bam {output.sublist} -GL {params.gl_model} -ref {input.ref} \
            -nThreads {threads} {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -anc {input.ref} \
            -out {params.out}) &> {log}
        """


rule downsample_1dSFS:
    """
    Generate a downsampled 1D site frequency spectrum.
    """
    input:
        saf="results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.saf.idx",
        others=multiext(
            "results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
    output:
        sfs=temp(
            ensure(
                "results/datasets/{dataset}/analyses/sfs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.sfs",
                non_empty=True,
            )
        ),
    container:
        angsd.angsd_container
    log:
        "logs/{dataset}/realSFS/1dSFS/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.log",
    group:
        "theta"
    benchmark:
        "benchmarks/{dataset}/realSFS/1dSFS/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.log"
    params:
        fold=config["params"]["realsfs"]["fold"],
    threads: 2
    resources:
        runtime="120m",
    shell:
        """
        realSFS {input.saf} -fold {params.fold} -P {threads} \
            > {output.sfs} 2> {log}
        """


rule downsample_2dSFS:
    """
    Generate a 2D site frequency spectrum.
    """
    input:
        saf1="results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population1}.N{samplesize}-rep{rep}_{sites}-filts.saf.idx",
        saf1_others=multiext(
            "results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population1}.N{samplesize}-rep{rep}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
        saf2="results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population2}.N{samplesize}-rep{rep}_{sites}-filts.saf.idx",
        saf2_others=multiext(
            "results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population2}.N{samplesize}-rep{rep}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
    output:
        sfs=temp(
            ensure(
                "results/datasets/{dataset}/analyses/sfs/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.sfs",
                non_empty=True,
            )
        ),
    container:
        angsd.angsd_container
    group:
        "fst"
    log:
        "logs/{dataset}/realSFS/2dSFS/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/2dSFS/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.log"
    wildcard_constraints:
        population1="|".join(
            [i for i in angsd.samples.index.tolist()]
            + [i for i in angsd.samples.population.values.tolist()]
        ),
        population2="|".join(
            [i for i in angsd.samples.index.tolist()]
            + [i for i in angsd.samples.population.values.tolist()]
        ),
    params:
        fold=config["params"]["realsfs"]["fold"],
    threads: 2
    resources:
        runtime="180m",
    shell:
        """
        realSFS {input.saf1} {input.saf2} -fold {params.fold} \
            -P {threads} > {output.sfs} 2> {log}
        """


rule downsample_saf2theta:
    """
    Generates a thetas index for each population from SAF and SFS.
    """
    input:
        safidx="results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.saf.idx",
        saf_others=multiext(
            "results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
        sfs="results/datasets/{dataset}/analyses/sfs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.sfs",
    output:
        thetasidx=temp(
            "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetas.idx"
        ),
        thetas=temp(
            "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetas.gz"
        ),
    container:
        angsd.angsd_container
    group:
        "theta"
    log:
        "logs/{dataset}/realSFS/saf2theta/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.log",
    params:
        out=lambda w, output: output.thetas.removesuffix(".thetas.gz"),
        fold=config["params"]["realsfs"]["fold"],
    resources:
        runtime="120m",
    threads: 2
    shell:
        """
        realSFS saf2theta {input.safidx} -sfs {input.sfs} -fold {params.fold} \
            -outname {params.out} &> {log}
        """


rule downsample_thetaStat:
    """
    Estimates diversity and neutrality statistics in sliding windows of a set size.
    """
    input:
        thetasidx="results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetas.idx",
        thetas="results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetas.gz",
    output:
        thetas=temp(
            "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.pestPG"
        ),
    container:
        angsd.angsd_container
    group:
        "theta"
    log:
        "logs/{dataset}/thetaStat/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.{win}_{step}.log",
    params:
        out=lambda w, output: os.path.splitext(output.thetas)[0],
    threads: 2
    resources:
        runtime="120m",
    shell:
        """
        thetaStat do_stat {input.thetasidx} -win {wildcards.win} -type 0 \
            -step {wildcards.step} -outnames {params.out} &> {log}
        """


rule average_downsampled_thetas:
    input:
        "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.pestPG",
    output:
        "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaMean.{win}_{step}.tsv",
    conda:
        "../envs/r.yaml"
    group:
        "theta"
    log:
        "logs/{dataset}/thetaStat/average/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaMean.{win}_{step}.log",
    params:
        minsites=config["params"]["thetas"]["minsites"],
    script:
        "../scripts/average_downsampled_thetas.R"


rule aggregate_downsampled_thetas:
    input:
        expand(
            "results/datasets/{{dataset}}/analyses/thetas/downsampled/{{dataset}}.{{ref}}_{population}.N{samplesize}-rep{rep}_{{sites}}-filts.thetaMean.{{win}}_{{step}}.tsv",
            population=angsd.pop_list,
            samplesize=config["downsample_sizes"],
            rep=list(range(0, config["downsample_reps"][0])),
        ),
    output:
        "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_all.downsampled_{sites}-filts.thetaMean.{win}_{step}.tsv",
    log:
        "logs/{dataset}/thetaStat/aggregate/downsampled/{dataset}.{ref}_all.downsampled_{sites}-filts.thetaMean.{win}_{step}.log",
    conda:
        "../envs/shell.yaml"
    shell:
        """
        (printf "pop\tdownsample.n\tdownsample.rep\tpi\twatterson\ttajima\n" > {output}
        cat {input} >> {output}) 2> {log}
        """


use rule realSFS_fst_index from angsd as downsample_fst_index with:
    input:
        saf1="results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population1}.N{samplesize}-rep{rep}_{sites}-filts.saf.idx",
        saf1_others=multiext(
            "results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population1}.N{samplesize}-rep{rep}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
        saf2="results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population2}.N{samplesize}-rep{rep}_{sites}-filts.saf.idx",
        saf2_others=multiext(
            "results/datasets/{dataset}/safs/downsampled/{dataset}.{ref}_{population2}.N{samplesize}-rep{rep}_{sites}-filts.saf",
            ".pos.gz",
            ".gz",
        ),
        sfs="results/datasets/{dataset}/analyses/sfs/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.sfs",
    output:
        fstidx=temp(
            "results/datasets/{dataset}/analyses/fst/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.fst.idx"
        ),
    group:
        "fst"
    log:
        "logs/{dataset}/realSFS/fst/index/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/fst/index/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.log"
    threads: 2
    resources:
        runtime="120m",


rule downsample_fst_stats:
    input:
        fstidx="results/datasets/{dataset}/analyses/fst/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.fst.idx",
    output:
        fstglob="results/datasets/{dataset}/analyses/fst/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.fst.global.tsv",
    container:
        angsd.angsd_container
    group:
        "fst"
    log:
        "logs/{dataset}/realSFS/fst/stats/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/fst/stats/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.log"
    threads: 2
    resources:
        runtime="60m",
    params:
        full_n1=lambda w: len(angsd.get_samples_from_pop(w.population1)),
        full_n2=lambda w: len(angsd.get_samples_from_pop(w.population2)),
    shell:
        r"""
        realSFS fst stats {input.fstidx} | \
            awk '{{print "{wildcards.population1}\t{wildcards.population2}\t"\
            $1"\t"$2"\t"{wildcards.samplesize}"\t"{wildcards.rep}"\t"\
            {params.full_n1}"\t"{params.full_n2}}}' \
            > {output.fstglob} 2> {log}
        """


def get_downsample_fst(wildcards):
    unit = angsd.pop_list
    combos = angsd.pairwise_combos(unit)
    pop1 = [pair[0] for pair in combos]
    pop2 = [pair[1] for pair in combos]
    return expand(
        expand(
            "results/datasets/{{dataset}}/analyses/fst/downsampled/{{dataset}}.{{ref}}_{population1}-{population2}.N{{samplesize}}-rep{{rep}}_{{sites}}-filts.fst.global.tsv",
            zip,
            population1=pop1,
            population2=pop2,
        ),
        samplesize=config["downsample_sizes"],
        rep=list(range(0, config["downsample_reps"][0])),
        allow_missing=True,
    )


rule aggregate_downsample_fst:
    input:
        get_downsample_fst,
    output:
        "results/datasets/{dataset}/analyses/fst/downsampled/{dataset}.{ref}_poppairs.downsampled_{sites}-filts.fst.global.tsv",
    log:
        "logs/{dataset}/realSFS/fst/aggregate/downsampled/{dataset}.{ref}_poppairs_{sites}-filts.global.log",
    benchmark:
        "benchmarks/{dataset}/realSFS/fst/aggregate/downsampled/{dataset}.{ref}_poppairs_{sites}-filts.global.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        (printf "pop1\tpop2\tunweight.fst\tweight.fst\tdownsample.n\tdownsample.rep\tpop1.full.samplesize\tpop2.full.samplesize\n" > {output}
        cat {input} >> {output}) 2> {log}
        """
