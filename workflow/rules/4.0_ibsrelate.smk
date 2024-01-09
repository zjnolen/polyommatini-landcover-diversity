rule angsd_doGlf1:
    """
    Generate a site allele frequency file for a given population and genome chunk.
    """
    input:
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=angsd.get_bamlist_bams,
        bais=angsd.get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
    output:
        glf=temp(
            "results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf.gz"
        ),
        arg="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.arg",
    log:
        "logs/{dataset}/angsd/doGlf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/doGlf1/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd.angsd_container
    wildcard_constraints:
        population="|".join(
            ["all"] + [i for i in angsd.samples.population.values.tolist()]
        ),
        dp=".{0}",
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: lambda wildcards, attempt: attempt * 2
    shell:
        """
        angsd -doGlf 1 -bam {input.bam} -GL {params.gl_model} -ref {input.ref} \
            -nThreads {threads} {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -rf {input.regions} \
            -out {params.out} &> {log}
        """


rule ibsrelate:
    input:
        "results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.glf.gz",
    output:
        "results/datasets/{dataset}/analyses/kinship/ibsrelate/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.ibspair",
    log:
        "logs/{dataset}/angsd/ibsrelate/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/ibsrelate/{dataset}.{ref}_{population}{dp}_chunk{chunk}_{sites}-filts.log"
    container:
        angsd.angsd_container
    wildcard_constraints:
        population="|".join(
            ["all"] + [i for i in angsd.samples.population.values.tolist()]
        ),
        dp=".{0}",
    params:
        nind=lambda w: len(angsd.get_samples_from_pop(w.population)),
        maxsites=config["chunk_size"],
        out=lambda w, output: os.path.splitext(output[0])[0],
    resources:
        runtime="7d",
    threads: lambda wildcards, attempt: attempt * 10
    shell:
        """
        ibs -glf {input} -model 0 -nInd {params.nind} -allpairs 1 \
            -m {params.maxsites} -outFileName {params.out}
        """
