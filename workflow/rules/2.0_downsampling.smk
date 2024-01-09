localrules:
    downsample_bamlist,


rule downsample_bamlist:
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}.bamlist",
        bams=angsd.get_bamlist_bams,
        bais=angsd.get_bamlist_bais,
    output:
        sublist="results/datasets/{dataset}/bamlists/downsampled/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.bamlist",
    log:
        "logs/{dataset}/downsampled/sublist/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/downsampled/sublist/{dataset}.{ref}_{population}{dp}.N{samplesize}-rep{rep}_{sites}-filts.log"
    conda:
        "../envs/shell.yaml"
    wildcard_constraints:
        dp=".{0}",
    resources:
        runtime="10m",
    shell:
        """
        shuf -n{wildcards.samplesize} {input.bamlist} > {output.sublist} 2> {log}
        """


rule downsample_thetas:
    input:
        sublist="results/datasets/{dataset}/bamlists/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.bamlist",
        ref="results/ref/{ref}/{ref}.fa",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
    output:
        sfs=ensure(
            "results/datasets/{dataset}/analyses/sfs/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.sfs",
            non_empty=True,
        ),
        thetas="results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.pestPG",
    container:
        angsd.angsd_container
    log:
        "logs/{dataset}/downsampled/thetas/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.log",
    group:
        "theta"
    benchmark:
        "benchmarks/{dataset}/downsampled/thetas/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.log"
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        fold=config["params"]["realsfs"]["fold"],
        out=lambda w, output: os.path.splitext(output.thetas)[0],
    threads: lambda w, attempt: attempt * 2
    resources:
        runtime=lambda w, attempt: attempt * 1200,
    shell:
        """
        (angsd -doSaf 1 -bam {input.sublist} -GL {params.gl_model} -ref {input.ref} \
            -nThreads {threads} {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -anc {input.ref} \
            -rf <(cut -f1 {input.sites} | uniq | sed -e 's/$/:/') \
            -out {resources.tmpdir}/saf

        realSFS {resources.tmpdir}/saf.saf.idx -fold {params.fold} \
            -P {threads} > {output.sfs}

        realSFS saf2theta {resources.tmpdir}/saf.saf.idx -sfs {output.sfs} \
            -fold {params.fold} -outname {resources.tmpdir}/thetas
        
        thetaStat do_stat {resources.tmpdir}/thetas.thetas.idx \
            -win {wildcards.win} -type 0 -step {wildcards.step} \
            -outnames {params.out}) &> {log}
        """


rule average_downsampled_thetas:
    input:
        "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaWindows.{win}_{step}.pestPG",
    output:
        ensure(
            "results/datasets/{dataset}/analyses/thetas/downsampled/{dataset}.{ref}_{population}.N{samplesize}-rep{rep}_{sites}-filts.thetaMean.{win}_{step}.tsv"
        ),
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
    resources:
        runtime=lambda w, attempt: attempt * 360,
    shell:
        """
        (printf "pop\tdownsample.n\tdownsample.rep\tpi\twatterson\ttajima\n" > {output}
        cat {input} >> {output}) 2> {log}
        """


rule downsample_fst:
    input:
        sublist1="results/datasets/{dataset}/bamlists/downsampled/{dataset}.{ref}_{population1}.N{samplesize}-rep{rep}_{sites}-filts.bamlist",
        sublist2="results/datasets/{dataset}/bamlists/downsampled/{dataset}.{ref}_{population2}.N{samplesize}-rep{rep}_{sites}-filts.bamlist",
        ref="results/ref/{ref}/{ref}.fa",
        sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites",
        idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}_{sites}-filts.sites.idx",
    output:
        sfs=temp(
            ensure(
                "results/datasets/{dataset}/analyses/sfs/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.sfs",
                non_empty=True,
            )
        ),
        fstglob=ensure(
            "results/datasets/{dataset}/analyses/fst/downsampled/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.fst.global.tsv"
        ),
    container:
        angsd.angsd_container
    log:
        "logs/{dataset}/downsampled/fst/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/downsampled/fst/{dataset}.{ref}_{population1}-{population2}.N{samplesize}-rep{rep}_{sites}-filts.log"
    wildcard_constraints:
        population1="|".join([i for i in angsd.samples.population.values.tolist()]),
        population2="|".join([i for i in angsd.samples.population.values.tolist()]),
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        fold=config["params"]["realsfs"]["fold"],
        fst=config["params"]["fst"]["whichFst"],
        full_n1=lambda w: len(angsd.get_samples_from_pop(w.population1)),
        full_n2=lambda w: len(angsd.get_samples_from_pop(w.population2)),
    threads: lambda w, attempt: attempt * 3
    resources:
        runtime=lambda w, attempt: attempt * 1320,
    shell:
        """
        (angsd -doSaf 1 -bam {input.sublist1} -GL {params.gl_model} -ref {input.ref} \
            -nThreads 1 {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -anc {input.ref} \
            -rf <(cut -f1 {input.sites} | uniq | sed -e 's/$/:/') \
            -out {resources.tmpdir}/saf1
        
        angsd -doSaf 1 -bam {input.sublist2} -GL {params.gl_model} -ref {input.ref} \
            -nThreads 1 {params.extra} -minMapQ {params.mapQ} \
            -minQ {params.baseQ} -sites {input.sites} -anc {input.ref} \
            -rf <(cut -f1 {input.sites} | uniq | sed -e 's/$/:/') \
            -out {resources.tmpdir}/saf2

        realSFS \
            {resources.tmpdir}/saf1.saf.idx \
            {resources.tmpdir}/saf2.saf.idx \
            -fold {params.fold} -P {threads} > {output.sfs}
        
        realSFS fst index -whichFst {params.fst} \
            {resources.tmpdir}/saf1.saf.idx \
            {resources.tmpdir}/saf2.saf.idx \
            -sfs {output.sfs} \
            -fstout {resources.tmpdir}/fst

        realSFS fst stats {resources.tmpdir}/fst.fst.idx | \
            awk '{{print "{wildcards.population1}\t{wildcards.population2}\t"\
            $1"\t"$2"\t"{wildcards.samplesize}"\t"{wildcards.rep}"\t"\
            {params.full_n1}"\t"{params.full_n2}}}' \
            > {output.fstglob}) &> {log}
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
    resources:
        runtime=lambda w, attempt: attempt * 360,
    shell:
        """
        (printf "pop1\tpop2\tunweight.fst\tweight.fst\tdownsample.n\tdownsample.rep\tpop1.full.samplesize\tpop2.full.samplesize\n" > {output}
        cat {input} >> {output}) 2> {log}
        """
