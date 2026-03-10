process SUMMARY {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path falco_pre_txt
    path falco_post_txt
    path fastp_json
    path star_logs
    path fc_summaries

    output:
    path "summary.tsv"

    script:
    """
    #!/usr/bin/env python3
    import glob
    
    data = {}
    samples = set()
    TAB = chr(9)
    NL = chr(10)
    
    def add_metric(sample, metric, value):
        if sample not in data:
            data[sample] = {}
        data[sample][metric] = value
        samples.add(sample)
    
    for f in glob.glob("*_R1_pre_data.txt"):
        sample = f.replace("_R1_pre_data.txt", "")
        with open(f) as fh:
            for line in fh:
                if line.startswith("Total Sequences"):
                    add_metric(sample, "1_Falco_Pre_Total_Seq", line.strip().split(TAB)[1])
                elif line.startswith("%GC"):
                    add_metric(sample, "1_Falco_Pre_GC_Percent", line.strip().split(TAB)[1])
    
    for f in glob.glob("*_R1_post_data.txt"):
        sample = f.replace("_R1_post_data.txt", "")
        with open(f) as fh:
            for line in fh:
                if line.startswith("Total Sequences"):
                    add_metric(sample, "2_Falco_Post_Total_Seq", line.strip().split(TAB)[1])
                elif line.startswith("%GC"):
                    add_metric(sample, "2_Falco_Post_GC_Percent", line.strip().split(TAB)[1])
    
    for f in glob.glob("*_Log.final.out"):
        sample = f.replace("_Log.final.out", "")
        with open(f) as fh:
            for line in fh:
                if "Number of input reads" in line:
                    add_metric(sample, "3_STAR_Input_Reads", line.split("|")[1].strip())
                elif "Uniquely mapped reads number" in line:
                    add_metric(sample, "3_STAR_Uniquely_Mapped", line.split("|")[1].strip())
                elif "Uniquely mapped reads %" in line:
                    add_metric(sample, "3_STAR_Uniquely_Mapped_Percent", line.split("|")[1].strip())
    
    for f in glob.glob("*_featurecounts.txt.summary"):
        sample = f.replace("_featurecounts.txt.summary", "")
        with open(f) as fh:
            for i, line in enumerate(fh):
                if i == 0: continue
                parts = line.strip().split(TAB)
                if len(parts) >= 2:
                    add_metric(sample, "4_FC_" + parts[0], parts[1])
    
    with open("summary.tsv", "w") as out:
        samples = sorted(list(samples))
        if not samples:
             out.write("No metrics gathered\\n")
             exit(0)
        out.write("Metric" + TAB + TAB.join(samples) + NL)
        seen_metrics = set()
        for s in data:
            seen_metrics.update(data[s].keys())
        for m in sorted(list(seen_metrics)):
            row = [m] + [str(data.get(s, {}).get(m, "NA")) for s in samples]
            out.write(TAB.join(row) + NL)
    """
}
