function isforwardstrand(record::BAM.Record)
    if BAM.flag(record) & SAM.FLAG_READ1 ≠ 0               # if it's read1
        if BAM.flag(record) & SAM.FLAG_REVERSE ≠ 0         # and if it's reverse complement
            return true
        end
    elseif BAM.flag(record) & SAM.FLAG_READ2 ≠ 0           # if it's read2
        if BAM.flag(record) & SAM.FLAG_REVERSE == 0        # and it's forward
            return true
        end
    elseif BAM.flag(record) & SAM.FLAG_REVERSE ≠ 0        # if it's single-end & reverse complement
            return true
    end
    return false
end

function pairpassesmuster(flags1,flags2)
    flags1 & SAM.FLAG_PROPER_PAIR == 0 && return false
    flags1 & SAM.FLAG_SECONDARY ≠ 0 && return false
    flags1 & SAM.FLAG_QCFAIL ≠ 0 && return false
    flags1 & SAM.FLAG_DUP ≠ 0 && return false
    flags1 & SAM.FLAG_SUPPLEMENTARY ≠ 0 && return false

    flags2 & SAM.FLAG_PROPER_PAIR == 0 && return false
    flags2 & SAM.FLAG_SECONDARY ≠ 0 && return false
    flags2 & SAM.FLAG_QCFAIL ≠ 0 && return false
    flags2 & SAM.FLAG_DUP ≠ 0 && return false
    flags2 & SAM.FLAG_SUPPLEMENTARY ≠ 0 && return false
    return true
end

function pairedcountbases(refseqs::Array{LongDNASeq},reader::BAM.Reader,mapQ_threshold::Int,baseQ_threshold::Int,contextwindow::Int,utoc::Bool)

    fwd_base_counts = Array{Array{Int}}(undef,length(refseqs))
    rev_base_counts = Array{Array{Int}}(undef,length(refseqs))

    for (index,ref) in enumerate(refseqs)
        fwd_base_counts[index] = zeros(Int64, length(ref), 4)
        rev_base_counts[index] = zeros(Int64, length(ref), 4)
    end

    fwd_bases = Dict(DNA_A=>1, DNA_C=>2, DNA_G=>3, DNA_T=>4)
    rev_bases = Dict(DNA_A=>4, DNA_C=>3, DNA_G=>2, DNA_T=>1)

    queue = Vector{BAM.Record}(undef,0)
    num_pairs = 0

    for next_record in reader
        !BAM.isfilled(next_record) && continue
        #println("reading ",BAM.tempname(next_record)," refID: ",BAM.flag(next_record))
        push!(queue,next_record)
        length(queue) != 2 && continue
        if BAM.tempname(queue[1]) != BAM.tempname(queue[2])
            popfirst!(queue)
            continue
        end

        record1 = queue[1]
        record2 = queue[2]
        #println("read1 ",BAM.tempname(record1)," ",BAM.flag(record1) & SAM.FLAG_READ1," ",BAM.flag(record1) & SAM.FLAG_READ2)
        #println("read2 ",BAM.tempname(record2)," ",BAM.flag(record2) & SAM.FLAG_READ1," ",BAM.flag(record2) & SAM.FLAG_READ2)

        # ignore bad read pairs (not properly mapped, duplicates etc)
        if !pairpassesmuster(BAM.flag(record1),BAM.flag(record2)) || BAM.mappingquality(record1) < mapQ_threshold ||BAM.mappingquality(record2) < mapQ_threshold
            empty!(queue)
            continue
        end

        #println("counting pair ",BAM.tempname(record1)," ",BAM.tempname(record2))

        refindex = BAM.refid(record1) #both reads should have same reference as they are a proper pair

        # progress tracking
        #num_pairs += 1
        #if num_pairs%1000 == 0; print(num_pairs,"\r");end

        refpos1 = BAM.position(record1)
        ops1, oplens1 = BAM.cigar_rle(record1)
        seq1 = AlignedSequence(BAM.sequence(record1), BAM.alignment(record1)).seq
        q1 = BAM.quality(record1)
        range1 = refpos1:first(seq2ref(BAM.alignment(record1),length(seq1)))
        #println(seq1," ",range1)

        refpos2 = BAM.position(record2)
        ops2, oplens2 = BAM.cigar_rle(record2)
        seq2 = AlignedSequence(BAM.sequence(record2), BAM.alignment(record2)).seq
        q2 = BAM.quality(record2)
        range2 = refpos2:first(seq2ref(BAM.alignment(record2),length(seq2)))
        #println(seq2," ",range2)

        forward = isforwardstrand(record1)
        matches = 0
        seq1pos = 1
        for (index, op) in enumerate(ops1), n in 1:oplens1[index]
            t2nt = DNA_N
            if refpos1 in range2
                t2nt = seq2[first(ref2seq(BAM.alignment(record2),refpos1))]
            end
            if isinsertop(op)
                refnt = DNA_Gap
                t1nt = seq1[seq1pos]
                seq1pos += 1
            elseif isdeleteop(op)
                refnt = refseqs[refindex][refpos1]
                refpos1 += 1
                t1nt = DNA_Gap
            else
                refnt = refseqs[refindex][refpos1]
                refpos1 += 1
                t1nt = seq1[seq1pos]
                seq1pos += 1
            end
            if t2nt ≠ DNA_N && t2nt ≠ t1nt; t1nt = DNA_N; end
            if ismatchop(op) && (refnt == t1nt || (refnt == DNA_C && forward && t1nt == DNA_T)
                                || (refnt == DNA_G && !forward && t1nt == DNA_A)
                                || (refnt == DNA_T && forward && utoc && t1nt == DNA_C)
                                || (refnt == DNA_A && !forward && utoc && t1nt == DNA_G))
                matches += 1
            else
                matches = 0
            end
            #println(op," ",refnt," ",t1nt," ",t2nt," ",matches)
            if matches >= contextwindow*2+1
                targetindex = seq1pos-contextwindow-1
                q1[targetindex] < baseQ_threshold && continue
                targetbase = seq1[targetindex]
                if forward
                    base_index = get(fwd_bases,targetbase,0)
                    if base_index > 0; fwd_base_counts[refindex][refpos1-contextwindow-1,base_index] += 1; end
                    #println(BAM.tempname(record1),"\t",refpos1-contextwindow-1,"\tF\t",targetbase)
                else
                    base_index = get(rev_bases,targetbase,0)
                    if base_index > 0; rev_base_counts[refindex][refpos1-contextwindow-1,base_index] += 1; end
                    #println(BAM.tempname(record1),"\t",refpos1-contextwindow-1,"\tR\t",targetbase)
                end
            end
        end
        forward = isforwardstrand(record2)
        matches = 0
        seq2pos = 1
        alreadycounted = false
        for (index, op) in enumerate(ops2), n in 1:oplens2[index]
            t1nt = DNA_N
            if refpos2 in range1
                t1nt = seq1[first(ref2seq(BAM.alignment(record1),refpos2))]
                alreadycounted = true
            end
            if isinsertop(op)
                refnt = DNA_Gap
                t2nt = seq2[seq2pos]
                seq2pos += 1
            elseif isdeleteop(op)
                refnt = refseqs[refindex][refpos2]
                refpos2 += 1
                t2nt = DNA_Gap
            else
                refnt = refseqs[refindex][refpos2]
                refpos2 += 1
                t2nt = seq2[seq2pos]
                seq2pos += 1
            end
            if t1nt ≠ DNA_N && t2nt ≠ t1nt; t2nt = DNA_N; end
            if ismatchop(op) && (refnt == t1nt || (refnt == DNA_C && forward && t1nt == DNA_T)
                                || (refnt == DNA_G && !forward && t1nt == DNA_A)
                                || (refnt == DNA_T && forward && utoc && t1nt == DNA_C)
                                || (refnt == DNA_A && !forward && utoc && t1nt == DNA_G))
                matches += 1
            else
                matches = 0
            end
            if matches >= contextwindow*2+1 && !alreadycounted
                targetindex = seq2pos-contextwindow-1
                q2[targetindex] < baseQ_threshold && continue
                targetbase = seq2[targetindex]
                if forward
                    base_index = get(fwd_bases,targetbase,0)
                    if base_index > 0; fwd_base_counts[refindex][refpos2-contextwindow-1,base_index] += 1; end
                    #println(BAM.tempname(record2),"\t",refpos1-contextwindow-1,"\tF\t",targetbase)
                else
                    base_index = get(rev_bases,targetbase,0)
                    if base_index > 0; rev_base_counts[refindex][refpos2-contextwindow-1,base_index] += 1; end
                    #println(BAM.tempname(record2),"\t",refpos1-contextwindow-1,"\tR\t",targetbase)
                end
            end
        end
        empty!(queue)
    end
    return fwd_base_counts,rev_base_counts
end

function unpairpassesmuster(flags)
    flags & SAM.FLAG_SECONDARY ≠ 0 && return false
    flags & SAM.FLAG_QCFAIL ≠ 0 && return false
    flags & SAM.FLAG_DUP ≠ 0 && return false
    flags & SAM.FLAG_SUPPLEMENTARY ≠ 0 && return false
    return true
end

function unpairedcountbases(refseqs::Array{LongDNASeq},reader::BAM.Reader,mapQ_threshold::Int,baseQ_threshold::Int,contextwindow::Int,utoc::Bool)
    fwd_base_counts = Array{Array{Int}}(undef,length(refseqs))
    rev_base_counts = Array{Array{Int}}(undef,length(refseqs))

    for (index,ref) in enumerate(refseqs)
        fwd_base_counts[index] = zeros(Int64, length(ref), 4)
        rev_base_counts[index] = zeros(Int64, length(ref), 4)
    end

    fwd_bases = Dict(DNA_A=>1, DNA_C=>2, DNA_G=>3, DNA_T=>4)
    rev_bases = Dict(DNA_A=>4, DNA_C=>3, DNA_G=>2, DNA_T=>1)

    record1 = BAM.Record()

    for record1 in reader
        !BAM.isfilled(record1) && continue

        # ignore bad read pairs (not properly mapped, duplicates etc)
        if !unpairpassesmuster(BAM.flag(record1)) || BAM.mappingquality(record1) < mapQ_threshold
            continue
        end

        refindex = BAM.refid(record1)

        refpos1 = BAM.position(record1)
        seq1 = AlignedSequence(BAM.sequence(record1), BAM.alignment(record1)).seq
        ops, oplens = BAM.cigar_rle(record1)
        q1 = BAM.quality(record1)
        range1 = refpos1:first(seq2ref(BAM.alignment(record1),length(seq1)))

        forward = isforwardstrand(record1)
        matches = 0
        seqpos = 1
        for (index, op) in enumerate(ops), n in 1:oplens[index]
            if isinsertop(op)
                refnt = DNA_Gap
                t1nt = seq1[seqpos]
                seqpos += 1
            elseif isdeleteop(op)
                refnt = refseqs[refindex][refpos1]
                refpos1 += 1
                t1nt = DNA_Gap
            else
                refnt = refseqs[refindex][refpos1]
                refpos1 += 1
                t1nt = seq1[seqpos]
                seqpos += 1
            end
            if ismatchop(op) && (refnt == t1nt || (refnt == DNA_C && forward && t1nt == DNA_T)
                                || (refnt == DNA_G && !forward && t1nt == DNA_A)
                                || (refnt == DNA_T && forward && utoc && t1nt == DNA_C)
                                || (refnt == DNA_A && !forward && utoc && t1nt == DNA_G))
                matches += 1
            else
                matches = 0
            end
            if matches >= contextwindow*2+1
                targetindex = seqpos-contextwindow-1
                q1[targetindex] < baseQ_threshold && continue
                targetbase = seq1[targetindex]
                if forward
                    base_index = get(fwd_bases,targetbase,0)
                    if base_index > 0; fwd_base_counts[refindex][refpos1-contextwindow-1,base_index] += 1; end
                    #println(BAM.tempname(record1),"\t",refpos1-contextwindow-1,"\tF\t",targetbase)
                else
                    base_index = get(rev_bases,targetbase,0)
                    if base_index > 0; rev_base_counts[refindex][refpos1-contextwindow-1,base_index] += 1; end
                    #println(BAM.tempname(record1),"\t",refpos1-contextwindow-1,"\tR\t",targetbase)
                end
            end
        end
    end
    return fwd_base_counts,rev_base_counts
end
