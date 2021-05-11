# Notes on Biojulia

[TOC]

## 1. Automa.jl

> The first step of any bioinformatics project is to define a new file format, incompatible with all previous ones.

Automa是构建文本解析器的julia包

TODO: 这个比较底层，后边要构建自己的解释器的时候再学吧！

## 2. 怎么写出更快的程序

TODO: 后边再看！

## 3. BioSequences.jl

### 3.1 BioSymbols

- `alphabet.([DNA, RNA, AminoAcid])`: 返回bioseq定义的BioSymbols常量
- `gap()`
- `iscompatible(x, y)`: x是否属于y（简并）
- `isambiguous()` : 是否是简并碱基/AA

### 3.2 BioSeq Types

#### BioSequence 抽象类型

`BioSequence{A<:Alphabet}`: BioSequence 是一个抽象类型，包括实体类型Alphabet

!!! 任意BioSeq相关的类型必须继承自BioSequence{alphabet}，必须包含如下方法：

- `BioSequences.encoded_data`
- `Base.length`
- `BioSequence.encoded_data_type`: seq type
- `BioSequences.encoded_data_eltype`: element type
- `BioSequences.Alphabet`: Alphabet type
- `BioSymbols.alphabet`
- `BioSequences.BitsPerSymbol`: ->Type: bits required to store the seq
- `BioSequences.bits_per_symbol` -> funcion

#### LongSequence 类型

> LongSequence{A<:Alphabet} <: BioSequence{A}

TIPS:
默认的LongDNASeq和LongRNASeq用的是`DNAAlphabet{4}`类型(或RNA)，这个类型允许有简并碱基，如果你确定你的数据中没有简并碱基，可以把类型改成`DNAAlphabet{2}`，可以用bitwise操作，更快！

#### Kmers & Skipmers

- Kmers：

    ```julia
    # Kmers are immutable
    # seq <: UInt64 && K <= 32
    Mer{A<:NucleicAcidAlphabet{2},K}

    # BigMer:
    # seq <: UInt128 && K <= 64 (最大用63，奇数mer可以避免回文)
    BigMer{A<:NucleicAcidAlphabet{2},K}

    # AbstractMer: Mer 和 BigMer 的抽象类型
    ```

    提供了方便的别名：

    Type alias | Type
    -----------|-----
    DNAMer{K} | Mer{DNAAlphabet{2},K}
    RNAMer{K} | Mer{RNAAlphabet{2},K}
    DNAKmer | DNAMer{31}
    RNAKmer | RNAMer{31}
    BigDNAMer{K} | BigMer{DNAAlphabet{2},K}
    BigRNAMer{K} | BigMer{RNAAlphabet{2},K}
    BigDNAKmer | BigMer{DNAAlphabet{2},63}
    BigRNAKmer | BigMer{RNAAlphabet{2},63}
    DNACodon | DNAMer{3}
    RNACodon | RNAMer{3}

- Skipmers

    - 连续的kmer容易受变异的影响，也更容易受重复的影响。
    
    - `Skipmers`: Skipmers are a generalisation of the concept of a kmer. They are created using a cyclic pattern of used-and-skipped positions which achieves increased entropy and tolerance to nucleotide substitution differences by following some simple rules.
    ![](https://juliahub.com/docs/BioSequences/i7zyu/2.0.5/skipmers.png)

    - Skipmer可以用迭代器从kmer中产生

#### Ref seq

`ReferenceSequence`类型只接受“ATCGN”，并对N进行了压缩，所以堪比只有ATCG的2bits类型，但是ReferenceSequence类型是不可变的。


```julia
seq = ReferenceSequence(dna"NNCGTATTTCN")

seq[1]
# DNA_N
seq[2:6]
# NCGTA
```

### 3.3 创建序列

```julia
LongDNASeq("TTANC")
LongSequence{DNAAlphabet{2}}("TTAGC")
LongRNASeq("UUANC")
DNAMer{8}("ATCGATCG")

# from vectors or arrays of BioSymbol:
LongDNASeq([DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])
DNAMer{3}([DNA_T, DNA_A, DNA_G])

# from other sequences
LongDNASeq(LongDNASeq("ATCG"), LongDNASeq("NNNN"), LongDNASeq("TCGA"))
LongDNASeq("ATCG") * LongDNASeq("NNNN")
repeat(LongDNASeq("TA"), 10)
LongDNASeq("TA") ^ 10

# mker <=> long seq
m = DNAMer{5}("ATCGA")
LongSequence(m)
DNAMer{5}(LongSequence(m))
```

### 3.4 改变序列类型


```julia
julia> dna = dna"ttangtagaccg"
12nt DNA Sequence:
TTANGTAGACCG

julia> dnastr = convert(String, dna)
"TTANGTAGACCG"

# 把dna放到String向量中，强制类型转换
julia> String[dna]
1-element Vector{String}:
 "TTANGTAGACCG"
```

### 3.5 序列字符串宏

```julia

dna"atcgn"
rna"aucgn"
aa"arnd"

# 可以加尾标 d 和 s 表示序列是动态变化的 还是静态的（默认是静态的）：

julia> function foo()
       s = dna"CCT"
       push!(s, DNA_A)
       end
foo (generic function with 1 method)

julia> foo()
4nt DNA Sequence:
CCTA

julia> foo()
5nt DNA Sequence:
CCTAA

julia> foo()
6nt DNA Sequence:
CCTAAA

# 不加d尾标，默认s是静态的，所以每次都累加了A
# 因为用的是宏定义的序列，宏在函数编译之前就先完成了分配（且只分配一次），所以每次调用s（虽然是局部变量），指向的都是同一个已经分配的序列

# 加上尾标d，就会在函数运行的同时分配序列信息，就不会累加了：
julia> function foo()
       s = dna"CTT"d
       push!(s, DNA_A)
       end
foo (generic function with 1 method)

julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

# QUST: 那我不用宏定义呢？
julia> function foo()
           s = LongDNASeq("CTT")
           push!(s, DNA_A)
           end
foo (generic function with 1 method)

julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

# 哈哈！ 果然就没事了！ 宏编程直呼内行！
```

kmer 字符串宏：

```julia
mer"ATCG"
mer"ATCG"dna
mer"AUCG"rna
bigmer"ATCG"dna
```

### 3.6 序列索引和修改

#### Indexing

```julia
seq = dna"ACGTACN"
seq[5]
seq[5] = DNA_N
seq
# ACGTNCN

# !!! kmer 类型不能用范围索引：
mer"ATCG"[2] # fine
mer"ATCG"[2:4] # oops!
# DEBUG: WHY ???

# !!! long seq 类型中，创建子序列并不会增加拷贝，因为实际上是创建了指向老序列的索引，但是编辑子序列并不会改变老序列， 神奇：
seq = dna"AAAA"
subseq = seq[1:2]
subseq[2] = DNA_T
julia> subseq
# 2nt DNA Sequence:
# AT
julia> seq
# 4nt DNA Sequence:
# AAAA
# NOTE: 其实是因为 long seq 序列编辑后 会先检查这个序列的指针是不是也指向了其他序列，只有当这个指针指向了其他序列的时候，才会新建一个新拷贝。流弊！

# 序列编辑：
push!
pop!
insert!
append!
union!
deleteat!
resize!
empty!

reverse!
complement!
reverse_complement!
ungap!
canonical!
#   A canonical sequence is the numerical lesser of a k-mer and its reverse complement.
#   This is useful in hashing/counting sequences in data that is not strand specific,
#   and thus observing the short sequence is equivalent to observing its reverse
#   complement.

# QUST: 还是不太懂，啥是canonical啊！！
# 懂了！ 就是序列 A 和它的反向互补序列 A' 中， 那个排序更靠前的序列：

julia> canonical(dna"AAATT")
5nt DNA Sequence:
AAATT

julia> canonical(dna"TAATT")
5nt DNA Sequence:
AATTA

translate
# see `ncbi_trans_table`

```

### 3.7 predicates （判断，断言）

```julia
isrepetitive(seq, n::Int) # ture if seq包含长度至少为n的重复
ispalindromic(seq) # true if seq是回文
hasambiguity(seq)
iscanonical(seq)
```

### 3.8 random seq


```julia

randseq([rng::AbstractRNG], A::Alphabet, sp::Sampler, len::Integer)

# rng：随机生成方式 (see Random.jl), 可选项
# A：序列集和
# sp: 随机规则，比如频率等等
# len: 长度

randdnaseq()
randrnaseq()
randaaseq()

# 生成一个随机样本库集合
SamplerUniform{T}
sp = SamplerUniform(rna"ACGU")

# 带权重的随机样本库集合
SampleWeighted{T}
# 依次分配权重：
# A -> 0.24
# C -> 0.24
# G -> 0.24
# U -> 0.24
# N -> 1 - (0.24 * 4) = 0.04
sp = SamplerWeighted(rna"ACGUN", fill(0.24, 4))

# 随机kmer: 直接用基础的随机函数
rand(DNAMer{7})
rand(BigDNAMer{63})

```

### 3.9 序列查找

- 精确查找

    ```julia

    seq = dna"ACGTAACTGGGGG"
    query = dna"ACGT"
    findfirst(query, seq)
    findlast(query, seq)
    occursin(query, seq)

    # 精确查找的时候支持简并碱基，但查单个碱基symbol的时候除外：
    findfirst(dna"CNT", dna"ACNT")
    # 2:4

    # DNA_N symbol 不会被展开, 而且匹配碱基symbol（实际上是char），返回的是单个index，而不是range：
    findfirst(DNA_N, dna"ACNT")
    # 3
    findfirst(dna"N", dna"ACNT")
    # 1:1

    # query在每次查找之前，都会被转换成 ExactSearchQuery 类型，所以如果同一个query要用很多次的话，提前赋值成ExactSearchQuery类型就可以快很多：
    query = ExactSearchQuery(dna"ACG")
    findfirst(query, seq)

    ```

- 近似查找

    允许特定数目的错配，参见 [Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance)

    ```julia
    approxsearch(seq, dna"ACGG", 0)
    approxsearch(seq, dna"ACGG", 1) # 最多一个错配
    approxsearch(seq, dna"ACGG", 2) # 最多两个错配

    # approxsearch 返回的是第一次匹配的位置a 和 a - 允许错配数（m） 生成的 range：（为啥这么整啊？）
    julia> for i in 0:3
        print("mismatch $i : ")
        println(approxsearch(dna"TTTATCAAA", dna"CTTA", i))
        end
    # mismatch 0 : 0:-1
    # mismatch 1 : 1:4
    # mismatch 2 : 1:2
    # mismatch 3 : 1:1
    ```

- 正则匹配

    ```julia
    # 用 biore 保留字符串：
    # 尾部加dna rna aa表示匹配类型，也可简写成d r a
    match(biore"A+C*"dna, dna"AAACCC")
    # 如果匹配成功，返回RegexMatch数据类型，否则返回nothing
    occursin(biore"A+C*"d, dna"ACCG")

    matched() # 获取匹配
    eachmatch() # 匹配迭代器
    collect() # 收集迭代器或collection中的元素成一个arry

    collect(matched(x) for x in eachmatch(biore"TATA*?"d, dna"TATTATAATTA"))
    ```

    Syntax | Description | Example
    -------|-------------|--------
    `|` | alternation | "A|T" matches "A" and "T"
    `*` | zero or more times repeat | "TA*" matches "T", "TA" and "TAA"
    `+` | one or more times repeat | "TA+" matches "TA" and "TAA"
    `?` | zero or one time | "TA?" matches "T" and "TA"
    `{n,}` | n or more times repeat | "A{3,}" matches "AAA" and "AAAA"
    `{n,m}` | n-m times repeat | "A{3,5}" matches "AAA", "AAAA" and "AAAAA"
    `^` | the start of the sequence | "^TAN*" matches "TATGT"
    `$` | the end of the sequence | "N*TA$" matches "GCTA"
    `(...)` | pattern grouping | "(TA)+" matches "TA" and "TATA"
    `[...]` | one of symbols | "[ACG]+" matches "AGGC"

    NOTE:
    1. 支持通配符： `biore"N"d == biore"[ATCG]"d`
    2. 忽略空格：`biore"A C    G"d == biore"ACG"d`

- 位置权重矩阵：PWM and PFM

    see [Position_weight_matrix](https://en.wikipedia.org/wiki/Position_weight_matrix)


    ```julia
    # seqs -> PFW -> PWM
    julia> kmers = DNAMer.(["TTA", "CTA", "ACA", "TCA", "GTA"])
    5-element Vector{Mer{DNAAlphabet{2}, 3}}:
    TTA
    CTA
    ACA
    TCA
    GTA

    julia> pfm = PFM(kmers)
    4×3 PFM{DNA, Int64}:
    A  1  0  5
    C  1  2  0
    G  1  0  0
    T  2  3  0

    julia> pwm = PWM(pfm)
    4×3 PWM{DNA, Float64}:
    A -0.321928 -Inf       2.0
    C -0.321928  0.678072 -Inf
    G -0.321928 -Inf      -Inf
    T  0.678072  1.26303  -Inf

    julia> pwm = PWM(pfm .+ 0.01, prior=[0.2, 0.3, 0.3, 0.2])  # 加0.01，防止Inf，设置权重匹配GC-rich prior
    4×3 PWM{DNA, Float64}:
    A  0.00285965 -6.65535   2.31331
    C -0.582103    0.410737 -7.24031
    G -0.582103   -7.24031  -7.24031
    T  0.9957      1.57827  -6.65535

    #~q 加0.01影响大不大？是不是要加更小的数值？
    ```

### 3.10 序列解标签 demuliplexing

```julia
julia> barcodes = LongDNASeq.(["ATGG", "CAGA", "GGAA", "TACG"]);

julia> dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:hamming)
Demultiplexer{LongSequence{DNAAlphabet{4}}}:
  distance: hamming
  number of barcodes: 4
  number of correctable errors: 1

foreach(x -> println(demultiplex(dplxr, x)), [dna"ATGGCGNT", dna"CAGGCGNT", dna"GGAACGNT", dna"TGACCGNT"])
(1, 0)
(2, 1)
(3, 0)
(0, -1)

julia> demultiplex(dplxr, dna"TGACCGNT", false)  # linear search off (default)
(0, -1)
#~n 设置第三个参数为true，会在index搜索完以后，如果找不到合适的barcode，再用lnear search搜索一个最近的结果返回，即便错配很多
julia> demultiplex(dplxr, dna"TGACCGNT", true)   # linear search on
(3, 2)
```

### 3.11 序列组成

```julia
# Composition type
# composition function
julia> a = dna"AAAAAAAATTTTTT"; b = a; c = a[5:10]
6nt DNA Sequence:
AAAATT

julia> composition([a, b, c])
Composition{LongDNASeq} with 2 entries:
  AAAAAAAATTTTTT => 2
  AAAATT         => 1

julia> comp = composition(dna"ACGAG");

julia> comp[DNA_A]
2
julia> comp[DNA_T]
0

```

Composition 看上去像具名数组或者字典，但是有一些区别：

1. composition最开始会初始化所有要用到的元素，即便最后没有结果（记为0）
2. merge! 方法会合并两个composition的结果，但在字典中会替换value

```julia
# 用merge!累计统计多条序列：

# initiaize empty composion
comp = composition(dna"")
for seq in seqs
    merge!(comp, composition(seq))
end

# 如果写成函数可以这样：
foldl(x, y) -> merge(x, composition(y)), composition(dna""), seqs)

# compositon kmer:
composition(each(DNAMer{4}, dna"ACGT" ^ 100))
# Composition{DNAMer{4}} with 4 entries:
#   ACGT => 100
#   CGTA => 99
#   GTAC => 99
#   TACG => 99
```

### 3.12 序列迭代器

```julia

# 序列是碱基的迭代器：
for nt in dna"ATCG"
    println(nt)
end

# kmer迭代： each 方法
each(::Type{T}, seq::BioSequence, step::Int) where {T<:AbstractMer}
each(DNAMer{27}, seq)

each(::Type{T}, seq::BioSequence, step::Integer) where T<:AbstractMer in BioSequences
each(DNAMer{27}, seq, 10) # step 10

# skipmer生成
each(::Type{T}, seq::BioSequence, bases_per_cycle::Int64, cycle_len::Int64) where T<:AbstractMer in BioSequences
each(DNAMer{27}, seq, 2, 3)

# 根据碱基位置和碱基进行过滤
each(f::Function, seq::BioSequence) = ConditionIterator(f, seq)
for (pos, nuc) in each(isambiguous, dna"NATTCGRATY")
    println("Position $pos is $nuc")
end

```

### 3.13 序列统计

```julia
# 拓展了count方法：
count(isambiguous, dna"ATCGM")
# match / mismatch
count(==, dna"ATCGM", dna"GCCGM")
count(!=, dna"ATCGM", dna"GCCGM")
#~q 长度不一样会咋样？
julia> count(==, dna"ATCGAT", dna"ATCGAT")
6
julia> count(==, dna"ATCGAT", dna"NATCGAT")
0
# 只会按位两两比较，不比对

# 大部分情况的count，是把function map给每个位置的碱基，但在某些情况下（BioSymbols.jl和BioSequences.jl中的函数），会用到“位并行”技术加速: BioSequences.count_*_bitpar; 如果想强制用基础的map方法，直接调用 BioSequences.count_naive

```

Alias function | Base.count call(s)
---------------|-------------------
n_ambiguous | count(isambiguous, seq), count(isambiguous, seqa, seqb)
n_certain | count(iscertain, seq), count(iscertain, seqa, seqb)
n_gap | count(isgap, seq), count(isgap, seqa, seqb)
matches | count(==, seqa, seqb)
mismatches | count(!=, seqa, seqb)


### 3.14 I / O

- FASTA

    ```julia
    using FASTX
    # as IOStream
    r = FASTA.Reader(open("seqin.fa", "r"))
    w = FASTA.Reader(open("seqout.fa", "w"))
    # as a method of open
    r = open(FASTA.Reader, "my-seqs.fasta")
    w = open(FASTA.Writer, "my-out.fasta")
    # read gz:
    reader = FASTA.Reader(GzipDecompressorStream(open("my-reads.fasta.gz")))
    # combine with loop for:
    for record in reader
        ## do something
    end
    close(reader)
    # while:
    reader = open(FASTA.Reader, "my-seqs.fasta")
    record = FASTA.Record()
    while !eof(reader)
        read!(reader, record)
        ## Do something.
    end

    # read fai index, and jump to specific sequence:
    reader = open(FASTA.Reader, "sacCer.fa", index = "sacCer.fa.fai")
    chrIV = reader["chrIV"]  # directly read sequences called chrIV.
    close(reader)

    # 读取的结果是FASTA.Record类型，有以下属性：
    FASTA.hasidentifier
    FASTA.identifier
    FASTA.hasdescription
    FASTA.description
    FASTA.hassequence
    FASTA.sequence

    # write a BioSequence to Fasta:
    # get a FASTA.Record type first:
    using BioSequences
    x = dna"aaaaatttttcccccggggg"
    rec = FASTA.Record("MySeq", x) # assign seq with name "MySeq", and cvrt to FASTA.Record
    w = open(FASTA.Writer, "my-out.fasta")
    write(w, rec)
    close(w)

    # combine filehandle with do-block ti avoid close manually:
    open(FASTA.Reader, "my-reads.fasta") do reader
        for record in reader
            ## Do something
        end
    end
    ```

- FASTQ

    ```julia
    FASTQ.Reader
    FASTQ.Writer
    FASTQ.Record # type
    # getters and setters：
    FASTQ.hasidentifier
    FASTQ.identifier
    FASTQ.hasdescription
    FASTQ.description
    FASTQ.hassequence
    FASTQ.sequence
    FASTQ.hasquality # 比fasta多的
    FASTQ.quality # 比fasta多的

    # Quality encoding：
    FASTQ.SANGER_QUAL_ENCODING
    FASTQ.SOLEXA_QUAL_ENCODING
    FASTQ.ILLUMINA13_QUAL_ENCODING
    FASTQ.ILLUMINA15_QUAL_ENCODING
    FASTQ.ILLUMINA18_QUAL_ENCODING
    # https://github.com/BioJulia/FASTX.jl/blob/9a69972da599ffac513e8d10afd4c7dfabc2dcb0/src/fastq/quality.jl#L29
    ```

- 2bit

    ```julia
    using TwoBit

    TwoBit.Reader(input::IO)
    TwoBit.seqnames(reader) # 提取序列名字
    TwoBitWriter(output::IO, names::AbstractVector)
    TowBit.Record #Type

    TwoBit.sequence
    TwoBit.hassequence
    TwoBit.maskedblocks #Get the masked blocks.

    ```

## 4. BioStructures.jl

主要是处理蛋白结构的，先略过吧

## 5.  GenomeGraphs

> THIS IS WHY I AM HERE !

### 5.1 Construct a DBG

```julia
# 1. prepare the reads
using GenomeGraphs, ReadDatastores
ws = WorkSpace() # create a workspace

# read test reads
fwq = open(FASTQ.Reader, "test/ecoli_tester_R1.fastq")
rvq = open(FASTQ.Reader, "test/ecoli_tester_R2.fastq")
# make them PE:
#~t PairedReads{A}(rdrx::FASTQ.Reader, rdry::FASTQ.Reader, outfile::String, name::Union{String,Symbol}, minsize::Integer, maxsize::Integer, fragsize::Integer, orientation::PairedReadOrientation) where {A<:DNAAlphabet}
# 1. new version need to declare the Alphabet type
# 2. FwRv -> common pe seq; some long-mate seq would be RvFw
ds = PairedReads{DNAAlphabet{2}}(fwq, rvq, "ecoli-test-paired", "my-ecoli-test", 250, 300, 0, FwRv)

# add data to workspace
add_paired_reads!(ws, ds)

# run the dbg process
#~q f**k! 找了半天才找到怎么用， 为啥不把readme写清楚！！

#~t 以下是dbg!的源码：
# """
#     dbg!(ws::WorkSpace, ds::String, ::Type{M}, min_freq::Integer, name::Symbol) where {M<:AbstractMer}
# """
# function dbg!(ws::WorkSpace, ::Type{M}, min_freq::Integer, name::Symbol) where {M<:AbstractMer}
#     reads = paired_reads(ws, name)
#     _dbg!(ws.sdg, reads, M, UInt8(min_freq))
#     return ws
# end
#~w 函数的说明写错了！应该是
# dbg!(ws::WorkSpace, ::Type{M}, min_freq::Integer, name::Symbol) where {M<:AbstractMer}
# Symbol是pe reads的名字的Symbol，使用的时候注意“-”字符不能直接解析，所以不能用:my-ecoli-test这样声明，要用Symbol()方法
# 真难顶啊！！
dbg!(ws, BigDNAMer{61}, 10, Symbol("my-ecoli-test"))
# [ Info: Counting kmers in datastore
# [ Info: CConstructing compressed de-bruijn graph from 0 61-mers
# [ Info: Constructing unitigs from 0 61-mers
# [ Info: Constructed 0 unitigs
# [ Info: Done cConstructing compressed de-bruijn graph from 0 61-mers
# Graph Genome workspace
#  Sequence distance graph (0 nodes)
#  Paired Read Datastores (1):
#   'my-ecoli-test': 20 reads (10 pairs)
#  Long Read Datastores (0):
#  Linked Read Datastores (0):
#  Mer Count Datastores (0):
#~q 怎么把dbg输出成gfa？
#~n GenomeGraph 的官方manual写的太草了， 直接看github中的doc学吧：
#~n https://github.com/BioJulia/GenomeGraphs.jl/blob/master/docs/src/DeBruijnGraph.md

# >>>>>>>>>>>>>>>>>>>>>>>> new notes of DBG >>>>>>>>>>>>>>>>>>>>>>>>

# <<<<<<<<<<<<<<<<<<<<<<<< new notes of DBG <<<<<<<<<<<<<<<<<<<<<<<<
```

## 6. BioAlignments

> https://biojulia.net/BioAlignments.jl/stable/

Alignment representation

```julia
# AlignmentAnchor类型
struct AlignmentAnchor
    seqpos::Int
    refpos::Int
    op::Operation
end
```

Alignment 的表示形式：

![](https://biojulia.net/BioAlignments.jl/stable/assets/alignment.svg)

这个比对结果：

                  0   4        9  12 15     19
                  |   |        |  |  |      |
        query:     TGGC----ATCATTTAACG---CAAG
    reference: AGGGTGGCATTTATCAG---ACGTTTCGAGAC
                  |   |   |    |     |  |   |
                  4   8   12   17    20 23  27

可以用Anchor类型表示成：

```julia
# query idx, ref idx, OP_TAGS
#~q 1. del会和前一个的match共用query idx；ins会和前一个的match共用ref idx
 Alignment(
       [
           AlignmentAnchor( 0,  4, OP_START),
           AlignmentAnchor( 4,  8, OP_MATCH),
           AlignmentAnchor( 4, 12, OP_DELETE),
           AlignmentAnchor( 9, 17, OP_MATCH),
           AlignmentAnchor(12, 17, OP_INSERT),
           AlignmentAnchor(15, 20, OP_MATCH),
           AlignmentAnchor(15, 23, OP_DELETE),
           AlignmentAnchor(19, 27, OP_MATCH),
       ]
    )
```

Aln的类型，跟SAM类似：

| Operation            | Operation Type     | Description                                                                     |
| :------------------- | :----------------- | :------------------------------------------------------------------------------ |
| `OP_MATCH`           | match              | non-specific match                                                              |
| `OP_INSERT`          | insert             | insertion into reference sequence                                               |
| `OP_DELETE`          | delete             | deletion from reference sequence                                                |
| `OP_SKIP`            | delete             | (typically long) deletion from the reference, e.g. due to RNA splicing          |
| `OP_SOFT_CLIP`       | insert             | sequence removed from the beginning or end of the query sequence but stored     |
| `OP_HARD_CLIP`       | insert             | sequence removed from the beginning or end of the query sequence and not stored |
| `OP_PAD`             | special            | not currently supported, but present for SAM/BAM compatibility                  |
| `OP_SEQ_MATCH`       | match              | match operation with matching sequence positions                                |
| `OP_SEQ_MISMATCH`    | match              | match operation with mismatching sequence positions                             |
| `OP_BACK`            | special            | not currently supported, but present for SAM/BAM compatibility                  |
| `OP_START`           | special            | indicate the start of an alignment within the reference and query sequence      |

```julia
# aln的操作都存到了Operation类型：
convert(Operation, 'M')
# OP_MATCH

# AlignedSequence 类型：
julia> AlignedSequence(  # pass an Alignment object
           dna"ACGTAT",
           Alignment([
               AlignmentAnchor(0, 0, OP_START),
               AlignmentAnchor(3, 3, OP_MATCH),
               AlignmentAnchor(6, 3, OP_INSERT)
           ])
       )
# ···---
# ACGTAT

julia> seq = dna"ACGT--AAT--"
11nt DNA Sequence:
ACGT--AAT--

julia> ref = dna"ACGTTTAT-GG"
11nt DNA Sequence:
ACGTTTAT-GG

julia> AlignedSequence(seq, ref)
········-··
ACGT--AAT--

# AlignedSequence 有5个方法，略
```

比对

```julia
# pairwise比对
julia> s1 = dna"CCTAGGAGGG";

julia> s2 = dna"ACCTGGTATGATAGCG";

julia> scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);

julia> res = pairalign(GlobalAlignment(), s1, s2, scoremodel)  # run pairwise alignment
PairwiseAlignmentResult{Int64,BioSequences.LongSequence{BioSequences.DNAAlphabet{4}},BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}:
  score: 13
  seq:  0 -CCTAGG------AGGG 10
           ||| ||      || |
  ref:  1 ACCT-GGTATGATAGCG 16


julia> score(res)  # get the achieved score of this alignment
13

julia> aln = alignment(res)
PairwiseAlignment{BioSequences.LongSequence{BioSequences.DNAAlphabet{4}},BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}:
  seq:  0 -CCTAGG------AGGG 10
           ||| ||      || |
  ref:  1 ACCT-GGTATGATAGCG 16

# 以上两个函数等价于直接 res.value  res.aln

julia> res.
aln     isscore  value
julia> res.isscore
true

julia> count_matches(aln)
8

julia> count_mismatches(aln)
1

julia> count_insertions(aln)
1

julia> count_deletions(aln)
7

julia> count_aligned(aln)
17

julia> collect(aln)  # pairwise alignment is iterable
17-element Array{Tuple{DNA,DNA},1}:
 (DNA_Gap, DNA_A)
 (DNA_C, DNA_C)
 (DNA_C, DNA_C)
 (DNA_T, DNA_T)
 (DNA_A, DNA_Gap)
 (DNA_G, DNA_G)
 (DNA_G, DNA_G)
 (DNA_Gap, DNA_T)
 (DNA_Gap, DNA_A)
 (DNA_Gap, DNA_T)
 (DNA_Gap, DNA_G)
 (DNA_Gap, DNA_A)
 (DNA_Gap, DNA_T)
 (DNA_A, DNA_A)
 (DNA_G, DNA_G)
 (DNA_G, DNA_C)
 (DNA_G, DNA_G)

julia> LongDNASeq([x for (x, _) in aln])  # create aligned `s1` with gaps
17nt DNA Sequence:
-CCTAGG------AGGG

julia> LongDNASeq([y for (_, y) in aln])  # create aligned `s2` with gaps
17nt DNA Sequence:
ACCT-GGTATGATAGCG

# “_” 是 throw away variable

# substitution matrix
# EDNAFULL -> DNA and RNA
EDNAFULL

# PAM (Point Accepted Mutation) and BLOSUM (BLOcks SUbstitution Matrix) -> pr

```

| Matrix   | Constants                                                  |
| :------- | :----------                                                |
| PAM      | `PAM30`, `PAM70`, `PAM250`                                 |
| BLOSUM   | `BLOSUM45`, `BLOSUM50`, `BLOSUM62`, `BLOSUM80`, `BLOSUM90` |

```jldoctest
# DichotomousSubstitutionMatrix 定义一种最简单的矩阵，只有match和mismathc的区别，当对精度没有要求，但是对速度有要求的时候，可以用这个
julia> submat = DichotomousSubstitutionMatrix(1, -1)
DichotomousSubstitutionMatrix{Int64}:
     match =  1
  mismatch = -1

julia> submat['A','A']  # match
1

julia> submat['A','B']  # mismatch
-1

```

几个有用的API

```julia
Alignment() # 有很多方法，其中一个可以 cigar --> aln
cigar(aln::Alignment)

AffineGapScoreModel(submat, gap_open, gap_extend)
AffineGapScoreModel(submat, gap_open=, gap_extend=)
AffineGapScoreModel(match=, mismatch=, gap_open=, gap_extend=)

CostModel(submat, insertion, deletion)
CostModel(submat, insertion=, deletion=)
CostModel(match=, mismatch=, insertion=, deletion=)

using BioSequences
using BioAlignments

# create affine gap scoring model
affinegap = AffineGapScoreModel(
    match=5,
    mismatch=-4,
    gap_open=-5,
    gap_extend=-3
)

# run global alignment between two DNA sequences
pairalign(GlobalAlignment(), dna"AGGTAG", dna"ATTG", affinegap)

# run local alignment between two DNA sequences
pairalign(LocalAlignment(), dna"AGGTAG", dna"ATTG", affinegap)

# you cannot specify a cost model in LevenshteinDistance
pairalign(LevenshteinDistance(), dna"AGGTAG", dna"ATTG")

#~q 没有PairwiseAlignmentResult 转 Alignment的方法；没有Cigar类型；没有MSA的方法
```

> 比对算法相关 https://github.com/BioJulia/BioAlignments.jl/tree/master/src/pairwise/algorithms

## 7. XAM

> https://github.com/BioJulia/XAM.jl/blob/develop/docs/src/man/hts-files.md

```julia
using XAM
reader = open(BAM.Reader, "data.bam")
record = BAM.Record() # 定义一个BAM.Reccord类型的空变量
# 因为bam文件一般都很大，所以最好随操作随扔，所以把bam记录的变量在循环外定义
# 最好不要用for rcd in reader 这种方式
while !eof(reader)
    empty!(record)
    read!(reader, record)
    # do something
end

# header() 获取 sam bam的头文件
find(header(reader), "SQ")

# sam bam 相关的方法：
#~w 太多了，用的时候再细看吧！
julia> XAM.BAM.
AuxData              cigar_rle             hasrefname            nextrefid
BAI                  compare_intervals     hasrightposition      nextrefname
FIXED_FIELDS_BYTES   data_size             hasseqlength          position
OverlapIterator      eval                  hassequence           quality
OverlapIteratorState extract_cigar_rle     hastemplength         read_bai
Reader               findauxtag            hastempname           refid
Record               flag                  include               reflen
Writer               getauxvalue           init_bam_reader       refname
_read!               hasalignment          ismapped              rightposition
alignlength          hasauxdata            isnextmapped          seqlength
alignment            hasflag               ispositivestrand      seqname_length
auxdata              hasmappingquality     isprimary             sequence
auxdata_position     hasnextposition       loadauxtype           templength
checkauxtag          hasnextrefid          loadauxvalue          tempname
checked_refid        hasnextrefname        mappingquality        write_header
checkfilled          hasposition           n_cigar_op
cigar                hasquality            next_tag_position
cigar_position       hasrefid              nextposition
julia> XAM.SAM.
FLAG_DUP                datarange                ismapped                 sam_initcode_record
FLAG_MREVERSE           eval                     ismissing                sam_loopcode_body
FLAG_MUNMAP             findauxtag               isnextmapped             sam_loopcode_header
FLAG_PAIRED             findkey                  isprimary                sam_loopcode_metainfo
FLAG_PROPER_PAIR        flag                     keyvalues                sam_loopcode_record
FLAG_QCFAIL             hasalignment             mappingquality           sam_machine
FLAG_READ1              hascigar                 nextposition             sam_machine_body
FLAG_READ2              hasflag                  nextrefname              sam_machine_header
FLAG_REVERSE            hasmappingquality        parse_hexarray           sam_machine_metainfo
FLAG_SECONDARY          hasnextposition          parse_typedarray         sam_machine_record
FLAG_SUPPLEMENTARY      hasnextrefname           position                 sam_returncode_body
FLAG_UNMAP              hasposition              quality                  sam_returncode_header
Header                  hasquality               readheader!              sam_returncode_metainfo
MetaInfo                hasrefname               readrecord!              sam_returncode_record
Reader                  hasrightposition         refname                  seqlength
Record                  hasseqlength             rightposition            sequence
Writer                  hassequence              sam_actions_body         tag
action_metainfo         hastemplength            sam_actions_header       templength
alignlength             hastempname              sam_actions_metainfo     tempname
alignment               include                  sam_actions_record       unsafe_parse_decimal
appendfrom!             index!                   sam_context              value
auxdata                 initialize!              sam_initcode_body
checkfilled             iscomment                sam_initcode_header
cigar                   isequalkey               sam_initcode_metainfo

# 获取tag的信息：
for record in open(BAM.Reader, "data.bam")
    nm = record["NM"]::UInt8
    # do something
end

# 支持索引
reader = open(BAM.Reader, "SRR1238088.sort.bam", index="SRR1238088.sort.bam.bai")
for record in eachoverlap(reader, "Chr2", 10000:11000)
    # `record` is a BAM.Record object
    # ...
end
close(reader)
```


## 8. IntervalTrees | GenomicFeatures | GFF3 | GenomicAnnotations

> https://juliahub.com/docs/GenomicFeatures/kSGNI/2.0.4/
> https://juliahub.com/docs/GenomicAnnotations/ckOyU/0.2.3/


```julia
struct Interval{T} <: IntervalTrees.AbstractInterval{Int64}
# The first three fields (seqname, first, and last) are mandatory arguments when constructing the Interval object.

# seqname::String: the sequence name associated with the interval.
# first::Int64: the leftmost position.
# last::Int64: the rightmost position.
# strand::Strand: the strand.
# metadata::T

coverage()
eachoverlap()
......
```

### GenomicAnnotations

- I/O


```julia
GeneBank.Reader # read gbk
GFF.Reader      # read gff3

open(GenBank.Reader, "example.gbk") do reader
    for record in reader
        # do something
    end
end

GeneBank.Writer()
GFF.Writer()

open(GFF.Writer, outfile) do writer
    write(writer, genome)
end

addgene!()
sort!()
delete!()
replace!()
pushproperty!()

sequence()
```



## 9. PopGen

https://biojulia.net/PopGen.jl/docs/

## 10. GeneticVariation

> VCF BCF 文件的处理

## 11. BioServices

> 学习怎么用julia与公共数据库联动

## 12. BioTools

> 目前只有blast，看他们的源码，学习如何用julia调用外部工具
> https://github.com/BioJulia/BioTools.jl/blob/master/src/blast/BLAST.jl
> https://github.com/BioJulia/BioTools.jl/blob/master/src/blast/blastcommandline.jl

```julia
#~w 定义 blastResult 结构，包含几个特征值
struct BLASTResult
    bitscore::Float64
    expect::Float64
    queryname::String
    hitname::String
    hit::BioSequence
    alignment::AlignedSequence
end

#~w 定义一个函数解析xml格式的blast输出结果
# 方法1：输入的是xml文件
function readblastXML(blastrun::AbstractString; seqtype="nucl")
    # blablalba 略过
        # 逐行解析blast结果，转成BLASTResult类型，输出到results array
        push!(results, BLASTResult(bitscore, expect, queryname, hitname, hseq, aln))
    return results
end
# 方法2：输入的是blast运行命令
function readblastXML(blastrun::Cmd; seqtype="nucl")
    # 就是先读blast运行的结果， 再传递给方法一
    return readblastXML(read(blastrun, String), seqtype=seqtype)
end

#~w 定义 blast wrapper 函数们：blastn blastp
# 以blastn为列：
# 方法1： 如果输入的都是字符串，认为是指定的本地文件，就调用本地命令，用readblastXML函数解析成blastResult类型的输出
#~q 最好要加一个判断文件是否存在吧？
function blastn(query::AbstractString, subject::AbstractString, flags=[]; db::Bool=false)
    if db
        results = readblastXML(`blastn -query $query -db $subject $flags -outfmt 5`)
    else
        results = readblastXML(`blastn -query $query -subject $subject $flags -outfmt 5`)
    end
    return results
end
# 方法2：输入都是DNASequence类型的序列，就先存成 fasta 格式的本地文件，然后再调用方法1
#~q 这样如果subject文件很大，或者要复用，不给建索引，很慢啊
function blastn(query::DNASequence, subject::DNASequence, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastn(querypath, subjectpath, flags)
end

# 方法3：同方法2， 如果ref是DNAseq的array
function blastn(query::DNASequence, subject::Vector{DNASequence}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    blastn(querypath, subjectpath, flags)
end
# 方法4：如果ref是本地文件或db，但query是DNAseq，就先把query转成本地文件，然后调用方法1
function blastn(query::DNASequence, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        return blastn(querypath, subject, flags, db=true)
    else
        return blastn(querypath, subject, flags)
    end
end
# 方法5：q t都是seq array
function blastn(query::Vector{DNASequence}, subject::Vector{DNASequence}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastn(querypath, subjectpath, flags)
end
# 方法6：q是seq array， t是path
function blastn(query::Vector{DNASequence}, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        return blastn(querypath, subject, flags, db=true)
    else
        return blastn(querypath, subject, flags)
    end
end
# 方法7：q是path， t是seq array
function blastn(query::AbstractString, subject::Vector{DNASequence}, flags=[])
    subjectpath = makefasta(subject)
    return blastn(query, subjectpath, flags)
end
#~t: 总结：
#~t: blast函数想支持三种输入类型： path，DNASeq，DNASeq array
#~t: 总共有两个变量，所以一共是 2 ^ 3 = 8 种参数排列
#~t: 先指定最基础的方法：path path；其他方法都是这个的变种

# Create temporary multi fasta-formated file for blasting.
#~q 临时目录会定时清理，所以就不管了
#~q 为啥不用FASTX.jl啊，是当时还没有么？
function makefasta(sequences::Vector{T}) where T <: BioSequence
    path, io = mktemp()
    counter = 1
    for sequence in sequences
        write(io, ">$path$counter\n$(convert(String, sequence))\n")
        counter += 1
    end
    close(io)
    return path
end

```

## 13. SubstitutionModels and PhyloTrees

> https://biojulia.net/SubstitutionModels.jl/stable/man/models/
> https://juliahub.com/docs/PhyloTrees/XIHOb/0.11.1/autodocs/

SubstitutionModels

```julia
NucleicAcidSubstitutionModel # type
SubstitutionModels.Q # function，生成NucleicAcidSubstitutionModel 类型的Q矩阵
SubstitutionModels.P # 生成P矩阵

# 预设model：
# JC69
# K80
# F81
# F84
# HKY85
# TN93
# GTR

# If this substitution model allows for unequal base frequencies, methods for _πA, _πC, _πG, and _πT will also need to be defined. With these, SubstitutionModels.jl will calculate the correct Q and P matrices.
```

# OpenMendel

https://github.com/OpenMendel

# 进阶

## Whippet.jl and JWAS.jl 学习如何写一个成熟的包

## Julia调用 minigraph 和 minimap

## MSA pipe

## VNTR identify

## multi-allelic analysis

- LD
- GWAS
- ...

## Database

