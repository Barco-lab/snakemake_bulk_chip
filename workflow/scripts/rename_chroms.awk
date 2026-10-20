# 把第 col 列的染色体名转换成目标命名（UCSC chr1/chrM <-> Ensembl 1/MT），目标里不存在的行丢弃。
#   awk -F'\t' -v OFS='\t' -v col=1 -f rename_chroms.awk 目标染色体列表 输入文件
function conv(n,    m) {
    if (n in T) return n
    if (("chr" n) in T) return "chr" n
    if (n == "MT" && ("chrM" in T)) return "chrM"
    if (n ~ /^chr/) {
        m = substr(n, 4)
        if (m == "M") m = "MT"
        if (m in T) return m
    }
    return ""
}
FNR == NR { T[$1] = 1; next }
/^#/ { next }
{
    c = conv($col)
    if (c == "") { dropped[$col] = 1; next }
    $col = c
    print
}
END {
    n = 0
    for (k in dropped) { n++; s = s " " k }
    if (n) print "rename_chroms: 目标中不存在、已丢弃的序列 (" n "):" s > "/dev/stderr"
}
