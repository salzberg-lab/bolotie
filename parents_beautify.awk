# awk -F',' -f .parents_beautify.awk data/gisaid_hcov-19_2020_09_02_16.fasta parents.clean_table.csv > output

function print_name_in_cols(str) {
    # country = substr(str, 1, index(str, "_") - 1) # cut off after first _
    country = countries[gisaid[1]]
    if (country == "SaudiArabia")
        country = "Saudi Arabia"
    match(str, /(EPI_ISL_[0-9]+)/, gisaid)
    # YYYY-MM-DD
    pos = match(str, /([0-9]{4}-[0-9]{2}-[0-9]{2}$)/, date)
    # YYYY-MM
    if (pos == 0)
        pos = match(str, /([0-9]{4}-[0-9]{2}$)/, date)
    # YYYY
    if (pos == 0)
        pos = match(str, /([0-9]{4}$)/, date)
    printf "%s\t%s\t%s\t", gisaid[1], country, date[1]
}
BEGIN {
    OFS="\t"
    # headline
    print "Recomb GISAID", "Country", "Date",
          "P1 GISAID", "Country", "Date", "Clade", "Segment",
          "P2 GISAID", "Country", "Date", "Clade", "Segment"
}
FNR==NR {
    if ($0 ~ /^>/)
    {
        suffix = substr($0, 10) # removes >hCov-19/
        country = substr(suffix, 1, index(suffix, "/") - 1) # Country/...
        match($0, /(EPI_ISL_[0-9]+)/, gisaid)
        countries[gisaid[1]] = country
    }
    next
}
{
    print_name_in_cols($1)
    printf "%s\t%s\t", $2, $4 # P1 clade and segment
    print_name_in_cols($3)
    printf "%s\t%s\t", $5, $7 # P2 clade and segment
    print_name_in_cols($6)
    print ""
}
