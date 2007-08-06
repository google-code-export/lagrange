def parse_aln(infile):
    ntax, nareas = map(int, infile.readline().split())
    taxon2dist = {}
    for line in [ x for x in infile if x.strip() ]:
        taxon, s = line.split()
        taxon2dist[taxon] = tuple([ bool(int(x)) for x in s ])
    return ntax, nareas, taxon2dist

## infile = file("./cyrtandra/AREACYR2.aln")
## for x in parse_aln(infile):
##     print x
