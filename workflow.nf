params.nboot = 100
params.seed=2000
params.outpath="results"
params.itolconfig= "data/itol_image_config.txt"
params.refseqid = "data/refseq_ids.txt"
params.ncbirenamefile = "data/ncbi_rename.txt"
params.itolkey="none"
params.itolproject="none"

nboot = params.nboot
seed = params.seed
outpath = file(params.outpath)
itolconfig=file(params.itolconfig)
gene2accession=file("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz")
geneinfo=file("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz")
refseqid=file(params.refseqid)
ncbirenamefile=file(params.ncbirenamefile)
itolkey=params.itolkey
itolproject=params.itolproject


// Get the NCBI & HGNC ID Related to the RefSeq protein ID given in 
// https://journals.plos.org/plosbiology/article/file?type=supplementary&id=info:doi/10.1371/journal.pbio.3000954.s008&rev=1
process getNCBIIds{
	publishDir "${outpath}", mode: 'copy'

	input:
	file refseqid
	file gene2accession
	file geneinfo

	output:
	file 'refseq_ids_xref.txt'
	file 'refseq_ids_hgnc.txt' into hgncid

	script:
	"""
	add_hgnc.pl $refseqid $gene2accession $geneinfo | sort -u > refseq_ids_xref.txt
	cut -f 3 refseq_ids_xref.txt | sort -u > refseq_ids_hgnc.txt
	"""
}


// Get OrthoDB sequences corresponding to NCBI ID
process getOrthoDBIds{
	maxForks 3

	input:
	val hgnc from hgncid.splitText( by: 1 ).map{ v -> v.trim() }
	
	output:
	file "ids.txt" into protids,protids2

	script:
	"""
	curl 'https://v100.orthodb.org/search?query=${hgnc}&ncbi=0&level=9443&species=9443&singlecopy=1&universal=0.9'|jq '.data | join(",")' | sed 's/"//g' | sed 's/,/\\n/g' > ids.txt
	sleep 1
	"""
}

protids2.collectFile(name: "all_orthoid.txt").subscribe{file -> file.copyTo(outpath.resolve(file.name))}

// Download sequences and metadata
process downloadSequences{
	maxForks 1
	tag "${id}"

	input:
	val id from protids.collectFile(name: 'allids.txt').splitText( by: 1 ).map{ v -> v.trim() }.unique().filter{ it.length() > 0 }
	
	output:
	tuple val(id), file("sequences.fasta") into sequences

	script:
	"""
	goalign rename -i "https://v100.orthodb.org/fasta?id=${id}"  --regexp '([^\\s]+).*' --replace '\$1' --unaligned > sequences.fasta
	sleep 2
	"""
}


// Download meta data
process getMapTable{
	maxForks 1
	tag "${id}"

	input:
	set val(id), file(seq) from sequences

	output:
	set val(id), file(seq), file("map.txt") into mapfile
	file "gene.txt" into genefile

	script:
	"""
	wget -O align.tab https://v100.orthodb.org/tab?query=${id}
	cut -f 5,6 align.tab > map.txt
	cut -f 1,2 align.tab | tail -n+2 | sort -u > gene.txt
	"""
}

genefile.collectFile(name: 'genes.txt').subscribe{file -> file.copyTo(outpath.resolve(file.name))}

// Renaming sequences
process renameSequences{
	tag "${id}"

	input:
	set val(id), file(sequences), file(mapfile) from mapfile

	output:
	file "renamed.fasta" into renamed

	script:
	"""
	goalign rename -r -m ${mapfile} -i ${sequences} --unaligned | goalign rename --regexp " " --replace "_" --unaligned  > renamed.fasta
	"""
}


// Cleaning sequences
process cleanSequences{
	input:
	file(sequences) from renamed

	output:
	file "cleaned.fasta" into cleaned

	script:
	"""
	goalign replace -s U -n X -i ${sequences} -o cleaned.fasta --unaligned
	"""
}


// Aligning sequences
process alignSequences{
	input:
	file cleaned

	output:
	file "aligned.fasta" into alignment

	script:
	"""
	mafft --quiet ${cleaned} > aligned.fasta
	"""
}


// Concatenating all multiple sequence alignments
process concatSequences {
	input:
	file 'align_fasta' from alignment.collect()

	output:
	file "concat.fasta" into concat

	script:
	"""
	goalign concat -o concat.fasta -i none align_fasta*
	"""
}

// Cleaning the concatenated multiple sequence alignment (BMGE)
process cleanAlign {
	input:
	file align from concat

	output:
	file "cleanalign.fasta" into cleanalign

	script:
	"""
	BMGE -i ${align} -t AA -m BLOSUM62 -w 3 -g 0.2 -h 0.5  -b 5 -of cleanalign.fasta
	"""
}


// Reformatting the alignment in Phylip
process reformatAlignment{
	input:
	file cleanalign

	output:
	file "aligned.phylip" into alignmentphylip

	script:
	"""
	goalign reformat phylip -i ${cleanalign} -o aligned.phylip
	"""
}

// Inferring the species tree from the concatenated 
// multiple sequence alignment
process inferTrueTree{
	publishDir "${outpath}", mode: 'copy'

	input:
	file align from alignmentphylip
	val seed

	output:
	file "tree.nw" into tree, tree2

	script:
	"""
	iqtree -s ${align} -seed ${seed} -m MFP -b 100 -nt ${task.cpus}
	mv *.treefile tree.nw
	"""
}

// Rerooting the tree
process rerootTree{
	publishDir "${outpath}", mode: 'copy'

	input:
	file tree from tree2

	output:
	file "rerooted.nw" into reroottree1, reroottree2

	script:
	"""
	gotree reroot outgroup -i ${tree} -o rerooted.nw Otolemur_garnettii Microcebus_murinus Propithecus_coquereli 
	"""
}

// Downloading the NCBI Taxonomy
process downloadNewickTaxonomy {
	output:
	file "ncbi.nw" into ncbitax

	script:
	"""
	gotree download ncbitax -o ncbi.nw
	"""
}


// Keep the desired species from the NCBI taxonomy
process pruneNCBITax {

	input:
	file tree from tree
	file map from ncbirenamefile
	file ncbi from ncbitax

	output:
	file "ncbi_pruned.nw" into ncbipruned

	script:
	"""
	gotree rename -i ${tree} -m ${map} -r -o tmp 
	gotree prune -i ${ncbi} -c tmp -o ncbi_pruned.nw
	"""
}

// Rename some tips from the NCBI Taxonomy to match 
// orthodb species names
process renameNCBITaxonomy {
	input:
	file ncbi from ncbipruned
	file map from ncbirenamefile

	output:
	file "ncbi_rename.nw" into ncbitaxrename

	script:
	"""
	gotree rename -i $ncbi -o ncbi_rename.nw -m $map
	"""
}

// Reroot the NCBI taxonomy
process rerootNCBITax {

	input:
	file tree from ncbitaxrename

	output:
	file "ncbi_rerooted.nw" into ncbirerooted1, ncbirerooted2

	script:
	"""
	gotree reroot outgroup -i ${tree} -o ncbi_rerooted.nw Otolemur_garnettii Microcebus_murinus Propithecus_coquereli 
	"""
}


// Annotating the inferred tree
process annotateTree{
	publishDir "${outpath}", mode: 'copy'

	input:
	file tree from reroottree1
	file ncbi from ncbirerooted1

	output:
	file "annotated.nw" into annotated

	script:
	"""
	gotree annotate -i ${tree} -c ${ncbi} -o annotated.nw
	"""
}


// Comparing inferred tree and NCBI taxonomy branches
process compareTrees{
	publishDir "${outpath}", mode: 'copy'

	input:
	file tree from reroottree2
	file ncbi from ncbirerooted2

	output:
	file "tree_comparison.txt" into comparison

	script:
	"""
	gotree compare trees -i ${tree} -c ${ncbi} > tree_comparison.txt
	"""
}

// Uploading tree to ITol &
// Downloading the resulting image
process uploadTree{
	publishDir "${outpath}", mode: 'copy'

	input:
	file tree from annotated
	file itolconfig
	val itolkey
	val itolproject

	output:
	file "tree_url.txt" into iTOLurl
	file "tree_image.svg" into iTOLimage

	script:
	"""
	# Upload the tree
	gotree upload itol --name "AnnotatedTree" -i ${tree} --user-id ${itolkey} --project ${itolproject} > tree_url.txt
	# We get the iTOL id
	ID=$(basename $(cat tree_url.txt ))
	# We Download the image with options defined in data/itol_image_config.txt
	gotree download itol -c ${itolconfig} -f svg -o tree_image.svg -i $ID
	"""
}
